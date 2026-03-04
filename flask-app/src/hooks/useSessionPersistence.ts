/**
 * Session persistence hook for auto-saving experiment state to MariaDB.
 * This provides automatic session save/restore functionality.
 */
import { useState, useEffect, useRef, useCallback, MutableRefObject } from 'react';
import {
  TreeNode,
  Edge,
  MetricHistoryItem,
  VisibleMetrics,
  SidebarMessage,
  VisibleSources,
} from '../types';

// Auto-save interval in milliseconds
const AUTO_SAVE_INTERVAL = 10000; // 10 seconds

// Session expiry in hours (ignore sessions older than this)
const SESSION_EXPIRY_HOURS = 24;

// Session state interface matching the backend API
export interface SessionState {
  // Project/Experiment identification
  projectId?: string;
  projectName?: string;
  experimentId?: string;
  experimentName?: string;

  // Core experiment state
  smiles?: string;
  problemType?: string;
  systemPrompt?: string;
  problemPrompt?: string;
  promptsModified?: boolean;
  autoZoom?: boolean;
  treeNodes?: TreeNode[];
  edges?: Edge[];
  offset?: { x: number; y: number };
  zoom?: number;
  metricsHistory?: MetricHistoryItem[];
  visibleMetrics?: VisibleMetrics;
  isComputing?: boolean;
  serverSessionId?: string | null;

  // Property optimization
  propertyType?: string;
  customPropertyName?: string;
  customPropertyDesc?: string;
  customPropertyAscending?: boolean;

  // Sidebar state
  sidebarMessages?: SidebarMessage[];
  sidebarSourceFilterOpen?: boolean;
  sidebarVisibleSources?: VisibleSources;
}

interface SessionResponse {
  sessionId: string;
  name: string;
  createdAt: string;
  lastModified: string;
  isRunning?: boolean;
  state: SessionState;
}

export interface SessionPersistenceState {
  dbSessionId: string | null;
  sessionLoaded: boolean;
  sessionRestored: boolean;
  lastSaved: Date | null;
  isSaving: boolean;
  saveError: string | null;
  serverSessionId: string | null;
  sessionWasComputing: boolean;
  resumeAttempted: boolean;
}

export interface SessionPersistenceActions {
  saveSession: (
    state: SessionState,
    force?: boolean,
    options?: { checkpoint?: boolean; name?: string }
  ) => Promise<void>;
  loadSession: () => Promise<SessionState | null>;
  clearSession: () => Promise<void>;
  setServerSessionId: (id: string | null) => void;
  setSessionWasComputing: (computing: boolean) => void;
  setResumeAttempted: (attempted: boolean) => void;
  getPendingResume: () => { sessionId: string | null; smiles: string; problemType: string } | null;
  setPendingResume: (
    data: { sessionId: string | null; smiles: string; problemType: string } | null
  ) => void;
  handleWebSocketMessage: (data: any) => boolean;
  sendSessionSaveNow: (
    state: SessionState,
    options?: { checkpoint?: boolean; name?: string }
  ) => void;
}

export const useSessionPersistence = (
  wsRef?: MutableRefObject<WebSocket | null>
): SessionPersistenceState & SessionPersistenceActions => {
  const [dbSessionId, setDbSessionId] = useState<string | null>(null);
  const [sessionLoaded, setSessionLoaded] = useState(false);
  const [sessionRestored, setSessionRestored] = useState(false);
  const [lastSaved, setLastSaved] = useState<Date | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const [saveError, setSaveError] = useState<string | null>(null);
  const [serverSessionId, setServerSessionId] = useState<string | null>(null);
  const [sessionWasComputing, setSessionWasComputing] = useState(false);
  const [resumeAttempted, setResumeAttempted] = useState(false);

  const isSavingRef = useRef(false);
  const dbSessionIdRef = useRef<string | null>(null);
  const pendingResumeRef = useRef<{
    sessionId: string | null;
    smiles: string;
    problemType: string;
  } | null>(null);
  const saveQueueRef = useRef<
    Array<{
      state: SessionState;
      checkpoint: boolean;
      name?: string;
      resolve: () => void;
    }>
  >([]);
  const processingQueueRef = useRef(false);
  const requestCounterRef = useRef(0);
  const pendingRequestsRef = useRef<
    Map<
      string,
      {
        resolve: (data: any) => void;
        reject: (error: Error) => void;
        timeoutId: number;
      }
    >
  >(new Map());

  useEffect(() => {
    dbSessionIdRef.current = dbSessionId;
  }, [dbSessionId]);

  const waitForOpenWebSocket = useCallback(
    async (timeoutMs = 10000): Promise<WebSocket> => {
      const startedAt = Date.now();

      while (Date.now() - startedAt < timeoutMs) {
        const socket = wsRef?.current;
        if (socket && socket.readyState === WebSocket.OPEN) {
          return socket;
        }
        await new Promise((resolve) => setTimeout(resolve, 100));
      }

      throw new Error('WebSocket not connected');
    },
    [wsRef]
  );

  const sendWsRequest = useCallback(
    async (action: string, payload: Record<string, any>, timeoutMs = 10000): Promise<any> => {
      const socket = await waitForOpenWebSocket(timeoutMs);
      const requestId = `${action}-${Date.now()}-${requestCounterRef.current++}`;

      return new Promise((resolve, reject) => {
        const timeoutId = window.setTimeout(() => {
          pendingRequestsRef.current.delete(requestId);
          reject(new Error(`${action} timed out`));
        }, timeoutMs);

        pendingRequestsRef.current.set(requestId, { resolve, reject, timeoutId });

        socket.send(
          JSON.stringify({
            action,
            requestId,
            ...payload,
          })
        );
      });
    },
    [waitForOpenWebSocket]
  );

  const sendSessionSaveNow = useCallback(
    (state: SessionState, options?: { checkpoint?: boolean; name?: string }): void => {
      const socket = wsRef?.current;
      if (!socket || socket.readyState !== WebSocket.OPEN) return;

      const checkpoint = options?.checkpoint ?? false;
      const effectiveSessionId = checkpoint ? null : state.experimentId || dbSessionIdRef.current;

      socket.send(
        JSON.stringify({
          action: 'session_save',
          sessionId: effectiveSessionId,
          name: options?.name,
          state,
        })
      );
    },
    [wsRef]
  );

  const processSaveQueue = useCallback(async (): Promise<void> => {
    if (processingQueueRef.current) return;
    const next = saveQueueRef.current.shift();
    if (!next) return;

    processingQueueRef.current = true;
    isSavingRef.current = true;
    setIsSaving(true);
    setSaveError(null);

    try {
      // Prefer the sidebar's experimentId so saves always update the
      // experiment the user is working in.  Fall back to the session
      // persistence tracking ID for edge cases.
      const effectiveSessionId = next.checkpoint
        ? null
        : next.state.experimentId || dbSessionIdRef.current;
      const data = await sendWsRequest('session_save', {
        sessionId: effectiveSessionId,
        name: next.name,
        state: next.state,
      });

      if (data?.session) {
        const sessionData: SessionResponse = data.session;
        setDbSessionId(sessionData.sessionId);
        setLastSaved(new Date(sessionData.lastModified));
        console.log('Session saved over WebSocket:', sessionData.sessionId);
        next.resolve();
      } else {
        setSaveError('Save failed: malformed response');
        next.resolve();
      }
    } catch (error) {
      console.error('Failed to save session over WebSocket:', error);
      setSaveError('Save failed: network error');
      next.resolve();
    } finally {
      isSavingRef.current = false;
      setIsSaving(false);
      processingQueueRef.current = false;
      if (saveQueueRef.current.length > 0) {
        void processSaveQueue();
      }
    }
  }, [sendWsRequest]);

  const saveSession = useCallback(
    async (
      state: SessionState,
      force = false,
      options?: { checkpoint?: boolean; name?: string }
    ): Promise<void> => {
      // Only save if there's meaningful data or force save
      if (!force && (!state.treeNodes || state.treeNodes.length === 0) && !state.smiles) return;

      const checkpoint = options?.checkpoint ?? false;
      const name = options?.name;

      return new Promise<void>((resolve) => {
        saveQueueRef.current.push({ state, checkpoint, name, resolve });

        if (!processingQueueRef.current) {
          void processSaveQueue();
        }
      });
    },
    [processSaveQueue]
  );

  const loadSession = useCallback(async (): Promise<SessionState | null> => {
    try {
      const response = await sendWsRequest('session_get_latest', {});
      const data: SessionResponse | null = response?.session ?? null;

      if (!data) {
        setSessionLoaded(true);
        return null;
      }

      const state = data.state;

      // Check if session is recent (within SESSION_EXPIRY_HOURS)
      const sessionTime = new Date(data.lastModified);
      const now = new Date();
      const hoursDiff = (now.getTime() - sessionTime.getTime()) / (1000 * 60 * 60);

      if (hoursDiff > SESSION_EXPIRY_HOURS) {
        console.log('Session too old, ignoring');
        setSessionLoaded(true);
        return null;
      }

      // Store the database session ID
      setDbSessionId(data.sessionId);

      // Check if computation was in progress - either by isComputing flag or by detecting incomplete edges
      const hasIncompleteEdges = state.edges && state.edges.some((e) => e.status === 'computing');
      const wasComputing = state.isComputing || hasIncompleteEdges;

      if (wasComputing && state.smiles) {
        console.log('Detected incomplete computation, will attempt to resume');
        setSessionWasComputing(true);
        if (state.serverSessionId) {
          setServerSessionId(state.serverSessionId);
        }
        pendingResumeRef.current = {
          sessionId: state.serverSessionId || null,
          smiles: state.smiles,
          problemType: state.problemType || 'retrosynthesis',
        };
      }

      setSessionLoaded(true);
      setSessionRestored(true);
      setLastSaved(sessionTime);

      console.log(
        'Session restored from database:',
        data.sessionId,
        state.serverSessionId ? `(server session: ${state.serverSessionId})` : ''
      );
      return state;
    } catch (error) {
      console.error('Failed to load session over WebSocket:', error);
      setSessionLoaded(true);
      return null;
    }
  }, [sendWsRequest]);

  const clearSession = useCallback(async (): Promise<void> => {
    if (dbSessionId) {
      try {
        await sendWsRequest('session_delete', { sessionId: dbSessionId });
        console.log('Session deleted over WebSocket:', dbSessionId);
      } catch (error) {
        console.error('Failed to delete session over WebSocket:', error);
      }
    }
    setDbSessionId(null);
    setSessionRestored(false);
    setLastSaved(null);
    setServerSessionId(null);
    setSessionWasComputing(false);
    setResumeAttempted(false);
    pendingResumeRef.current = null;
    console.log('Session cleared');
  }, [dbSessionId, sendWsRequest]);

  const handleWebSocketMessage = useCallback((data: any): boolean => {
    const responseTypes = new Set([
      'session_save_response',
      'session_get_latest_response',
      'session_get_response',
      'session_list_response',
      'session_delete_response',
      'session_error',
    ]);

    if (!data?.type || !responseTypes.has(data.type)) {
      return false;
    }

    const requestId = data.requestId;
    if (!requestId) {
      return true;
    }

    const pending = pendingRequestsRef.current.get(requestId);
    if (!pending) {
      return true;
    }

    window.clearTimeout(pending.timeoutId);
    pendingRequestsRef.current.delete(requestId);

    if (data.type === 'session_error') {
      pending.reject(new Error(data.error || 'Session operation failed'));
      return true;
    }

    pending.resolve(data);
    return true;
  }, []);

  useEffect(() => {
    return () => {
      pendingRequestsRef.current.forEach((pending) => {
        window.clearTimeout(pending.timeoutId);
        pending.reject(new Error('Session request cancelled'));
      });
      pendingRequestsRef.current.clear();
    };
  }, []);

  const getPendingResume = useCallback(() => pendingResumeRef.current, []);

  const setPendingResume = useCallback(
    (data: { sessionId: string | null; smiles: string; problemType: string } | null) => {
      pendingResumeRef.current = data;
    },
    []
  );

  return {
    // State
    dbSessionId,
    sessionLoaded,
    sessionRestored,
    lastSaved,
    isSaving,
    saveError,
    serverSessionId,
    sessionWasComputing,
    resumeAttempted,
    // Actions
    saveSession,
    loadSession,
    clearSession,
    setServerSessionId,
    setSessionWasComputing,
    setResumeAttempted,
    getPendingResume,
    setPendingResume,
    handleWebSocketMessage,
    sendSessionSaveNow,
  };
};

/**
 * Hook to set up checkpoint save on page unload.
 * Checkpoints are now triggered on-demand (e.g., when new molecules are generated)
 * rather than at fixed intervals.
 */
export const useCheckpointOnUnload = (
  sessionPersistence: SessionPersistenceState & SessionPersistenceActions,
  getState: () => SessionState
): void => {
  const { saveSession, sendSessionSaveNow } = sessionPersistence;
  const getStateRef = useRef(getState);

  // Keep the getState ref updated
  useEffect(() => {
    getStateRef.current = getState;
  }, [getState]);

  // Save session before page unload (best-effort websocket send)
  useEffect(() => {
    const handleBeforeUnload = () => {
      const state = getStateRef.current();
      if ((!state.treeNodes || state.treeNodes.length === 0) && !state.smiles) return;

      sendSessionSaveNow(state);
    };

    window.addEventListener('beforeunload', handleBeforeUnload);
    return () => {
      window.removeEventListener('beforeunload', handleBeforeUnload);
      // Save on unmount
      saveSession(getStateRef.current(), true);
    };
  }, [saveSession, sendSessionSaveNow]);
};

// Keep the old name for backward compatibility but it now just sets up unload handler
export const useAutoSave = useCheckpointOnUnload;
