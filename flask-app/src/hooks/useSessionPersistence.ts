/**
 * Session persistence hook for auto-saving experiment state to MariaDB.
 * This provides automatic session save/restore functionality.
 */
import { useState, useEffect, useRef, useCallback } from 'react';
import { TreeNode, Edge, MetricHistoryItem, VisibleMetrics, SidebarMessage, VisibleSources } from '../types';
import { HTTP_SERVER } from '../config';

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
  saveSession: (state: SessionState, force?: boolean, options?: { checkpoint?: boolean; name?: string }) => Promise<void>;
  loadSession: () => Promise<SessionState | null>;
  clearSession: () => Promise<void>;
  setServerSessionId: (id: string | null) => void;
  setSessionWasComputing: (computing: boolean) => void;
  setResumeAttempted: (attempted: boolean) => void;
  getPendingResume: () => { sessionId: string | null; smiles: string; problemType: string } | null;
  setPendingResume: (data: { sessionId: string | null; smiles: string; problemType: string } | null) => void;
}

export const useSessionPersistence = (): SessionPersistenceState & SessionPersistenceActions => {
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
  const pendingResumeRef = useRef<{ sessionId: string | null; smiles: string; problemType: string } | null>(null);
  const saveQueueRef = useRef<Array<{
    state: SessionState;
    checkpoint: boolean;
    name?: string;
    resolve: () => void;
  }>>([]);
  const processingQueueRef = useRef(false);

  useEffect(() => {
    dbSessionIdRef.current = dbSessionId;
  }, [dbSessionId]);

  const processSaveQueue = useCallback(async (): Promise<void> => {
    if (processingQueueRef.current) return;
    const next = saveQueueRef.current.shift();
    if (!next) return;

    processingQueueRef.current = true;
    isSavingRef.current = true;
    setIsSaving(true);
    setSaveError(null);

    try {
      const response = await fetch(`${HTTP_SERVER}/api/sessions/save`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          sessionId: next.checkpoint ? null : dbSessionIdRef.current,
          name: next.name,
          state: next.state,
        }),
      });

      if (response.ok) {
        const data: SessionResponse = await response.json();
        setDbSessionId(data.sessionId);
        setLastSaved(new Date(data.lastModified));
        console.log('Session saved to database:', data.sessionId);
        next.resolve();
      } else {
        console.error('Failed to save session to database:', response.status);
        setSaveError(`Save failed: ${response.status}`);
        next.resolve();
      }
    } catch (error) {
      console.error('Failed to save session to database:', error);
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
  }, []);

  const saveSession = useCallback(async (
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
  }, [processSaveQueue]);

  const loadSession = useCallback(async (): Promise<SessionState | null> => {
    try {
      const response = await fetch(`${HTTP_SERVER}/api/sessions/latest`);
      
      if (!response.ok) {
        if (response.status === 404) {
          console.log('No previous session found');
        }
        setSessionLoaded(true);
        return null;
      }
      
      const data: SessionResponse | null = await response.json();
      
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
      const hasIncompleteEdges = state.edges && state.edges.some(e => 
        e.status === 'computing'
      );
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
          problemType: state.problemType || 'retrosynthesis'
        };
      }
      
      setSessionLoaded(true);
      setSessionRestored(true);
      setLastSaved(sessionTime);
      
      console.log('Session restored from database:', data.sessionId, state.serverSessionId ? `(server session: ${state.serverSessionId})` : '');
      return state;
    } catch (error) {
      console.error('Failed to load session from database:', error);
      setSessionLoaded(true);
      return null;
    }
  }, []);

  const clearSession = useCallback(async (): Promise<void> => {
    if (dbSessionId) {
      try {
        await fetch(`${HTTP_SERVER}/api/sessions/${dbSessionId}`, {
          method: 'DELETE',
        });
        console.log('Session deleted from database:', dbSessionId);
      } catch (error) {
        console.error('Failed to delete session from database:', error);
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
  }, [dbSessionId]);

  const getPendingResume = useCallback(() => pendingResumeRef.current, []);
  
  const setPendingResume = useCallback((data: { sessionId: string | null; smiles: string; problemType: string } | null) => {
    pendingResumeRef.current = data;
  }, []);

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
  const { dbSessionId, saveSession } = sessionPersistence;
  const getStateRef = useRef(getState);
  
  // Keep the getState ref updated
  useEffect(() => {
    getStateRef.current = getState;
  }, [getState]);

  // Save session before page unload (using sendBeacon for reliability)
  useEffect(() => {
    const handleBeforeUnload = () => {
      const state = getStateRef.current();
      if ((!state.treeNodes || state.treeNodes.length === 0) && !state.smiles) return;
      
      const payload = JSON.stringify({
        sessionId: dbSessionId,
        state: state,
      });
      
      // sendBeacon is more reliable than fetch for unload events
      navigator.sendBeacon(`${HTTP_SERVER}/api/sessions/save`, new Blob([payload], { type: 'application/json' }));
    };
    
    window.addEventListener('beforeunload', handleBeforeUnload);
    return () => {
      window.removeEventListener('beforeunload', handleBeforeUnload);
      // Save on unmount
      saveSession(getStateRef.current(), true);
    };
  }, [dbSessionId, saveSession]);
};

// Keep the old name for backward compatibility but it now just sets up unload handler
export const useAutoSave = useCheckpointOnUnload;
