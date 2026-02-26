/* eslint-disable react-hooks/exhaustive-deps */
import React, { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import {
  Loader2,
  FlaskConical,
  TestTubeDiagonal,
  Network,
  Play,
  RotateCcw,
  X,
  Send,
  RefreshCw,
  Sparkles,
  MessageCircleQuestion,
  StepForward,
  MessageSquareShare,
  Brain,
  PanelRightOpen,
  Sliders,
  Wrench,
  Settings,
  Bug,
  CheckCircle,
  Minus,
} from 'lucide-react';
import 'recharts';
import 'react-markdown';
import 'remark-gfm';
import 'react-syntax-highlighter';
import 'react-syntax-highlighter/dist/esm/styles/prism';

import { WS_SERVER, HTTP_SERVER, VERSION } from './config';
import { DEFAULT_CUSTOM_SYSTEM_PROMPT, PROPERTY_NAMES } from './constants';
import {
  TreeNode,
  Edge,
  ContextMenuState,
  Tool,
  WebSocketMessageToServer,
  WebSocketMessage,
  SelectableTool,
  Experiment,
  FlaskOrchestratorSettings as OrchestratorSettings,
  OptimizationCustomization,
  ReactionAlternative,
} from './types';

import { loadRDKit } from './components/molecule';
import { MoleculeGraph, useGraphState } from './components/graph';
import {
  ProjectSidebar,
  useProjectSidebar,
  useProjectManagement,
} from './components/project_sidebar';
import {
  SidebarMessage,
  SettingsButton,
  ReasoningSidebar,
  useSidebarState,
  MarkdownText,
  BACKEND_OPTIONS,
} from 'lc-conductor';
import { CombinedCustomizationModal } from './components/combined_customization_modal';
import { Modal } from './components/modal';

import { findAllDescendants, isRootNode, relayoutTree } from './tree_utils';
import { copyToClipboard } from './utils';

import './animations.css';
import { MetricsDashboard, useMetricsDashboardState } from './components/metrics';
import { useProjectData } from './hooks/useProjectData';
import { ReactionAlternativesSidebar } from './components/reaction_alternatives';
import { useSessionPersistence, useAutoSave, SessionState } from './hooks/useSessionPersistence';


const ChemistryTool: React.FC = () => {
  const [smiles, setSmiles] = useState<string>('');
  const [problemType, setProblemType] = useState<string>('retrosynthesis');
  const [propertyType, setPropertyType] = useState<string>('density');
  const [systemPrompt, setSystemPrompt] = useState<string>(DEFAULT_CUSTOM_SYSTEM_PROMPT);
  const [problemPrompt, setProblemPrompt] = useState<string>('');
  const [editPromptsModal, setEditPromptsModal] = useState<boolean>(false);
  const [editPropertyModal, setEditPropertyModal] = useState<boolean>(false);
  const [showCustomizationModal, setShowCustomizationModal] = useState<boolean>(false);
  const [customPropertyName, setCustomPropertyName] = useState<string>('');
  const [customPropertyDesc, setCustomPropertyDesc] = useState<string>('');
  const [customPropertyAscending, setCustomPropertyAscending] = useState<boolean>(true);
  const [isComputing, setIsComputing] = useState<boolean>(false);
  const [autoZoom, setAutoZoom] = useState<boolean>(true);
  const [treeNodes, setTreeNodes] = useState<TreeNode[]>([]);
  const [edges, setEdges] = useState<Edge[]>([]);

  // Track which experiments have running computations: experimentId → serverSessionId
  const computingExperimentsRef = useRef<Map<string, string>>(new Map());
  const [contextMenu, setContextMenu] = useState<ContextMenuState>({node: null, isReaction: false, x: 0, y: 0});
  const [customQueryModal, setCustomQueryModal] = useState<TreeNode | null>(null);
  const [customQueryText, setCustomQueryText] = useState<string>('');
  const [customQueryType, setCustomQueryType] = useState<string | null>(null);
  const [wsConnected, setWsConnected] = useState<boolean>(false);
  const [saveDropdownOpen, setSaveDropdownOpen] = useState<boolean>(false);
  const [loadDropdownOpen, setLoadDropdownOpen] = useState<boolean>(false);
  const [checkpointDropdownOpen, setCheckpointDropdownOpen] = useState<boolean>(false);
  const [autoCheckpointEnabled, setAutoCheckpointEnabled] = useState<boolean>(true);
  const [wsError, setWsError] = useState<string>('');
  const [wsReconnecting, setWsReconnecting] = useState<boolean>(false);
  const rdkitModule = loadRDKit();
  const [sidebarOpen, setSidebarOpen] = useState<boolean>(false);
  const [copiedField, setCopiedField] = useState<string | null>(null);
  const [availableTools, setAvailableTools] = useState<Tool[]>([]);
  const [wsTooltipPinned, setWsTooltipPinned] = useState<boolean>(false);
  const [username, setUsername] = useState<string>('<LOCAL USER>');

  const wsRef = useRef<WebSocket | null>(null);

  // Customization state
  const [customization, setCustomization] = useState<OptimizationCustomization>({
    enableConstraints: false,
    molecularSimilarity: 0.7,
    diversityPenalty: 0.0,
    explorationRate: 0.5,
    additionalConstraints: [],
  });

  // AI debugging
  const [debugMode, setDebugMode] = useState<boolean>(false);
  const [promptBreakpoint, setPromptBreakpoint] = useState<{
    prompt: string;
    metadata?: any;
  } | null>(null);
  const [editedPrompt, setEditedPrompt] = useState<string>('');
  const [debugModalMinimized, setDebugModalMinimized] = useState<boolean>(false);

  // Function to refresh tools list from backend
  const refreshToolsList = useCallback(() => {
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      console.log('Refreshing tools list from backend');
      wsRef.current.send(JSON.stringify({ action: 'list-tools' }));
    }
  }, []);

  const getContextRef = useRef<() => Experiment>(() => {
    throw new Error("getContext called before initialization");
  });

  const graphState = useGraphState();
  const sidebarState = useSidebarState();
  const metricsDashboardState = useMetricsDashboardState();
  const projectSidebar = useProjectSidebar();
  const projectData = useProjectData();
  const projectManagement = useProjectManagement(projectData);
  const sessionPersistence = useSessionPersistence();

  const treeNodesRef = useRef(treeNodes);
  const edgesRef = useRef(edges);
  const sidebarStateRef = useRef(sidebarState);
  const saveCheckpointRef = useRef<() => void>(() => {});
  const saveStateToExperimentRef = useRef<() => boolean>(() => false);
  const restoredSessionRef = useRef<SessionState | null>(null);
  const sessionLoadCompleteRef = useRef(false);
  const isExperimentCompleteRef = useRef(false);
  const autoCheckpointEnabledRef = useRef(true);

  const [selectedTools, setSelectedTools] = useState<number[]>([]);
  const [availableToolsMap, setAvailableToolsMap] = useState<SelectableTool[]>([]);

  // Reaction alternatives sidebar state
  const [reactionSidebarOpen, setReactionSidebarOpen] = useState<boolean>(false);
  const [selectedReactionNode, setSelectedReactionNode] = useState<TreeNode | null>(null);
  const [isComputingTemplates, setIsComputingTemplates] = useState<boolean>(false);

  // Keep selectedReactionNode in sync with treeNodes updates
  useEffect(() => {
    if (selectedReactionNode) {
      const updatedNode = treeNodes.find(n => n.id === selectedReactionNode.id);
      if (updatedNode && updatedNode !== selectedReactionNode) {
        setSelectedReactionNode(updatedNode);
      }
    }
  }, [treeNodes, selectedReactionNode]);

  // Auto-select all tools when they first become available
  useEffect(() => {
    if (availableToolsMap.length > 0 && selectedTools.length === 0) {
      const allToolIds = availableToolsMap.map(tool => tool.id);
      setSelectedTools(allToolIds);
      console.log('Auto-selected all tools:', allToolIds);
    }
  }, [availableToolsMap, selectedTools.length]);

  // Update refs whenever state changes
  useEffect(() => {
    treeNodesRef.current = treeNodes;
  }, [treeNodes]);

  useEffect(() => {
    edgesRef.current = edges;
  }, [edges]);

  useEffect(() => {
    sidebarStateRef.current = sidebarState;
  }, [sidebarState]);

  useEffect(() => {
    autoCheckpointEnabledRef.current = autoCheckpointEnabled;
  }, [autoCheckpointEnabled]);


  // Load initial settings from localStorage
  const getInitialSettings = () => {
    const saved = localStorage.getItem('orchestratorSettings');
    if (saved) {
      try {
        return JSON.parse(saved);
      } catch (e) {
        console.error('Error parsing settings:', e);
      }
    }
    return {
      backend: 'vllm',
      customUrl: 'http://localhost:8000/v1',
      model: 'gpt-oss',
      apiKey: '',
      backendLabel: 'vLLM'
    };
  };
  const [orchestratorSettings, setOrchestratorSettings] = useState<OrchestratorSettings>(getInitialSettings());

  // Add this helper function near the top of the ChemistryTool component
  const getDisplayUrl = (): string => {
    if (orchestratorSettings.useCustomUrl && orchestratorSettings.customUrl) {
      return orchestratorSettings.customUrl;
    }
    const backendOption = BACKEND_OPTIONS.find(opt => opt.value === orchestratorSettings.backend);
    return backendOption?.defaultUrl || 'Not configured';
  };

  // Callback function to send selected tools to backend
  const handleToolSelectionConfirm = async (
    selectedIds: number[],
    selectedItemsData: SelectableTool[]
  ): Promise<void> => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    console.log(`Set Task Tool Selection`);

    if (wsRef.current && wsRef.current.readyState == WebSocket.OPEN) {
      const message: WebSocketMessageToServer = {
        action: "select-tools-for-task",
          enabledTools: {
            selectedIds: selectedIds,
            selectedTools: selectedItemsData,
          }
      };

      wsRef.current.send(JSON.stringify(message));
      console.log('Sending data:', JSON.stringify(message));
    }

    // Optional: Add any additional processing or API calls here
    // await fetch(HTTP_SERVER + '/api/save-selection', { method: 'POST', body: JSON.stringify(payload) });
  };

  // Callback function to handle molecule name preference changes
  const handleMoleculeNameSave = async (moleculeName: string): Promise<void> => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    console.log(`Updated Molecule Name Preference: ${moleculeName}`);

    // Update local orchestrator settings with new molecule name
    const updatedSettings = {
      ...orchestratorSettings,
      runSettings: { moleculeName: moleculeName },
    };
    setOrchestratorSettings(updatedSettings);
    localStorage.setItem('orchestratorSettings', JSON.stringify(updatedSettings));

    // Note: The molecule name is sent to backend as part of runSettings
    // in each compute action, so no separate websocket message needed here.
  };

  // Callback function to send updated settings to backend
  const handleSettingsUpdateConfirm = async (
    settings: OrchestratorSettings,
  ): Promise<void> => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    console.log(`Updated Settings Saved`);

    // Update local state immediately
    setOrchestratorSettings(settings);
    localStorage.setItem('orchestratorSettings', JSON.stringify(settings));

    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      const message = {
        action: 'ui-update-orchestrator-settings',
        backend: settings.backend,
        customUrl: settings.customUrl,
        model: settings.model,
        apiKey: settings.apiKey,
        moleculeName: settings.moleculeName
      };
      wsRef.current.send(JSON.stringify(message));

      // Refresh tools list after updating settings
      console.log('Refreshing tools list after settings update');
      refreshToolsList();
     }
  };

  useEffect(() => {
    const handleClickOutside = (): void => {
      setContextMenu({node: null, isReaction: false, x:0, y:0});
      setSaveDropdownOpen(false);
      setLoadDropdownOpen(false);
      setCheckpointDropdownOpen(false);
      sidebarState.setSourceFilterOpen(false);
      setWsTooltipPinned(false);
      setCopiedField(null);
    };
    if (contextMenu || saveDropdownOpen || loadDropdownOpen || checkpointDropdownOpen || sidebarState.sourceFilterOpen || wsTooltipPinned) {
      window.addEventListener('mousedown', handleClickOutside);
      return () => window.removeEventListener('mousedown', handleClickOutside);
    }
  }, [contextMenu, saveDropdownOpen, loadDropdownOpen, checkpointDropdownOpen, sidebarState, wsTooltipPinned, projectSidebar]);


  // State management
  const getContext = (): Experiment => {
    return getContextRef.current();
  }

  const loadContextFromExperiment = (projectId: string, experimentId: string | null): void => {
    console.log('Loading context:', { projectId, experimentId, projectCount: projectData.projectsRef.current.length });
    const project = projectData.projectsRef.current.find(p => p.id === projectId);
    if (project) {
      const experiment = project.experiments.find(e => e.id === experimentId);
      if (experiment) {
        loadContext(experiment);
      } else {
        console.warn('Experiment not found in project, triggering refresh:', experimentId);
        // Force a DB refresh then retry once
        projectData.refreshProjects().then(() => {
          const p2 = projectData.projectsRef.current.find(p => p.id === projectId);
          const e2 = p2?.experiments.find(e => e.id === experimentId);
          if (e2) loadContext(e2);
          else console.error('Experiment still not found after refresh:', experimentId);
        });
      }
    } else {
      console.warn('Project not found in projectsRef, triggering refresh:', projectId);
      // Force a DB refresh then retry once
      projectData.refreshProjects().then(() => {
        const p2 = projectData.projectsRef.current.find(p => p.id === projectId);
        if (p2) {
          const e2 = p2.experiments.find(e => e.id === experimentId);
          if (e2) loadContext(e2);
          else console.error('Experiment not found after refresh:', experimentId);
        } else {
          console.error('Project still not found after refresh:', projectId);
        }
      });
    }

    // Multi-experiment: check if the experiment we're switching to has a
    // running session.  If so, resume it; otherwise, detach the current
    // server session so background computation continues headless.
    if (experimentId) {
      const sessionId = computingExperimentsRef.current.get(experimentId);
      if (sessionId && wsRef.current?.readyState === WebSocket.OPEN) {
        // If the experiment already completed headless (the 2-second DB
        // poll picked up isRunning=false), skip the resume and load
        // directly from the in-memory cache which already has the
        // sidebar_state from DB.
        const cachedExp = project?.experiments.find(e => e.id === experimentId);
        if (cachedExp && !cachedExp.isRunning) {
          console.log('Experiment completed headless, loading from cache:', experimentId);
          computingExperimentsRef.current.delete(experimentId);
          loadContext(cachedExp);
        } else {
          // Switching to a running experiment – resume the session stream
          setIsComputing(true);
          wsRef.current.send(JSON.stringify({
            action: 'resume_session',
            sessionId,
          }));
          console.log('Resuming session for experiment:', experimentId, sessionId);
          return;
        }
      }
    }

    // Switching to a non-computing experiment – detach the current
    // session on the server so it continues headless without streaming
    // here.  This is a no-op if no session is active.
    if (wsRef.current?.readyState === WebSocket.OPEN) {
      wsRef.current.send(JSON.stringify({ action: 'detach' }));
    }
  }

  const loadStateFromCurrentExperiment = (): void => {
    const { projectId, experimentId } = projectSidebar.selectionRef.current;
    console.log('loadStateFromCurrentExperiment:', { projectId, experimentId });
    if (projectId && experimentId) {
      loadContextFromExperiment(projectId, experimentId);
    }
  };

  const loadContext = (data: Experiment): void => {
    // Always reset to defaults first, then apply experiment data.
    // This ensures switching to a fresh/empty experiment clears old state.
    setSmiles(data.smiles ?? '');
    setProblemType(data.problemType ?? 'retrosynthesis');
    setSystemPrompt(data.systemPrompt ?? DEFAULT_CUSTOM_SYSTEM_PROMPT);
    setProblemPrompt(data.problemPrompt ?? '');
    setPropertyType(data.propertyType ?? 'density');
    setCustomPropertyName(data.customPropertyName ?? '');
    setCustomPropertyDesc(data.customPropertyDesc ?? '');
    setCustomPropertyAscending(data.customPropertyAscending ?? true);
    setCustomization(data.customization ?? {
      enableConstraints: false,
      molecularSimilarity: 0.7,
      diversityPenalty: 0.0,
      explorationRate: 0.5,
      additionalConstraints: [],
    });
    setTreeNodes(data.treeNodes ?? []);
    setEdges(data.edges ?? []);
    metricsDashboardState.setMetricsHistory(data.metricsHistory ?? []);
    metricsDashboardState.setVisibleMetrics(data.visibleMetrics ?? { cost: false, bandgap: false, sascore: false, yield: true, density: false });
    if (data.graphState) {
      graphState.setZoom(data.graphState.zoom);
      graphState.setOffset(data.graphState.offset);
    } else {
      graphState.setZoom(1);
      graphState.setOffset({ x: 50, y: 50 });
    }
    setAutoZoom(data.autoZoom ?? true);
    if (data.sidebarState) {
      sidebarState.setMessages(data.sidebarState.messages);
      sidebarState.setVisibleSources(data.sidebarState.visibleSources);
    } else {
      sidebarState.setMessages([]);
    }
    if (data.experimentContext) {
      sendMessageToServer('load-context', {experimentContext: data.experimentContext});
    }
  }

  const saveStateToExperiment = useCallback((): boolean => {
    // Use the ref directly to always get the latest selection
    const projectId = projectSidebar.selectionRef.current.projectId;
    const experimentId = projectSidebar.selectionRef.current.experimentId;
    console.log("Saving experiments", projectId, experimentId);
    if (projectId && experimentId) {
      try {
        const ctx = getContext();

        // Guard: don't save empty local state over non-empty DB data.
        // This prevents a session-restore + experiment-switch cycle from
        // overwriting data that save_session_to_db persisted while the
        // browser was closed.
        const localHasData = (ctx.treeNodes && ctx.treeNodes.length > 0) || !!ctx.smiles;
        if (!localHasData) {
          const project = projectData.projectsRef.current.find(p => p.id === projectId);
          const dbExp = project?.experiments.find(e => e.id === experimentId);
          const dbHasData = dbExp && ((dbExp.treeNodes && dbExp.treeNodes.length > 0) || !!dbExp.smiles);
          if (dbHasData) {
            console.warn('saveStateToExperiment: skipping save — local state is empty but DB has data for', experimentId);
            return false;
          }
        }

        // Guard: preserve sidebar messages from DB if local has none
        // but DB has reasoning messages (e.g. after browser reopen).
        const localSidebarMsgs = ctx.sidebarState?.messages ?? [];
        if (localSidebarMsgs.length === 0) {
          const project = projectData.projectsRef.current.find(p => p.id === projectId);
          const dbExp = project?.experiments.find(e => e.id === experimentId);
          const dbSidebarMsgs = dbExp?.sidebarState?.messages ?? [];
          if (dbSidebarMsgs.length > 0) {
            console.log(`saveStateToExperiment: local has 0 sidebar msgs but DB has ${dbSidebarMsgs.length} — preserving DB sidebar`);
            ctx.sidebarState = dbExp!.sidebarState;
          }
        }

        projectManagement.updateExperiment(projectId, ctx);
        return true;
      } catch (e) {
        console.warn('saveStateToExperiment: getContext failed, skipping save:', e);
        return false;
      }
    }
    return false;
  }, [projectSidebar.selectionRef, projectManagement, getContext, projectData]);

  const runComputation = async (): Promise<void> => {
    setSidebarOpen(true);

    // Default experiment names
    let experimentName = null;
    if (problemType === "optimization") {
      const propertyName = propertyType === "custom" ? customPropertyName : PROPERTY_NAMES[propertyType];
      experimentName = `Optimizing ${propertyName} for ${smiles}`;
    } else if (problemType === "retrosynthesis") {
      experimentName = `Synthesizing ${smiles}`;
    }

    // Check if we need to create project and/or experiment
    if (!projectSidebar.selectionRef.current.projectId) {
      // No project at all - create both project and experiment
      const now = new Date();
      const month = String(now.getMonth() + 1).padStart(2, '0');
      const day = String(now.getDate()).padStart(2, '0');
      const year = String(now.getFullYear()).slice(-2);
      const hours = String(now.getHours()).padStart(2, '0');
      const minutes = String(now.getMinutes()).padStart(2, '0');
      const timestamp = `${month}/${day}/${year} ${hours}:${minutes}`;

      const projectName = `Project ${timestamp}`;
      if (experimentName === null) {
        experimentName = `Experiment 1`;
      }
      try {
        const { projectId, experimentId } = await projectManagement.createProjectAndExperiment(
          projectName,
          experimentName
        );

        projectSidebar.setSelection({ projectId, experimentId });
        await new Promise(resolve => setTimeout(resolve, 100));
      } catch (error) {
        console.error('Error creating project:', error);
        alert('Failed to create project');
        return;
      }
    } else if (!projectSidebar.selectionRef.current.experimentId) {
      // Project exists but no experiment - create just an experiment
      const projectId = projectSidebar.selectionRef.current.projectId!;

      // Find the project to count existing experiments
      const project = projectData.projectsRef.current.find(p => p.id === projectId);
      if (experimentName === null) {
        const experimentCount = project ? project.experiments.length + 1 : 1;
        experimentName = `Experiment ${experimentCount}`;
      }

      try {
        const experiment = await projectManagement.createExperiment(projectId, experimentName);
        projectSidebar.setSelection({ projectId, experimentId: experiment.id });
        await new Promise(resolve => setTimeout(resolve, 100));
      } catch (error) {
        console.error('Error creating experiment:', error);
        alert('Failed to create experiment');
        return;
      }
    }

    setIsComputing(true);
    setTreeNodes([]);
    setEdges([]);
    graphState.setOffset({ x: 50, y: 50 });
    graphState.setZoom(1);

    const message: WebSocketMessageToServer = {
      action: 'compute',
      smiles,
      problemType,
      propertyType,
      customPropertyName,
      customPropertyDesc,
      customPropertyAscending,
      systemPrompt,
      userPrompt: problemPrompt,
      customization,
      experimentId: projectSidebar.selectionRef.current.experimentId || undefined,
    };

    wsRef.current?.send(JSON.stringify(message));
  };

  const reconnectingRef = useRef(false);

  const reconnectWS = (): void => {
    if (reconnectingRef.current) return; // Prevent overlapping reconnects
    reconnectingRef.current = true;

    // Close any existing socket (may have been orphaned by StrictMode cleanup)
    if (wsRef.current &&
        (wsRef.current.readyState === WebSocket.OPEN ||
         wsRef.current.readyState === WebSocket.CONNECTING)) {
      wsRef.current.close();
    }
    wsRef.current = null;  // Ensure old callbacks become no-ops

    setWsReconnecting(true);

    const socket = new WebSocket(WS_SERVER);
    wsRef.current = socket;

    socket.onopen = () => {
      reconnectingRef.current = false; // Clear guard
      console.log('WebSocket connected');
      setWsConnected(true);
      setWsReconnecting(false);
      setWsError('');
      resetLocalState();  // Clear UI state without killing server-side sessions
      isExperimentCompleteRef.current = false; // Reset experiment complete flag

      // Apply restored session from database if available (first connect only)
      if (restoredSessionRef.current) {
        console.log('Applying restored session after reset');
        applyRestoredSession(restoredSessionRef.current);
        restoredSessionRef.current = null; // Clear after first apply
      } else if (sessionLoadCompleteRef.current) {
        // Session load finished but nothing to restore (reconnect case,
        // or session was too old).  Load from current sidebar selection.
        loadStateFromCurrentExperiment();
      }
      // else: initial connect with session load still in flight.
      // The loadSavedSession callback will apply once it resolves.

      socket.send(JSON.stringify({ action: 'list-tools' }));
      socket.send(JSON.stringify({ action: 'get-username' }));
    };

    socket.onmessage = (event: MessageEvent) => {
      if (wsRef.current !== socket) return; // Ignore messages from old sockets

      const data: WebSocketMessage = JSON.parse(event.data);

      // Multi-experiment routing: determine if this message belongs to
      // the currently-selected experiment.  Messages without an
      // experimentId (e.g. tool lists, username) are always applied.
      const messageExperimentId = data.experimentId;
      const currentExperimentId = projectSidebar.selectionRef.current.experimentId;
      const isForCurrentExperiment =
        !messageExperimentId || messageExperimentId === currentExperimentId;

      if (data.type === 'node') {
        if (!isForCurrentExperiment) return; // Ignore background experiment nodes
        // Prevent duplicate nodes (can happen during session resume)
        setTreeNodes(prev => {
          const existingNode = prev.find(n => n.id === data.node!.id);
          if (existingNode) {
            console.log('Skipping duplicate node:', data.node!.id);
            return prev;
          }
          return [...prev, data.node!];
        });
        // Save state when new molecule is generated (if auto-checkpoint is enabled).
        // Only saveStateToExperiment is called here.  saveCheckpoint uses
        // getSessionState() which historically applied level filtering,
        // creating a race where filtered data could overwrite full data in
        // the DB (both endpoints write to the same experiment row).
        if (autoCheckpointEnabledRef.current) {
          setTimeout(() => {
            saveStateToExperimentRef.current();
          }, 100); // Small delay to ensure state is updated
        }
      } else if (data.type === 'stopped') {
        // Handle explicit stop from backend
        console.log('Computation stopped by backend', messageExperimentId);
        // Clean up computing map
        if (messageExperimentId) {
          computingExperimentsRef.current.delete(messageExperimentId);
          // Find project containing this experiment and update isRunning
          const project = projectData.projectsRef.current.find(p =>
            p.experiments.some(e => e.id === messageExperimentId)
          );
          if (project) {
            projectData.setExperimentRunning(project.id, messageExperimentId, false);
          }
        }
        if (!isForCurrentExperiment) return;
        setIsComputing(false);
        setIsComputingTemplates(false);
        unhighlightNodes();
        saveStateToExperiment();
      } else if (data.type === 'node_update') {
        if (!isForCurrentExperiment) return;
        const { id, ...restData } = data.node!;

        if (restData) {
           setIsComputing(true);
        }else {
           setIsComputing(false);
           setIsComputingTemplates(false);
        }
        setTreeNodes(prev => prev.map(n =>
          n.id === data.node!.id ? { ...n, ...restData } : n
        ));
      } else if (data.type === 'node_delete') {
        if (!isForCurrentExperiment) return;
        setTreeNodes(prev => {
          const descendants = findAllDescendants(data.node!.id, prev);
          return prev.filter(n => !descendants.has(n.id) && n.id !== data.node!.id);
        });
        setEdges(prev => prev.filter(e =>
          e.fromNode !== data.node!.id && e.toNode !== data.node!.id
        ));
      } else if (data.type === 'subtree_update') {
        if (!isForCurrentExperiment) return;
        const withNode = data.withNode || false;
        const { id, ...restData } = data.node!;
        setTreeNodes(prev => {
          const descendants = findAllDescendants(data.node!.id, prev);
          return prev.map(n =>
            (descendants.has(n.id) || (withNode && n.id === data.node!.id))
              ? { ...n, ...restData }
              : n
          );
        });
      } else if (data.type === 'edge') {
        if (!isForCurrentExperiment) return;
        // Prevent duplicate edges (can happen during session resume)
        setEdges(prev => {
          const existingEdge = prev.find(e => e.id === data.edge!.id);
          if (existingEdge) {
            console.log('Skipping duplicate edge:', data.edge!.id);
            return prev;
          }
          return [...prev, data.edge!];
        });
      } else if (data.type === 'edge_update') {
        if (!isForCurrentExperiment) return;
        const { id, ...restData } = data.edge!;
        setEdges(prev => prev.map(e =>
          e.id === data.edge!.id ? { ...e, ...restData } : e
         ));
      } else if (data.type === 'subtree_delete') {
        if (!isForCurrentExperiment) return;
        let descendantsSet: Set<string>;
        setTreeNodes(prev => {
          descendantsSet = findAllDescendants(data.node!.id, prev);
          return prev.filter(n => !descendantsSet.has(n.id));
        });
        setEdges(prev => prev.filter(e =>
          !descendantsSet!.has(e.fromNode) && !descendantsSet!.has(e.toNode)
        ));
      } else if (data.type === 'complete') {
        // Clean up computing map regardless of which experiment completed
        if (messageExperimentId) {
          computingExperimentsRef.current.delete(messageExperimentId);
          // Update isRunning flag for background experiment
          const project = projectData.projectsRef.current.find(p =>
            p.experiments.some(e => e.id === messageExperimentId)
          );
          if (project) {
            projectData.setExperimentRunning(project.id, messageExperimentId, false);
          }
        }
        if (!isForCurrentExperiment) {
          console.log('Background experiment completed:', messageExperimentId);
          // Refresh project data from DB so the sidebar_state
          // (including "Retrosynthesis Complete" / "Optimization Complete")
          // is in the in-memory cache when the user switches to this
          // experiment.  The server persists sidebar_state to DB before
          // sending 'complete', so the refresh will find the data.
          projectData.refreshProjects();
          return;
        }
        setIsComputing(false);
        isExperimentCompleteRef.current = true;  // Mark experiment as complete to stop checkpointing
        unhighlightNodes();
        saveStateToExperiment();  // Keep experiment up to date (localStorage / project data)

        // Persist final state to MariaDB so other browsers see the completed tree.
        // Use force=true to bypass the "no data" guard.
        sessionPersistence.saveSession(
          getSessionStateRef.current(),
          true, // force
        );

        // Clear server session ID on completion
        sessionPersistence.setServerSessionId(null);
        sessionPersistence.setSessionWasComputing(false);
      } else if (data.type === 'session_started') {
        // Store the server session ID for resume capability
        const serverSessionId = data.sessionId || (data as any).sessionId;
        console.log('Computation session started:', serverSessionId);
        // Track in computing experiments map — prefer experimentId from
        // the server message (includes it in the response), fall back to
        // whatever the user currently has selected.
        const expId = messageExperimentId || currentExperimentId;
        if (expId && serverSessionId) {
          computingExperimentsRef.current.set(expId, serverSessionId);
        }
        sessionPersistence.setServerSessionId(serverSessionId);
        sessionPersistence.setSessionWasComputing(true);
      } else if (data.type === 'session_resumed') {
        // Successfully resumed a session
        const serverSessionId = data.sessionId || (data as any).sessionId;
        console.log('Computation session resumed:', serverSessionId);
        // Track in computing experiments map
        const expId = messageExperimentId || currentExperimentId;
        if (expId && serverSessionId) {
          computingExperimentsRef.current.set(expId, serverSessionId);
        }
        sessionPersistence.setServerSessionId(serverSessionId);
        sessionPersistence.setSessionWasComputing(true);
        setIsComputing(true);
      } else if (data.type === 'session_status') {
        // Handle session status updates during resume
        const status = (data as any).status;
        console.log('Session status:', status);
        if (status === 'complete' || status === 'cancelled') {
          setIsComputing(false);
          // Clean up computing map
          if (messageExperimentId) {
            computingExperimentsRef.current.delete(messageExperimentId);
          }
          sessionPersistence.setServerSessionId(null);
          sessionPersistence.setSessionWasComputing(false);
        }
      } else if (data.type === 'session_not_found') {
        // Server couldn't find the session to resume
        console.log('Session not found, cannot resume');
        setIsComputing(false);
        // Clean up computing map for this experiment
        if (messageExperimentId) {
          computingExperimentsRef.current.delete(messageExperimentId);
        }
        sessionPersistence.setServerSessionId(null);
        sessionPersistence.setSessionWasComputing(false);
        sessionPersistence.setPendingResume(null);
      } else if (data.type === 'response') {
        // Only add reasoning messages for the currently-viewed experiment.
        // Background experiment messages are tracked server-side and
        // persisted to sidebar_state, so they appear when the user
        // switches to that experiment.
        if (!isForCurrentExperiment) return;
        addSidebarMessage(data.message!);
        console.log('Server response:', data.message);
      } else if (data.type === 'available-tools-response') {
        const newTools = data.tools || [];
        setAvailableTools(newTools);
        setAvailableToolsMap(
          newTools.map((server: Tool, index: number) => ({
            id: index,
            tool_server: server
          }))
        );
        setSelectedTools([]);
      } else if (data.type === 'server-update-orchestrator-settings') {
        // Handle orchestrator settings updates from server
        const newSettings: OrchestratorSettings = {
          backend: data.orchestratorSettings!.backend,
          useCustomUrl: data.orchestratorSettings!.useCustomUrl,
          customUrl: data.orchestratorSettings!.customUrl,
          model: data.orchestratorSettings!.model,
          // Don't take the use custom model field from the backend
          // Check the model against the list of models in copilot
          // useCustomModel: data.orchestratorSettings.useCustomModel,
          apiKey: data.orchestratorSettings!.apiKey,
          backendLabel: data.orchestratorSettings!.backendLabel
        };
        setOrchestratorSettings(newSettings);
        console.log('Updating the orchestrator settings ', newSettings);
        localStorage.setItem('orchestratorSettings', JSON.stringify(newSettings));
      } else if (data.type === 'error') {
        console.error(data.message);
        alert("Server error: " + data.message);
      } else if (data.type === 'save-context-response') {
        saveFullContext(data.experimentContext!);
      } else if (data.type === 'get-username-response') {
        setUsername(data.username!);
      }
    };

    socket.onerror = (error: Event) => {
      reconnectingRef.current = false; // Clear guard on error
      console.error('WebSocket error:', error);
      // Only update state if this is the current socket
      if (wsRef.current === socket) {
        setWsReconnecting(false);
        setIsComputing(false);
        setIsComputingTemplates(false);
        setWsError((error as any).message || 'Connection failed');
        setAvailableTools([]);
        setSelectedTools([]);
        setAvailableToolsMap([]);
        // All sessions lose their WS on error; server continues headless
        computingExperimentsRef.current.clear();
      }
    };

    socket.onclose = () => {
      reconnectingRef.current = false; // Clear guard on close
      console.log('WebSocket closed');
      // Only clear state if this is the current socket
      if (wsRef.current === socket) {
        wsRef.current = null;
        setWsConnected(false);
        setIsComputing(false);
        setIsComputingTemplates(false);
        setWsReconnecting(false);
        setAvailableTools([]);
        setSelectedTools([]);
        setAvailableToolsMap([]);
        // All sessions lose their WS on close; server continues headless
        computingExperimentsRef.current.clear();
      }
    };
  };

  // Connect WebSocket on mount
  useEffect(() => {
    reconnectWS();

    return () => {
      // Clear the guard so a StrictMode re-mount can reconnect.
      // The onclose handler only fires asynchronously, so the guard would
      // otherwise still be set when the second mount calls reconnectWS().
      reconnectingRef.current = false;
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    };
  }, []);

  // Reset local UI state only (no server message).
  // Used on WebSocket open so a session restore or resume can follow.
  const resetLocalState = (): void => {
    setTreeNodes([]);
    setEdges([]);
    setIsComputing(false);
    isExperimentCompleteRef.current = false;
    graphState.setOffset({ x: 50, y: 50 });
    graphState.setZoom(1);
    setContextMenu({node: null, isReaction: false, x:0, y:0});
    setCustomQueryModal(null);
    metricsDashboardState.setMetricsHistory([]);
    sidebarState.setMessages([]);
    setSaveDropdownOpen(false);
    setCheckpointDropdownOpen(false);
    sidebarState.setSourceFilterOpen(false);
    setWsTooltipPinned(false);
    // Clear and close reaction alternatives sidebar
    setReactionSidebarOpen(false);
    setSelectedReactionNode(null);
    setIsComputingTemplates(false);
  };

  // Full reset: clear local state AND tell the server to discard its session.
  // Only used for explicit user actions (Reset button, new computation).
  const reset = (): void => {
    resetLocalState();
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      wsRef.current.send(JSON.stringify({ action: 'reset' }));
    }
  };

  // Apply restored session state from database (called after reset in onopen)
  const applyRestoredSession = (savedState: SessionState): void => {
    console.log('Applying restored session state:', savedState);
    
    // Restore sidebar selection so saveStateToExperiment works after restore
    if (savedState.projectId || savedState.experimentId) {
      projectSidebar.setSelection({
        projectId: savedState.projectId || null,
        experimentId: savedState.experimentId || null,
      });
    }

    // If the DB (projectsRef) has a richer version of this experiment
    // (e.g. save_session_to_db ran after the browser closed but before
    // this page load), prefer the DB data over the stale session.
    if (savedState.projectId && savedState.experimentId) {
      const proj = projectData.projectsRef.current.find(p => p.id === savedState.projectId);
      const dbExp = proj?.experiments.find(e => e.id === savedState.experimentId);
      if (dbExp) {
        const dbNodeCount = dbExp.treeNodes?.length ?? 0;
        const sessionNodeCount = savedState.treeNodes?.length ?? 0;
        const dbSidebarMsgCount = dbExp.sidebarState?.messages?.length ?? 0;
        const sessionSidebarMsgCount = savedState.sidebarMessages?.length ?? 0;
        if (dbNodeCount > sessionNodeCount || (dbNodeCount === sessionNodeCount && dbSidebarMsgCount > sessionSidebarMsgCount)) {
          console.log(
            `applyRestoredSession: DB has ${dbNodeCount} nodes/${dbSidebarMsgCount} sidebar msgs vs session's ${sessionNodeCount}/${sessionSidebarMsgCount}. Loading from DB instead.`
          );
          loadContext(dbExp);
          // If the DB experiment has no problemType but the session
          // does, prefer the session's value (the DB row may not have
          // been updated yet when the projects were fetched).
          if (!dbExp.problemType && savedState.problemType) {
            setProblemType(savedState.problemType);
          }
          return;
        }
      }
    }
    
    // Restore core state
    if (savedState.smiles) setSmiles(savedState.smiles);
    if (savedState.problemType) setProblemType(savedState.problemType);
    if (savedState.systemPrompt) setSystemPrompt(savedState.systemPrompt);
    if (savedState.problemPrompt) setProblemPrompt(savedState.problemPrompt);
    if (savedState.autoZoom !== undefined) setAutoZoom(savedState.autoZoom);
    
    // Restore tree nodes and edges, applying relayout to ensure correct positions
    if (savedState.treeNodes && savedState.treeNodes.length > 0) {
      const nodesToRestore = savedState.treeNodes;
      const edgesToRestore = savedState.edges || [];
      
      // Apply relayout to fix any position issues
      const [relayoutedNodes, relayoutedEdges] = relayoutTree(nodesToRestore, edgesToRestore);
      
      setTreeNodes(relayoutedNodes);
      setEdges(relayoutedEdges);
    } else if (savedState.edges && savedState.edges.length > 0) {
      setEdges(savedState.edges);
    }
    
    // Restore graph state (offset/zoom will be overridden by auto-zoom if enabled)
    if (savedState.offset) graphState.setOffset(savedState.offset);
    if (savedState.zoom) graphState.setZoom(savedState.zoom);
    
    // Restore metrics
    if (savedState.metricsHistory) metricsDashboardState.setMetricsHistory(savedState.metricsHistory);
    if (savedState.visibleMetrics) metricsDashboardState.setVisibleMetrics(savedState.visibleMetrics);
    
    // Restore property optimization settings
    if (savedState.propertyType) setPropertyType(savedState.propertyType);
    if (savedState.customPropertyName) setCustomPropertyName(savedState.customPropertyName);
    if (savedState.customPropertyDesc) setCustomPropertyDesc(savedState.customPropertyDesc);
    if (savedState.customPropertyAscending !== undefined) setCustomPropertyAscending(savedState.customPropertyAscending);
    
    // Restore sidebar state
    if (savedState.sidebarMessages && savedState.sidebarMessages.length > 0) {
      sidebarState.setMessages(savedState.sidebarMessages);
      setSidebarOpen(true);
    }
    if (savedState.sidebarVisibleSources) {
      sidebarState.setVisibleSources(savedState.sidebarVisibleSources);
    }
    
    // Trigger auto-zoom fit after restore if autoZoom is enabled
    // Use multiple requestAnimationFrame to ensure DOM is updated
    if (savedState.autoZoom && savedState.treeNodes && savedState.treeNodes.length > 0) {
      setTimeout(() => {
        // Toggle autoZoom off and on to trigger the fitToView effect in graph.tsx
        setAutoZoom(false);
        requestAnimationFrame(() => {
          setAutoZoom(true);
        });
      }, 100);
    }
    
    console.log('Session state restored successfully');
  };

  const unhighlightNodes = (): void => {
    setTreeNodes(prev => prev.map(n =>
      n.highlight === "yellow" ? { ...n, highlight: "normal" } : n
    ));
  };

  const stop = (): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    // Include the sessionId so the server can stop the correct session
    // (important when multiple experiments are running concurrently).
    const experimentId = projectSidebar.selectionRef.current.experimentId;
    const sessionId = experimentId
      ? computingExperimentsRef.current.get(experimentId)
      : undefined;
    wsRef.current.send(JSON.stringify({ action: 'stop', sessionId }));
  };

  useEffect(() => {
    getContextRef.current = () => {
      const projectId = projectSidebar.selectionRef.current.projectId;
      const experimentId = projectSidebar.selectionRef.current.experimentId;
      const project = projectData.projectsRef.current.find(p => p.id === projectId);
      if (project) {
        const experiment = project.experiments.find(e => e.id === experimentId);
        if (experiment) {
          return {
            ...experiment,
            smiles,
            problemType,
            systemPrompt,
            problemPrompt,
            propertyType,
            customPropertyName,
            customPropertyDesc,
            customPropertyAscending,
            customization,
            treeNodes: treeNodesRef.current,
            edges: edgesRef.current,
            metricsHistory: metricsDashboardState.metricsHistory,
            visibleMetrics: metricsDashboardState.visibleMetrics,
            graphState,
            autoZoom,
            sidebarState: sidebarStateRef.current
          };
        }
      }
      throw "No experiment found";
    };
  }, [smiles, problemType, graphState, metricsDashboardState, autoZoom,
      systemPrompt, problemPrompt, propertyType, customPropertyName, customPropertyDesc, customPropertyAscending, customization, projectData, projectSidebar]);

  // Function to get current session state for auto-save to MariaDB
  const getSessionState = useCallback((): SessionState => {
    const projectId = projectSidebar.selectionRef.current.projectId;
    const experimentId = projectSidebar.selectionRef.current.experimentId;
    const project = projectData.projectsRef.current.find(p => p.id === projectId);
    const experiment = project?.experiments.find(e => e.id === experimentId);
    
    // Get current tree nodes and edges.  Save the full tree without
    // level filtering.  Previously, the deepest level was excluded during
    // computation to avoid saving an incomplete frontier.  However,
    // duplicate-node dedup on resume (the node handler checks
    // prev.find(n => n.id === data.node.id)) already handles this, and
    // the filtering created a race condition: the session-save endpoint
    // wrote filtered data while the experiment-save endpoint wrote full
    // data to the same DB row.  Whichever settled last determined the
    // final state, causing intermittent data loss.
    const nodesToSave = treeNodesRef.current;
    // Filter out transient 'computing' edges that have no value on restore.
    const edgesToSave = edgesRef.current.filter(e => e.status !== 'computing');
    
    return {
      // Project/Experiment identification
      projectId: projectId || undefined,
      projectName: project?.name,
      experimentId: experimentId || undefined,
      experimentName: experiment?.name,
      
      // Core experiment state
      smiles,
      problemType,
      systemPrompt,
      problemPrompt,
      autoZoom,
      treeNodes: nodesToSave,
      edges: edgesToSave,
      offset: graphState.offset,
      zoom: graphState.zoom,
      metricsHistory: metricsDashboardState.metricsHistory,
      visibleMetrics: metricsDashboardState.visibleMetrics,
      isComputing,
      serverSessionId: sessionPersistence.serverSessionId,
      
      // Property optimization
      propertyType,
      customPropertyName,
      customPropertyDesc,
      customPropertyAscending,
      
      // Sidebar state
      sidebarMessages: sidebarStateRef.current.messages,
      sidebarSourceFilterOpen: sidebarStateRef.current.sourceFilterOpen,
      sidebarVisibleSources: sidebarStateRef.current.visibleSources,
    };
  }, [smiles, problemType, systemPrompt, problemPrompt, autoZoom, graphState, 
      metricsDashboardState, isComputing, sessionPersistence.serverSessionId,
      propertyType, customPropertyName, customPropertyDesc, customPropertyAscending,
      projectSidebar, projectData]);

  // Keep a ref to getSessionState for use in event handlers
  const getSessionStateRef = useRef(getSessionState);
  useEffect(() => {
    getSessionStateRef.current = getSessionState;
  }, [getSessionState]);

  // Function to save checkpoint (called when new molecules are generated)
  const saveCheckpoint = useCallback(() => {
    // Don't save checkpoint if experiment is complete
    if (isExperimentCompleteRef.current) {
      console.log('Skipping checkpoint - experiment is complete');
      return;
    }
    // Update the current session in place (checkpoint: false) so we don't
    // create a new experiment entry for every batch of molecules.
    sessionPersistence.saveSession(
      getSessionStateRef.current(),
      false
    );
  }, [sessionPersistence]);

  // Function to manually trigger a checkpoint (always saves, ignores experiment complete status)
  const triggerManualCheckpoint = useCallback(() => {
    console.log('Manual checkpoint triggered');
    const state = getSessionStateRef.current();
    // For manual checkpoint, save all current state regardless of experiment status
    sessionPersistence.saveSession(
      state,
      true,
      { checkpoint: true, name: `Manual Checkpoint ${new Date().toISOString().replace('T', ' ').slice(0, 19)}` }
    );
    setCheckpointDropdownOpen(false);
  }, [sessionPersistence]);

  // Update the saveCheckpoint ref so WebSocket handler can use it
  useEffect(() => {
    saveCheckpointRef.current = saveCheckpoint;
  }, [saveCheckpoint]);

  // Update the saveStateToExperiment ref so WebSocket handler can persist to localStorage
  useEffect(() => {
    saveStateToExperimentRef.current = saveStateToExperiment;
  }, [saveStateToExperiment]);

  // Set up checkpoint on page unload (saves to database via sendBeacon)
  useAutoSave(sessionPersistence, getSessionState);

  // Also save experiment state to localStorage on page unload so it survives
  // browser restarts even when the database is unavailable.
  // This must be synchronous -- async functions may not complete during unload.
  useEffect(() => {
    const handleUnloadSaveToLocalStorage = () => {
      try {
        const projectId = projectSidebar.selectionRef.current.projectId;
        const experimentId = projectSidebar.selectionRef.current.experimentId;
        if (!projectId || !experimentId) return;

        const context = getContextRef.current();
        const storageKey = 'flask_copilot_projects';
        const raw = localStorage.getItem(storageKey);
        if (!raw) return;

        const projects = JSON.parse(raw);
        const project = projects.find((p: { id: string }) => p.id === projectId);
        if (!project) return;

        const expIndex = project.experiments.findIndex((e: { id: string }) => e.id === experimentId);
        if (expIndex === -1) return;

        project.experiments[expIndex] = { ...context, lastModified: new Date().toISOString() };
        project.lastModified = new Date().toISOString();
        localStorage.setItem(storageKey, JSON.stringify(projects));
      } catch {
        // Silently ignore -- experiment may not exist yet
      }
    };
    window.addEventListener('beforeunload', handleUnloadSaveToLocalStorage);
    return () => {
      window.removeEventListener('beforeunload', handleUnloadSaveToLocalStorage);
    };
  }, []);

  // Load session from database on initial mount (runs only once)
  const sessionLoadAttemptedRef = useRef(false);
  
  useEffect(() => {
    // Only attempt to load once
    if (sessionLoadAttemptedRef.current) return;
    sessionLoadAttemptedRef.current = true;
    
    const loadSavedSession = async () => {
      const savedState = await sessionPersistence.loadSession();
      sessionLoadCompleteRef.current = true;

      if (savedState) {
        console.log('Session loaded from database, storing for restore after WebSocket connects:', savedState);
        
        // Store the session in ref so it can be applied after WebSocket connects
        restoredSessionRef.current = savedState;
        
        // Check WebSocket readyState directly via ref -- the wsConnected state
        // variable is stale in this closure since both effects fire on mount.
        if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
          console.log('WebSocket already connected, applying restored session now');
          applyRestoredSession(savedState);
          restoredSessionRef.current = null;
        }
        // Otherwise, it will be applied in the onopen handler
      } else {
        // No saved session.  If the WS already opened (won the race) and
        // skipped loadStateFromCurrentExperiment because the load was
        // pending, do it now.
        if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
          loadStateFromCurrentExperiment();
        }
      }
    };
    
    loadSavedSession();

    return () => {
      // Reset guard so StrictMode re-mount can load the session again
      sessionLoadAttemptedRef.current = false;
      sessionLoadCompleteRef.current = false;
      restoredSessionRef.current = null;
    };
  }, []); // Empty dependency array - only run once on mount

  // Handle computation resume after WebSocket connects.
  // Only attempt resume once we have both a live WS connection AND the
  // session data from the DB (sessionLoaded).  Without the sessionLoaded
  // guard, the effect may fire before loadSession() resolves, see no
  // pending resume, and prematurely mark resumeAttempted = true.
  useEffect(() => {
    if (!wsConnected || sessionPersistence.resumeAttempted) return;

    // Wait until the session has been loaded from the database so
    // pendingResume has been populated (or confirmed absent).
    if (!sessionPersistence.sessionLoaded) return;
    
    const pendingResume = sessionPersistence.getPendingResume();
    // Only attempt resume if we have a valid server session ID
    if (pendingResume && pendingResume.sessionId && pendingResume.sessionId.length > 0) {
      console.log('Attempting to resume computation with session:', pendingResume.sessionId);
      sessionPersistence.setResumeAttempted(true);
      
      // Send resume request to server
      if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
        wsRef.current.send(JSON.stringify({
          action: 'resume_session',
          sessionId: pendingResume.sessionId
        }));
        setIsComputing(true);
      }
    } else {
      // No valid session to resume, mark as attempted to prevent future attempts
      sessionPersistence.setResumeAttempted(true);
    }
  }, [wsConnected, sessionPersistence]);


  const saveTree = (): void => {
    const data = { version: '1.0', type: 'tree', timestamp: new Date().toISOString(), smiles, nodes: treeNodes, edges };
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `molecule-tree-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
    setSaveDropdownOpen(false);
  };

  const requestSaveContext = (): void => {
    // We need to request the up-to-date Project object from the server before saving
    sendMessageToServer('save-context');
  };

  const saveFullContext = (experimentContext: string): void => {
    const data = { lastModified: new Date().toISOString(), smiles, problemType, systemPrompt, problemPrompt, propertyType, customPropertyName,
                   customPropertyDesc, customPropertyAscending, customization, treeNodes, edges, graphState, metricsDashboardState, sidebarState, experimentContext };
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `experiment-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
    setSaveDropdownOpen(false);
  };

  const loadContextFromFile = (): void => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.json';
    input.onchange = (e: Event) => {
      const file = (e.target as HTMLInputElement).files?.[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (e: ProgressEvent<FileReader>) => {
        try {
          const data = JSON.parse(e.target?.result as string) as Experiment;
          loadContext(data);
        } catch (error) {
          alert('Error loading file: ' + (error as Error).message);
        }
      };
      reader.readAsText(file);
    };
    input.click();
  };

  const savePrompts = (newSystemPrompt: string, newProblemPrompt: string): void => {
    setSystemPrompt(newSystemPrompt);
    setProblemPrompt(newProblemPrompt);
    setEditPromptsModal(false);
  };

  const saveCustomProperty = (newPropertyName: string, newPropertyDesc: string, newPropertyAscending: boolean): void => {
    setCustomPropertyName(newPropertyName);
    setCustomPropertyDesc(newPropertyDesc);
    setCustomPropertyAscending(newPropertyAscending);
    setEditPropertyModal(false);
  };

  const resetProblemType = (problem_type: string): void => {
    setSystemPrompt('');
    setProblemPrompt('');
    // Set default metrics
    if (problem_type === 'retrosynthesis') {
      metricsDashboardState.setVisibleMetrics({cost: false, bandgap: false, sascore: false, yield: true, density: false});
    } else if (problem_type === 'optimization') {
      metricsDashboardState.setVisibleMetrics({cost: false, bandgap: false, sascore: true, yield: false, density: true});
    } else if (problem_type === 'custom') {
      metricsDashboardState.setVisibleMetrics({cost: false, bandgap: false, sascore: false, yield: false, density: false});
      setSystemPrompt(DEFAULT_CUSTOM_SYSTEM_PROMPT);
    }
    setProblemType(problem_type);
  };

  const createNewRetrosynthesisExperiment = async (startingSmiles: string): Promise<void> => {
    // Close context menu
    setContextMenu({node: null, isReaction: false, x: 0, y: 0});

    saveStateToExperiment();

    // Find the project to create a new experiment in
    const projectId = projectSidebar.selectionRef.current.projectId!;
    const experimentName = `Synthesizing ${startingSmiles}`;

    let experiment = undefined;
    try {
      experiment = await projectManagement.createExperiment(projectId, experimentName);
      projectSidebar.setSelection({ projectId, experimentId: experiment.id });
      await new Promise(resolve => setTimeout(resolve, 100));
    } catch (error) {
      console.error('Error creating experiment:', error);
      alert('Failed to create experiment');
      return;
    }

    reset();
    setSmiles(startingSmiles);
    resetProblemType("retrosynthesis");
    saveStateToExperiment();
  }

  const handleNodeClick = (e: React.MouseEvent<HTMLDivElement>, node: TreeNode): void => {
      e.stopPropagation();
      if (isComputing) return; // Don't open menu while computing
      setContextMenu({
          node,
          isReaction: false,
          x: e.clientX,
          y: e.clientY
      });
  };

  const handleReactionClick = (e: React.MouseEvent<HTMLDivElement>, node: TreeNode): void => {
      e.stopPropagation();
      if (isComputing) return; // Don't open menu while computing
      setContextMenu({
          node,
          isReaction: true,
          x: e.clientX,
          y: e.clientY
      });
  };

  // Handle minimizing the prompt breakpoint modal
  const handleMinimizePromptModal = useCallback(() => {
    setDebugModalMinimized(true);
    // Keep isComputing true to show that we're still in a computing state
    // Don't clear promptBreakpoint or editedPrompt - preserve the state
  }, []);

  // Handle reopening the minimized prompt breakpoint modal
  const handleReopenPromptModal = useCallback(() => {
    setDebugModalMinimized(false);
  }, []);

  // Handle prompt approval/modification
  const handlePromptBreakpointResponse = useCallback(() => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }

    wsRef.current.send(
      JSON.stringify({
        action: 'prompt-breakpoint-response',
        prompt: editedPrompt,
        metadata: promptBreakpoint?.metadata,
      })
    );

    setPromptBreakpoint(null);
    setEditedPrompt('');
    setDebugModalMinimized(false);
    setIsComputing(true); // Resume computation
  }, [editedPrompt, promptBreakpoint]);

  const sendMessageToServer = useCallback((message: string, data?: Omit<WebSocketMessageToServer, 'action'>): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    const msg: WebSocketMessageToServer = {
      action: message,
      ...data
    };
    wsRef.current.send(JSON.stringify(msg));
    setContextMenu({node: null, isReaction: false, x: 0, y: 0});
  }, []);

  const handleReactionCardClick = useCallback((node: TreeNode) => {
    if (isComputing) return;
    setSelectedReactionNode(node);
    setReactionSidebarOpen(true);
  }, [isComputing]);  // Only depend on isComputing boolean

  const handleSelectAlternative = useCallback((alt: ReactionAlternative) => {
    const nodeId = selectedReactionNode?.id;
    if (!nodeId) return;
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) return;

    // Tell backend that the alternative subtree has been chosen
    wsRef.current.send(JSON.stringify({
      action: "set-reaction-alternative",
      nodeId: nodeId,
      alternativeId: alt.id
    }));

    // Don't close the sidebar - let user see the active status update
    setIsComputing(true);
  }, [selectedReactionNode?.id]);  // Only depend on the ID, not the whole node

  const handleCloseReactionAlternativesSidebar = useCallback(() => {
    setReactionSidebarOpen(false);
  }, []);  // No dependencies - just a state setter

  const handleComputeTemplates = useCallback(() => {
    const nodeId = selectedReactionNode?.id;
    const smiles = selectedReactionNode?.smiles;
    if (!nodeId || !smiles) return;

    setIsComputingTemplates(true);
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      wsRef.current.send(JSON.stringify({
        action: "compute-reaction-templates",
        nodeId: nodeId,
        smiles: smiles
      }));
    }
  }, [selectedReactionNode?.id, selectedReactionNode?.smiles]);  // Only depend on primitive values

  const handleComputeFlaskAI = useCallback((customPrompt: boolean) => {
    const nodeId = selectedReactionNode?.id;
    if (!nodeId) return;

    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      if (customPrompt) {
        handleCustomQuery(selectedReactionNode!, "compute-reaction-from")
      } else {
        wsRef.current.send(JSON.stringify({
          action: "compute-reaction-from",
          nodeId: nodeId
        }));
      }
    }
    setIsComputing(true);
  }, [selectedReactionNode?.id]);  // Only depend on the ID

  const stableAlternatives = useMemo(() => {
    return selectedReactionNode?.reaction?.alternatives || [];
  }, [selectedReactionNode?.id, selectedReactionNode?.reaction?.alternatives]);


  const handleCustomQuery = (node: TreeNode, queryType: string | null): void => {
    setCustomQueryModal(node);
    setCustomQueryText('');
    setCustomQueryType(queryType);
    setContextMenu({node: null, isReaction: false, x: 0, y: 0});
  };

  const submitCustomQuery = (): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }

    let propertyDetails = {};
    if (problemType === "optimization") {
      propertyDetails = {
        propertyType,
        customPropertyName,
        customPropertyDesc,
        customPropertyAscending,
        smiles: customQueryModal?.smiles,
        xpos: customQueryModal?.x,
        customization
      };
    }

    const message: WebSocketMessageToServer = {
      action: customQueryType ?? (problemType === "optimization" ? "optimize-from" : "recompute-reaction"),
      nodeId: customQueryModal?.id,
      smiles: customQueryModal?.smiles,
      query: customQueryText,
      ...propertyDetails
    };
    wsRef.current.send(JSON.stringify(message));

    setCustomQueryModal(null);
    setCustomQueryText('');
    setCustomQueryType(null);
    setIsComputing(true); // If expecting new nodes
  };

  const addSidebarMessage = (message: SidebarMessage): void => {
    message.id = Date.now();
    message.timestamp = new Date().toISOString();
    if (!message.source) {
      message.source = "Backend";
    }

    // Update the ref synchronously so that getSessionState() (called from
    // the 'complete' handler in the same event-loop tick) always sees the
    // latest messages, even before React commits the batched setState.
    const updatedMessages = [...sidebarStateRef.current.messages, message];
    sidebarStateRef.current = { ...sidebarStateRef.current, messages: updatedMessages };

    sidebarState.setMessages(prev => [...prev, message]);
    setSidebarOpen(true);

    sidebarState.setVisibleSources(prev => {
      if (!(message.source in prev)) {
        return { ...prev, [message.source]: true };
      }
      return prev;
    });
  };

  return (
    <div className="app-background">
      <div className="main-container">
        <ProjectSidebar
          projectData={projectData}
          isOpen={projectSidebar.isOpen}
          onToggle={projectSidebar.toggleSidebar}
          selection={projectSidebar.selection}
          onSelectionChange={projectSidebar.setSelection}
          onLoadContext={loadContextFromExperiment}
          onSaveContext={saveStateToExperiment}
          onReset={resetLocalState}
          isComputing={isComputing}
        />
        <div className="content-wrapper">
          <div className="w-full">
            <div className="content-header">
              { /* Left logos */ }
              <div className="text-white">
                <div className="w-full flex app-logo">
                  <svg width="40" height="40" viewBox="0 0 28 28" fill="none">
                    <path d="M13.967 0C6.65928 0 0.646366 5.60212 0 12.7273H2.77682C3.25522 11.7261 4.27793 11.0303 5.46222 11.0303C7.10365 11.0303 8.43891 12.3624 8.43891 14C8.43891 15.6376 7.10365 16.9697 5.46222 16.9697C4.27793 16.9697 3.25522 16.2739 2.77682 15.2727H0C0.646366 22.3979 6.65928 28 13.967 28C21.7043 28 28 21.7191 28 14C28 6.28091 21.7043 0 13.967 0ZM5.46222 19.5152C8.5112 19.5152 10.9904 17.0418 10.9904 14C10.9904 10.9582 8.5112 8.48485 5.46222 8.48485C5.32189 8.48485 5.18156 8.49121 5.04336 8.50182C6.3042 7.43273 7.935 6.78788 9.71463 6.78788C13.7013 6.78788 16.9437 10.0227 16.9437 14C16.9437 17.9773 13.7013 21.2121 9.71463 21.2121C7.935 21.2121 6.3042 20.5652 5.04336 19.4982C5.18156 19.5088 5.32189 19.5152 5.46222 19.5152ZM13.967 25.4545C11.6112 25.4545 9.42122 24.7418 7.59693 23.5242C8.27944 23.6749 8.98747 23.7576 9.71463 23.7576C15.1067 23.7576 19.4952 19.3794 19.4952 14C19.4952 8.62061 15.1067 4.24242 9.71463 4.24242C8.98747 4.24242 8.27944 4.32515 7.59693 4.47576C9.42122 3.25818 11.6112 2.54545 13.967 2.54545C20.2989 2.54545 25.4486 7.68303 25.4486 14C25.4486 20.317 20.2989 25.4545 13.967 25.4545Z" fill="white"/>
                  </svg>
                  <p className="text-center font-['Geist',sans-serif] text-[32px] leading-[1.3] font-medium text-nowrap whitespace-pre text-white noselect"> Genesis Mission</p>
                </div>
              </div>

              { /* Right logos */ }
              <div className="app-logo-right group flex">
                  <svg version="1.1" id="Layer_1" className="logo-svg" viewBox="0 0 40 40">
                    <g>
                      <rect x="1.73" y="0.01" fill="#FFFFFF" width="34.19" height="34.19"/>
                      <path fill="#1E59AE" d="M35.92,0.01v17.53H18.95V0.01H35.92z M15.88,21.82c-1.12-0.07-1.72-0.78-1.79-2.1V0.01h-0.76v19.73
              c0.09,1.72,1,2.75,2.53,2.84h15.28l-4.83,5l-11.79,0c-3.04-0.36-6.22-2.93-6.14-6.98V0.01H7.6V20.6c-0.09,4.49,3.45,7.34,6.86,7.75
              h11.09l-4.59,4.75h-6.68C9.71,32.93,3.19,29.44,2.99,21.13V0.01H0.05v37.3h35.87V17.62l-4.05,4.19L15.88,21.82z"/>
                    </g>
                  </svg>
                  {orchestratorSettings?.backend === "alcf" && (
                  <svg className="logo-svg" viewBox="87 0 26 24">
                    <path fill="#007934" d="M95.9 15.3h-8.1l4 7z"></path>
                    <path d="M103.9 15.3h-8.1l-4 7H108l-4.1-7z" fill="#0082ca"></path>
                    <path fill="#101e8e" d="M112 15.3h-8.1l4.1 7z"></path>
                    <path fill="#fff" d="M103.9 15.3h-8l4-7z"></path>
                    <path fill="#a22a2e" d="M103.9 1.3h-8l4 7z"></path>
                    <path fill="#d9272e" d="M103.9 1.3l-4 7 4 7h8.1z"></path>
                    <path d="M95.9 15.3l4-7-4-7-8.1 14h8.1z" fill="#82bc00"></path>
                  </svg>
                  )}
                </div>
            </div>

            <div className="app-header">
              <div className="app-header-content">
                <FlaskConical className="w-10 h-10 text-muted" />
                <h1 className="app-title">FLASK Copilot</h1>
              </div>
              <p className="app-subtitle">Real-time molecular assistant</p>
              <p className="app-subtitle">Connected to simulators at <code>llnl.gov</code> (LLNL)</p>
              <p className="app-subtitle">
                Connected to orchestrator <code>{orchestratorSettings.model || 'Not configured'}</code> at
                &nbsp;<code>{getDisplayUrl()}</code> ({orchestratorSettings.backendLabel})
              </p>
            </div>


            <div className="flex justify-end gap-2 mb-4">
              <div className="dropdown">
                <button onClick={(e) => { e.stopPropagation(); setLoadDropdownOpen(!loadDropdownOpen); }} onMouseDown={(e) => e.stopPropagation()} disabled={isComputing} className="btn btn-secondary btn-sm">
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
                  </svg>
                  Load
                  <svg className={`w-3 h-3 transition-transform ${loadDropdownOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                  </svg>
                </button>
                {loadDropdownOpen && (
                  <div className="dropdown-menu" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                    <button onClick={() => { loadContextFromFile(); setLoadDropdownOpen(false); }} className="dropdown-item">
                      <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
                      </svg>
                      Load from File
                    </button>
                    <div className="dropdown-divider"></div>
                    <button
                      onClick={() => {
                        if (window.confirm('Delete all projects and experiments from the database and clear cached data? This cannot be undone.')) {
                          // Clear session persistence first so no stale dbSessionId
                          // can trigger re-creation of deleted experiments.
                          sessionPersistence.clearSession();
                          restoredSessionRef.current = null;
                          projectData.clearAllProjects().then(() => {
                            projectSidebar.setSelection({ projectId: null, experimentId: null });
                          });
                        }
                        setLoadDropdownOpen(false);
                      }}
                      className="dropdown-item text-red-400"
                    >
                      <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
                      </svg>
                      Clear All Projects
                    </button>
                  </div>
                )}
              </div>
              <div className="dropdown">
                <button onClick={(e) => { e.stopPropagation(); setSaveDropdownOpen(!saveDropdownOpen); }} onMouseDown={(e) => e.stopPropagation()} className="btn btn-secondary btn-sm">
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 7H5a2 2 0 00-2 2v9a2 2 0 002 2h14a2 2 0 002-2V9a2 2 0 00-2-2h-3m-1 4l-3 3m0 0l-3-3m3 3V4" />
                  </svg>
                  Save
                  <svg className={`w-3 h-3 transition-transform ${saveDropdownOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                  </svg>
                </button>
                {saveDropdownOpen && (
                  <div className="dropdown-menu" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                    <div className="dropdown-item-static">
                      <label className="flex items-center gap-2 cursor-pointer">
                        <input 
                          type="checkbox" 
                          checked={autoCheckpointEnabled} 
                          onChange={(e) => setAutoCheckpointEnabled(e.target.checked)}
                          className="w-4 h-4 rounded border-purple-400/50 bg-white/20 text-purple-600 focus:ring-purple-500 focus:ring-offset-0"
                        />
                        <span className="text-sm">Auto-checkpoint on new molecule</span>
                      </label>
                    </div>
                    <div className="dropdown-divider"></div>
                    <button 
                      onClick={triggerManualCheckpoint} 
                      disabled={treeNodes.length === 0}
                      className="dropdown-item"
                    >
                      <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 7H5a2 2 0 00-2 2v9a2 2 0 002 2h14a2 2 0 002-2V9a2 2 0 00-2-2h-3m-1 4l-3 3m0 0l-3-3m3 3V4" />
                      </svg>
                      Save Checkpoint Now
                    </button>
                    <div className="dropdown-divider"></div>
                    <button onClick={saveTree} disabled={treeNodes.length === 0} className="dropdown-item">Save Tree Only</button>
                    <button onClick={requestSaveContext} disabled={treeNodes.length === 0 || !wsConnected} className="dropdown-item">Save Full Context</button>
                  </div>
                )}
              </div>
              <button onClick={() => setSidebarOpen(!sidebarOpen)} className="btn btn-secondary btn-sm">
                <Brain className="w-4 h-4" />
                Reasoning
              </button>
              <SettingsButton
                initialSettings={orchestratorSettings}
                onSettingsChange={handleSettingsUpdateConfirm}
                onServerAdded={refreshToolsList}
                onServerRemoved={refreshToolsList}
                username={username}
                httpServerUrl={HTTP_SERVER}
              />

              {/* WebSocket Status Indicator */}
                <div
                  className="ws-status-indicator group"
                  onClick={(e) => {
                    e.stopPropagation();
                    if (!wsConnected) {
                      reconnectWS();
                    } else {
                      setWsTooltipPinned(!wsTooltipPinned);
                    }
                  }}
                  onDoubleClick={(e) => {
                    e.stopPropagation();
                    reconnectWS();
                  }}
                  onMouseDown={(e) => e.stopPropagation()}
                  title={wsTooltipPinned ? "" : "Click for details • Double-click to reconnect"}
                >
                  <div className="relative cursor-pointer">
                    <div className={`status-indicator absolute ${
                      wsReconnecting ? 'status-indicator-ping bg-yellow-400' :
                      wsConnected ? '' :
                      ''
                    }`} />
                    <div className={`status-indicator ${
                      wsReconnecting ? 'status-indicator-reconnecting' :
                      wsConnected ? 'status-indicator-connected' :
                      'status-indicator-disconnected'
                    } ${wsTooltipPinned ? 'ring-2 ring-white ring-offset-2 ring-offset-slate-900' : ''}`} />
                  </div>

                  <div className={`ws-tooltip transition-opacity z-50 ${
                    wsTooltipPinned ? 'opacity-100' : 'opacity-0 group-hover:opacity-100 pointer-events-none'
                  }`}>
                    <div
                      onClick={(e) => e.stopPropagation()}
                      onMouseDown={(e) => e.stopPropagation()}
                    >
                      <div className="flex items-center justify-between mb-1">
                        <div className={`font-semibold ${
                          wsReconnecting ? 'status-reconnecting' :
                          wsConnected ? 'status-connected' :
                          'status-disconnected'
                        }`}>
                          {wsReconnecting ? '● Reconnecting...' :
                          wsConnected ? '● Connected' :
                          '● Disconnected'}
                          {wsConnected && username !== "<LOCAL USER>" && ` as ${username}`}
                        </div>
                        {wsTooltipPinned && (
                          <button
                            onClick={(e) => {
                              e.stopPropagation();
                              setWsTooltipPinned(false);
                            }}
                            className="text-muted hover:text-primary transition-colors cursor-pointer"
                          >
                            <X className="w-3 h-3" />
                          </button>
                        )}
                      </div>
                      <div className="text-secondary text-xs mt-1">
                        {WS_SERVER}
                        {wsError && (
                          <div className="text-tertiary text-xs mt-1">
                            {wsError}
                          </div>
                        )}
                        {wsConnected && availableTools.length > 0 && (
                          <div className="mt-3 pt-2 border-t border-secondary">
                            <div className="text-tertiary text-xs font-semibold mb-1.5">
                              Available Tool Servers ({availableTools.length})
                            </div>
                            <div className="custom-scrollbar max-h-60 overflow-y-auto pr-1">
                              {availableTools.map((tool, idx) => (
                                <div key={idx} className="text-xs bg-secondary rounded px-2 py-1 mb-1">
                                  <div className="text-secondary font-medium">{tool.server || "server" as string}</div>
                                  {tool.names && (
                                    <div className="text-tertiary mt-0.5 text-[10px] leading-tight">
                                      {tool.names.join(", ")}
                                    </div>
                                  )}
                                </div>
                              ))}
                            </div>
                          </div>
                        )}
                        {wsConnected && availableTools.length === 0 && (
                          <div className="mt-2 text-tertiary text-xs italic">
                            No MCP tool servers detected!
                          </div>
                        )}
                      </div>
                      {!wsConnected && !wsReconnecting && (
                        <div className="mt-2">
                          <div className="text-tertiary text-xs italic">
                            Backend server required for computation
                          </div>
                          {wsTooltipPinned && (
                            <button
                              onClick={(e) => {
                                e.stopPropagation();
                                reconnectWS();
                              }}
                              className="mt-2 w-full btn btn-secondary btn-sm"
                            >
                              <RefreshCw className="w-3 h-3" />
                              Reconnect
                            </button>
                          )}
                        </div>
                      )}
                      {!wsTooltipPinned && (
                        <div className="text-muted text-[10px] mt-2 italic text-center border-t border-secondary pt-1.5">
                          {wsConnected ?
                            "Click to pin, double-click to reconnect" :
                            "Click to reconnect"
                            }
                        </div>
                      )}
                    </div>
                  </div>
                </div>
            </div>

            <div className="card card-padding mb-6">
              <div className="input-row">
                <div className="flex-1">
                  <label className="form-label">Starting Molecule (SMILES)</label>
                  <input type="text" value={smiles} onChange={(e) => setSmiles(e.target.value)} disabled={isComputing} placeholder="Enter SMILES notation" className="form-input text-lg" />
                </div>
                <div>
                  <label className="form-label flex items-center gap-2 cursor-pointer">
                    <input type="checkbox" checked={autoZoom} onChange={(e) => setAutoZoom(e.target.checked)} className="form-checkbox" />
                    Auto-zoom to fit
                  </label>
                </div>
              </div>

              <div className="flex items-center gap-4">
                <div className="input-row-controls">
                  <div>
                    <label className="form-label">
                      Problem Type
                        {problemType === "custom" && (!systemPrompt || !problemPrompt) && (
                          <span className="warning-tooltip">
                            [!]
                          <div className="warning-tooltip-content">
                            <div className="warning-tooltip-box">
                              Custom problem description not given
                            </div>
                          </div>
                          </span>
                        )}
                    </label>
                    <select value={problemType} onChange={(e) => {reset(); resetProblemType(e.target.value)}} disabled={isComputing} className="form-select">
                      <option value="retrosynthesis">Retrosynthesis</option>
                      <option value="optimization">Lead Molecule Optimization</option>
                      <option value="custom">Custom</option>
                    </select>
                  </div>
                  { /* Problem-specific UI */ }
                  {problemType === "optimization" &&
                    <div>
                      <label className="form-label">
                        Property
                        {propertyType === "custom" && (!customPropertyName || !customPropertyDesc) && (
                          <span className="warning-tooltip">
                            [!]
                          <div className="warning-tooltip-content">
                            <div className="warning-tooltip-box">
                              Property name or description not given
                            </div>
                          </div>
                          </span>
                        )}
                      </label>
                      <select value={propertyType} onChange={(e) => {setPropertyType(e.target.value)}} disabled={isComputing} className="form-select w-48">
                        <option value="density">{PROPERTY_NAMES["density"]}</option>
                        <option value="hof">{PROPERTY_NAMES["hof"]}</option>
                        <option value="bandgap">{PROPERTY_NAMES["bandgap"]}</option>
                        <option value="custom">Other</option>
                      </select>
                    </div>
                  }
                  {problemType === "optimization" && propertyType === "custom" &&
                    <button onClick={() => setEditPropertyModal(true)} disabled={isComputing} className="btn btn-tertiary mt-5">
                      Property...
                    </button>
                  }
                  {problemType === "custom" &&
                    <button onClick={() => setEditPromptsModal(true)} disabled={isComputing || problemType !== "custom"} className="btn btn-tertiary mt-5">
                      Edit
                    </button>
                  }
                  <button
                    onClick={() => setShowCustomizationModal(true)}
                    disabled={isComputing}
                    className="btn btn-tertiary mt-5"
                  >
                    <Sliders className="w-4 h-4" />
                    Customize
                    {(selectedTools.length > 0 || (problemType === "optimization" && customization.enableConstraints)) && (
                      <span className="ml-1 px-2 py-0.5 text-xs rounded-full bg-blue-500/20 text-blue-400 flex items-center gap-1">
                        {selectedTools.length > 0 && (
                          <>
                            {selectedTools.length}
                            <Wrench className="w-3 h-3" />
                          </>
                        )}
                        {selectedTools.length > 0 && problemType === "optimization" && customization.enableConstraints && <span className="mx-0.5">•</span>}
                        {problemType === "optimization" && customization.enableConstraints && (
                          <>
                            ON
                            <Settings className="w-3 h-3" />
                          </>
                        )}
                      </span>
                    )}
                  </button>
                </div>

                <div className="input-row-actions">
                  <div className="relative group">
                    <div>
                      <label className="form-label">Actions</label>
                      <button
                        onClick={() => {
                          // If modal is minimized, reopen it
                          if (debugModalMinimized) {
                            handleReopenPromptModal();
                            return;
                          }
                          if (treeNodes.length > 0) {
                            if (
                              !window.confirm(
                                'Are you sure you want to rerun this computation? This will clear all previous progress.'
                              )
                            ) {
                              return;
                            }
                          }
                          runComputation();
                        }}
                        disabled={
                          !wsConnected ||
                          (isComputing && !debugModalMinimized) ||
                          (!smiles && !debugModalMinimized)
                        }
                        className="btn btn-primary"
                      >
                        {debugModalMinimized ? (
                          <>
                            <Bug className="w-5 h-5" />
                            Continue
                          </>
                        ) : isComputing ? (
                          <>
                            <Loader2 className="w-5 h-5 animate-spin" />
                            Computing
                          </>
                        ) : treeNodes.length === 0 ? (
                          <>
                            <Play className="w-5 h-5" />
                            Run
                          </>
                        ) : (
                          <>
                            <RefreshCw className="w-5 h-5" />
                            Rerun
                          </>
                        )}
                      </button>
                      {(!wsConnected ||
                        (isComputing && !debugModalMinimized) ||
                        (!smiles && !debugModalMinimized)) && (
                        <div className="tooltip absolute bottom-full mb-2 left-1/2 transform -translate-x-1/2 opacity-0 group-hover:opacity-100 transition-opacity pointer-events-none">
                          <div className="tooltip-content whitespace-nowrap">
                            {!wsConnected ? 'Backend server not connected' :
                            isComputing ? 'Computation already running' :
                            'Enter a SMILES string first'}
                          </div>
                          {!wsConnected && (
                            <div className="text-tertiary text-xs mt-1">
                              Start the backend server at {WS_SERVER}
                            </div>
                          )}
                        </div>
                      )}
                    </div>
                  </div>
                  <div>
                    <label className="form-label">&nbsp;</label>
                    <button onClick={() => {
                      const [updatedNodes, updatedEdges] = relayoutTree(treeNodes, edges);
                      setTreeNodes(updatedNodes);
                      setEdges(updatedEdges);
                    }} disabled={isComputing || treeNodes.length === 0} className="btn btn-secondary">
                      <Sparkles className="w-5 h-5" />Relayout
                    </button>
                  </div>
                  <div>
                    <label className="form-label">&nbsp;</label>
                    <button onClick={() => {
                      if (window.confirm('Are you sure you want to reset this window? This will clear all molecules.')) {
                        reset();
                      }
                    }} disabled={isComputing} className="btn btn-tertiary">
                      <RotateCcw className="w-5 h-5" />Reset
                    </button>
                  </div>
                  <div>
                    <label className="form-label">&nbsp;</label>
                    <button onClick={stop} disabled={!isComputing} className="btn btn-tertiary">
                      <X className="w-5 h-5" />Stop
                    </button>
                  </div>
                </div>
              </div>
            </div>

            <div className="card relative" style={{ height: '600px' }}>
              {treeNodes.length === 0 && !isComputing ? (
                <div className="empty-state">
                  <FlaskConical className="empty-state-icon" />
                  <p className="empty-state-text">
                    {wsConnected ?
                      `Click "Run" to start ${problemType === "optimization" ? "molecular discovery" : "the molecular computation tree"}` :
                      "Waiting for backend connection..."
                    }
                  </p>
                  <p className="empty-state-subtext">
                    {autoZoom ? 'Auto-zoom will fit all molecules' : 'Drag to pan • Scroll to zoom'}
                  </p>
                  {!wsConnected && (
                    <div className="alert alert-warning max-w-md mt-4">
                      <div className="alert-warning-text text-center">
                        <strong>Backend Required:</strong> Start your Python backend server at <code className="bg-black/30 px-2 py-1 rounded">{WS_SERVER}</code> to enable molecular computations.
                      </div>
                    </div>
                  )}
                </div>
              ) : (
                <MoleculeGraph
                  {...graphState}
                  nodes={treeNodes}
                  edges={edges}
                  autoZoom={autoZoom}
                  setAutoZoom={setAutoZoom}
                  ctx={contextMenu}
                  handleNodeClick={handleNodeClick}
                  handleReactionClick={handleReactionClick}
                  handleReactionCardClick={handleReactionCardClick}
                  selectedReactionNodeId={selectedReactionNode?.id}
                  reactionSidebarOpen={reactionSidebarOpen}
                  rdkitModule={rdkitModule}
                />
              )}

              {reactionSidebarOpen && selectedReactionNode?.reaction && (
                <ReactionAlternativesSidebar
                  isOpen={reactionSidebarOpen}
                  onClose={handleCloseReactionAlternativesSidebar}
                  productMolecule={selectedReactionNode.label} // Strip HTML
                  productSmiles={selectedReactionNode.smiles}
                  alternatives={stableAlternatives}
                  onSelectAlternative={handleSelectAlternative}
                  onComputeTemplates={handleComputeTemplates}
                  onComputeFlaskAI={handleComputeFlaskAI}
                  wsConnected={wsConnected}
                  isComputing={isComputing}
                  isComputingTemplates={isComputingTemplates}
                  templatesSearched={selectedReactionNode.reaction.templatesSearched}
                  rdkitModule={rdkitModule}
                />
              )}
            </div>

            {isComputing && (
              <div className="alert alert-info">
                <div className="alert-info-text">
                  <Loader2 className="w-5 h-5 animate-spin" />
                  <span className="font-medium">Streaming molecules... {treeNodes.length} nodes discovered</span>
                </div>
              </div>
            )}

            {!isComputing && treeNodes.length > 0 && (
              <div className="alert alert-success">
                <div className="alert-success-text">
                  <span className="font-medium">
                    Computation complete! Generated {treeNodes.length} molecules
                  </span>
                </div>
              </div>
            )}

            {/* Metrics Dashboard */}
            {(problemType === "optimization") && (treeNodes.length > 0) && (
              <MetricsDashboard {...metricsDashboardState} treeNodes={treeNodes} />
            )}


          <div className="app-footer">
            <p>This work was performed under the auspices of the U.S. Department of Energy
            by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344
            (LLNL-CODE-2006345).</p>
            {VERSION && (
              <p>Server version: {VERSION}</p>
            )}
          </div>
          </div>
        </div>

        <ReasoningSidebar
          {...sidebarState}
          setSidebarOpen={setSidebarOpen}
          rdkitModule={rdkitModule}
          isOpen={sidebarOpen}
          onToggle={() => setSidebarOpen(!sidebarOpen)}
        />
      </div>

      {contextMenu && contextMenu.node && (
        <div className="context-menu" style={{ left: `${contextMenu.x + 10}px`, top: `${contextMenu.y + 10}px` }} onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
          <div className="context-menu-header">
            <div className="context-menu-label">Actions for {(contextMenu.isReaction ? "Reaction Resulting in" : "Molecule")}</div>
            <div className="context-menu-title" dangerouslySetInnerHTML={{__html: contextMenu.node.label}}></div>
          </div>

          { !contextMenu.isReaction && (
            <>
              <button
                onClick={() => copyToClipboard(contextMenu.node!.smiles, 'smiles', setCopiedField)}
                className="context-menu-item"
              >
                {copiedField === 'smiles' ? (
                  <>✓ Copied!</>
                ) : (
                  <>
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 16H6a2 2 0 01-2-2V6a2 2 0 012-2h8a2 2 0 012 2v2m-6 12h8a2 2 0 002-2v-8a2 2 0 00-2-2h-8a2 2 0 00-2 2v8a2 2 0 002 2z" />
                    </svg>
                    Copy SMILES
                  </>
                )}
              </button>
              <button onClick={() => handleCustomQuery(contextMenu.node!, "query-molecule")} className="context-menu-item">
                <MessageCircleQuestion className="w-4 h-4" /> Ask about molecule...
              </button>
            </>
          )}

          { contextMenu.isReaction && (
            <button
              onClick={() => copyToClipboard(contextMenu.node!.reaction!.hoverInfo, 'reaction', setCopiedField)}
              className="context-menu-item"
            >
              {copiedField === 'reaction' ? (
                <>✓ Copied!</>
              ) : (
                <>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 16H6a2 2 0 01-2-2V6a2 2 0 012-2h8a2 2 0 012 2v2m-6 12h8a2 2 0 002-2v-8a2 2 0 00-2-2h-8a2 2 0 00-2 2v8a2 2 0 002 2z" />
                  </svg>
                  Copy details
                </>
              )}
            </button>
          )}

          {
            /* Lead Molecule Optimization */ problemType === 'optimization' && (
              <>
                <button
                  onClick={() => {
                    const nodeId = contextMenu.node!.id;
                    // Delete all nodes from this point on
                    setTreeNodes((prev) => {
                      return prev.filter((n) => n.x <= contextMenu.node!.x);
                    });
                    sendMessageToServer('optimize-from', {
                      nodeId: nodeId,
                      propertyType,
                      customPropertyName,
                      customPropertyDesc,
                      customPropertyAscending,
                      smiles: contextMenu.node!.smiles,
                      xpos: contextMenu.node!.x,
                    });
                    setIsComputing(true);
                  }}
                  className="context-menu-item context-menu-divider"
                >
                  <StepForward className="w-4 h-4" />
                  Refine search from here
                </button>
                <button
                  onClick={() => {
                    // const nodeId = contextMenu.node!.id;
                    // Delete all nodes from this point on
                    setTreeNodes((prev) => {
                      return prev.filter((n) => n.x <= contextMenu.node!.x);
                    });
                    handleCustomQuery(contextMenu.node!, 'optimize-from');
                  }}
                  className="context-menu-item"
                >
                  <MessageSquareShare className="w-4 h-4" />
                  Refine search (with prompt)
                </button>
                <button
                  onClick={() => {
                    createNewRetrosynthesisExperiment(contextMenu.node!.smiles);
                  }}
                  className="context-menu-item"
                >
                  <FlaskConical className="w-4 h-4" />
                  Plan synthesis pathway
                </button>
              </>
            )
          }

          {/* Retrosynthesis (Molecule) */ (problemType == "retrosynthesis" && !contextMenu.isReaction) && (
            <>
              {!contextMenu.node.reaction && (
                <button onClick={() => {sendMessageToServer("compute-reaction-from", {nodeId: contextMenu.node!.id});}} className="context-menu-item context-menu-divider">
                  <TestTubeDiagonal className="w-4 h-4" />How do I make this?
                </button>
              )}
              {!isRootNode(contextMenu.node.id, treeNodes) && (
                <button onClick={() => {sendMessageToServer("recompute-parent-reaction", {nodeId: contextMenu.node!.id});}} className="context-menu-item context-menu-divider">
                  <Network className="w-4 h-4" />Substitute Molecule
                </button>
              )}
            </>
          )}

          {/* Retrosynthesis (Reaction) */ (problemType == "retrosynthesis" && contextMenu.isReaction) && (
            <>
              <button onClick={() => handleCustomQuery(contextMenu.node!, "query-reaction")} className="context-menu-item">
                <MessageCircleQuestion className="w-4 h-4" /> Ask about reaction...
              </button>
              <button onClick={() => {handleReactionCardClick(contextMenu.node!); setContextMenu({node: null, isReaction: false, x: 0, y: 0});}} className="context-menu-item context-menu-divider">
                <PanelRightOpen className="w-4 h-4" />Other Reactions...
              </button>
            </>
          )}

          { contextMenu.isReaction && (
              <div className="context-menu-details custom-scrollbar">
                <MarkdownText text={contextMenu.node!.reaction!.hoverInfo} />
              </div>
          )}
        </div>
      )}

      {customQueryModal && (
        <div className="modal-overlay">
          <div className="modal-content modal-content-lg">
            <div className="modal-header">
              <div>
                <h2 className="modal-title">Custom Query</h2>
                <div className="modal-subtitle" dangerouslySetInnerHTML={{__html: "for " + customQueryModal.label}}></div>
              </div>
              <button onClick={() => setCustomQueryModal(null)} className="btn-icon">
                <X className="w-6 h-6" />
              </button>
            </div>

            <textarea
              value={customQueryText}
              onChange={(e) => setCustomQueryText(e.target.value)}
              placeholder="Enter your custom query here..."
              className="form-textarea h-40"
            />

            <div className="modal-footer">
              <button onClick={submitCustomQuery} disabled={!customQueryText.trim()} className="btn btn-primary flex-1">
                <Send className="w-5 h-5" />
                Submit Query
              </button>

              <button onClick={() => setCustomQueryModal(null)} className="btn btn-tertiary">
                Cancel
              </button>
            </div>
          </div>
        </div>
      )}

      {editPromptsModal && (
        <div className="modal-overlay">
          <div className="modal-content modal-content-lg">
            <div className="modal-header">
              <div>
                <h2 className="modal-title">Edit Prompts</h2>
                <p className="modal-subtitle">Configure system and problem-specific prompts</p>
              </div>
              <button onClick={() => setEditPromptsModal(false)} className="btn-icon">
                <X className="w-6 h-6" />
              </button>
            </div>

            <div className="modal-body space-y-4">
              <div className="form-group">
                <label className="form-label">System Prompt</label>
                <textarea value={systemPrompt} onChange={(e) => setSystemPrompt(e.target.value)} placeholder="Enter system-level instructions..." className="form-textarea h-32" />
              </div>

              <div className="form-group">
                <label className="form-label">Problem Prompt</label>
                <textarea value={problemPrompt} onChange={(e) => setProblemPrompt(e.target.value)} placeholder="Enter problem-specific instructions..." className="form-textarea h-32" />
              </div>
            </div>

            <div className="modal-footer">
              <button onClick={() => { savePrompts(DEFAULT_CUSTOM_SYSTEM_PROMPT, ''); }} className="btn btn-tertiary">
                <RotateCcw className="w-4 h-4" />
                Reset
              </button>
              <button onClick={() => {savePrompts(systemPrompt, problemPrompt); setProblemType('custom');}} className="btn btn-primary flex-1">
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                </svg>
                Save Prompts
              </button>
            </div>
          </div>
        </div>
      )}

      {editPropertyModal && (
        <div className="modal-overlay">
          <div className="modal-content modal-content-lg">
            <div className="modal-header">
              <div>
                <h2 className="modal-title">Custom Property</h2>
                <p className="modal-subtitle">Configure custom property type and units</p>
              </div>
              <button onClick={() => setEditPropertyModal(false)} className="btn-icon">
                <X className="w-6 h-6" />
              </button>
            </div>
            <div className="modal-body space-y-4">
              <div className="form-group">
                <label className="form-label">Property Name</label>
                <input
                    type="text"
                    value={customPropertyName}
                    onChange={(e) => setCustomPropertyName(e.target.value)}
                    placeholder="Enter a name for the property"
                    className="form-input form-input-text"
                  />
              </div>
              <div className="form-group">
                <label className="form-label">Property Description</label>
                <textarea value={customPropertyDesc} onChange={(e) => setCustomPropertyDesc(e.target.value)} placeholder="Enter a description of the property and its units..." className="form-textarea h-32" />
              </div>
            </div>
            <div className="flex-center gap-md py-4">
              <span className="text-sm text-secondary">Higher is better</span>
              <button
                onClick={() => setCustomPropertyAscending(!customPropertyAscending)}
                className={`toggle-switch ${customPropertyAscending ? 'toggle-switch-off' : 'toggle-switch-on'}`}
              >
                <div
                  className={`toggle-switch-handle ${customPropertyAscending ? 'toggle-switch-handle-off' : 'toggle-switch-handle-on'}`}
                />
              </button>
              <span className="text-sm text-secondary">Lower is better</span>
            </div>
            <div className="modal-footer">
              <button onClick={() => { saveCustomProperty('', '', true); }} className="btn btn-tertiary">
                <RotateCcw className="w-4 h-4" />
                Reset
              </button>
              <button onClick={() => {saveCustomProperty(customPropertyName, customPropertyDesc, customPropertyAscending);}} className="btn btn-primary flex-1">
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                </svg>
                Save
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Combined Customization Modal */}
      <CombinedCustomizationModal
        isOpen={showCustomizationModal}
        onClose={() => setShowCustomizationModal(false)}
        availableToolsMap={availableToolsMap}
        selectedTools={selectedTools}
        onToolSelectionChange={setSelectedTools}
        onToolConfirm={handleToolSelectionConfirm}
        initialCustomization={customization}
        onCustomizationSave={setCustomization}
        showOptimizationTab={problemType === "optimization"}
      />

      {/* Prompt Debugging Modal */}
      <Modal
        isOpen={!!promptBreakpoint && !debugModalMinimized}
        onClose={handleMinimizePromptModal}
        closeIcon={<Minus className="w-6 h-6" />}
        title="🔍 AI Prompt Breakpoint"
        subtitle="Review and modify the prompt before sending to the AI"
        size="lg"
        footer={
          <>
            <button
              onClick={() => handlePromptBreakpointResponse()}
              className="btn btn-primary flex-1"
            >
              <CheckCircle className="w-5 h-5" />
              Approve & Send
            </button>
          </>
        }
      >
        <div className="space-y-4">
          {promptBreakpoint?.metadata && (
            <div className="glass-panel">
              <div className="text-sm font-semibold text-secondary mb-2">Context Information:</div>
              <pre className="text-xs text-tertiary overflow-x-auto">
                {JSON.stringify(promptBreakpoint.metadata, null, 2)}
              </pre>
            </div>
          )}

          <div className="form-group">
            <label className="form-label-block">
              Prompt Content
              <span className="text-xs text-tertiary ml-2">(Edit as needed before approving)</span>
            </label>
            <textarea
              value={editedPrompt}
              onChange={(e) => setEditedPrompt(e.target.value)}
              className="form-textarea"
              style={{ minHeight: '300px', fontFamily: 'monospace' }}
              placeholder="Prompt content will appear here..."
            />
          </div>

          <div className="alert alert-info">
            <div className="text-sm text-secondary">
              You can review and modify this prompt before it is sent to the AI model. Changes will
              be used for this request only.
            </div>
          </div>
        </div>
      </Modal>
    </div>
  );
};

export default ChemistryTool;
