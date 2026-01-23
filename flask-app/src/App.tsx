/* eslint-disable react-hooks/exhaustive-deps */
import React, { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { Loader2, FlaskConical, TestTubeDiagonal, Network, Play, RotateCcw, X, Send, RefreshCw, Sparkles, MessageCircleQuestion, StepForward, MessageSquareShare, Brain, PanelRightOpen } from 'lucide-react';
import 'recharts';
import 'react-markdown';
import 'remark-gfm';
import 'react-syntax-highlighter';
import 'react-syntax-highlighter/dist/esm/styles/prism';

import { WS_SERVER, VERSION } from './config';
import { DEFAULT_CUSTOM_SYSTEM_PROMPT, PROPERTY_NAMES } from './constants';
import { TreeNode, Edge, ContextMenuState, SidebarMessage, Tool, WebSocketMessageToServer, WebSocketMessage, SelectableTool, Experiment, OrchestratorSettings, ReactionAlternative } from './types';

import { loadRDKit } from './components/molecule';
import { ReasoningSidebar, useSidebarState } from './components/sidebar';
import { MoleculeGraph, useGraphState } from './components/graph';
import { MultiSelectToolModal } from './components/multi_select_tools';
import { ProjectSidebar, useProjectSidebar, useProjectManagement } from './components/project_sidebar';
import { SettingsButton, BACKEND_OPTIONS } from './components/settings_button';

import { findAllDescendants, hasDescendants, isRootNode, relayoutTree } from './tree_utils';
import { copyToClipboard } from './utils';

import './animations.css';
import { MetricsDashboard, useMetricsDashboardState } from './components/metrics';
import { useProjectData } from './hooks/useProjectData';
import { MarkdownText } from './components/markdown';
import { ReactionAlternativesSidebar } from './components/reaction_alternatives';


const ChemistryTool: React.FC = () => {
  const [smiles, setSmiles] = useState<string>('');
  const [problemType, setProblemType] = useState<string>('retrosynthesis');
  const [propertyType, setPropertyType] = useState<string>('density');
  const [systemPrompt, setSystemPrompt] = useState<string>(DEFAULT_CUSTOM_SYSTEM_PROMPT);
  const [problemPrompt, setProblemPrompt] = useState<string>('');
  const [editPromptsModal, setEditPromptsModal] = useState<boolean>(false);
  const [editPropertyModal, setEditPropertyModal] = useState<boolean>(false);
  const [customPropertyName, setCustomPropertyName] = useState<string>('');
  const [customPropertyDesc, setCustomPropertyDesc] = useState<string>('');
  const [customPropertyAscending, setCustomPropertyAscending] = useState<boolean>(true);
  const [isComputing, setIsComputing] = useState<boolean>(false);
  const [autoZoom, setAutoZoom] = useState<boolean>(true);
  const [treeNodes, setTreeNodes] = useState<TreeNode[]>([]);
  const [edges, setEdges] = useState<Edge[]>([]);
  const [contextMenu, setContextMenu] = useState<ContextMenuState>({node: null, isReaction: false, x: 0, y: 0});
  const [customQueryModal, setCustomQueryModal] = useState<TreeNode | null>(null);
  const [customQueryText, setCustomQueryText] = useState<string>('');
  const [customQueryType, setCustomQueryType] = useState<string | null>(null);
  const [wsConnected, setWsConnected] = useState<boolean>(false);
  const [saveDropdownOpen, setSaveDropdownOpen] = useState<boolean>(false);
  const [wsError, setWsError] = useState<string>('');
  const [wsReconnecting, setWsReconnecting] = useState<boolean>(false);
  const rdkitModule = loadRDKit();
  const [sidebarOpen, setSidebarOpen] = useState<boolean>(false);
  const [copiedField, setCopiedField] = useState<string | null>(null);
  const [availableTools, setAvailableTools] = useState<Tool[]>([]);
  const [wsTooltipPinned, setWsTooltipPinned] = useState<boolean>(false);
  const [username, setUsername] = useState<string>('<LOCAL USER>');

  const wsRef = useRef<WebSocket | null>(null);

  // Function to refresh tools list from backend
  const refreshToolsList = useCallback(() => {
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      console.log('🔄 Refreshing tools list from backend');
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

  const treeNodesRef = useRef(treeNodes);
  const edgesRef = useRef(edges);
  const sidebarStateRef = useRef(sidebarState);

  const [showToolSelectionModal, setShowToolSelectionModal] = useState<boolean>(false);
  const [selectedTools, setSelectedTools] = useState<number[]>([]);
  const [availableToolsMap, setAvailableToolsMap] = useState<SelectableTool[]>([]);

  // Reaction alternatives sidebar state
  const [reactionSidebarOpen, setReactionSidebarOpen] = useState<boolean>(false);
  const [selectedReactionNode, setSelectedReactionNode] = useState<TreeNode | null>(null);
  const [isComputingTemplates, setIsComputingTemplates] = useState<boolean>(false);

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
      console.log('🔄 Refreshing tools list after settings update');
      refreshToolsList();
     }
  };

  useEffect(() => {
    const handleClickOutside = (): void => {
      setContextMenu({node: null, isReaction: false, x:0, y:0});
      setSaveDropdownOpen(false);
      sidebarState.setSourceFilterOpen(false);
      setWsTooltipPinned(false);
      setCopiedField(null);
    };
    if (contextMenu || saveDropdownOpen || sidebarState.sourceFilterOpen || wsTooltipPinned) {
      window.addEventListener('mousedown', handleClickOutside);
      return () => window.removeEventListener('mousedown', handleClickOutside);
    }
  }, [contextMenu, saveDropdownOpen, sidebarState, wsTooltipPinned, projectSidebar]);


  // State management
  const getContext = (): Experiment => {
    return getContextRef.current();
  }

  const loadContextFromExperiment = (projectId: string, experimentId: string | null): void => {
    console.log('Loading context:', { projectId, experimentId });
    const project = projectData.projectsRef.current.find(p => p.id === projectId);
    if (project) {
      const experiment = project.experiments.find(e => e.id === experimentId);
      if (experiment) {
        loadContext(experiment);
      }
    }
    return;
  }

  const loadStateFromCurrentExperiment = (): void => {
    const { projectId, experimentId } = projectSidebar.selectionRef.current;
    if (projectId && experimentId) {
      loadContextFromExperiment(projectId, experimentId);
    }
  };

  const loadContext = (data: Experiment): void => {
    // Conditionally set everything that is in the context
    (data.smiles !== undefined) && setSmiles(data.smiles);
    (data.problemType !== undefined) && setProblemType(data.problemType);
    (data.systemPrompt !== undefined) && setSystemPrompt(data.systemPrompt);
    (data.problemPrompt !== undefined) && setProblemPrompt(data.problemPrompt);
    (data.propertyType !== undefined) && setPropertyType(data.propertyType);
    (data.customPropertyName !== undefined) && setCustomPropertyName(data.customPropertyName);
    (data.customPropertyDesc !== undefined) && setCustomPropertyDesc(data.customPropertyDesc);
    (data.customPropertyAscending !== undefined) && setCustomPropertyAscending(data.customPropertyAscending);
    data.treeNodes && setTreeNodes(data.treeNodes);
    data.edges && setEdges(data.edges);
    data.metricsHistory && metricsDashboardState.setMetricsHistory(data.metricsHistory);
    data.visibleMetrics && metricsDashboardState.setVisibleMetrics(data.visibleMetrics);
    if (data.graphState) {
      graphState.setZoom(data.graphState.zoom);
      graphState.setOffset(data.graphState.offset);
    }
    (data.autoZoom !== undefined) && setAutoZoom(data.autoZoom);
    if (data.sidebarState) {
      sidebarState.setMessages(data.sidebarState.messages);
      sidebarState.setVisibleSources(data.sidebarState.visibleSources);
    }
    data.experimentContext && sendMessageToServer('load-context', {experimentContext: data.experimentContext});
  }

  const saveStateToExperiment = useCallback((): boolean => {
    // Use the ref directly to always get the latest selection
    const projectId = projectSidebar.selectionRef.current.projectId;
    const experimentId = projectSidebar.selectionRef.current.experimentId;
    console.log("Saving experiments", projectId, experimentId);
    if (projectId && experimentId) {
      projectManagement.updateExperiment(projectId, getContext());
      return true;
    }
    return false;
  }, [projectSidebar.selectionRef, projectManagement, getContext]);

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
      userPrompt: problemPrompt
    };

    wsRef.current?.send(JSON.stringify(message));
  };

  const reconnectingRef = useRef(false);

  const reconnectWS = (): void => {
    if (reconnectingRef.current) return; // Prevent overlapping reconnects
    reconnectingRef.current = true;

    if (wsRef.current &&
        (wsRef.current.readyState === WebSocket.OPEN ||
         wsRef.current.readyState === WebSocket.CONNECTING)) {
      wsRef.current.close();
    }

    setWsReconnecting(true);

    const socket = new WebSocket(WS_SERVER);
    wsRef.current = socket;

    socket.onopen = () => {
      reconnectingRef.current = false; // Clear guard
      console.log('WebSocket connected');
      setWsConnected(true);
      setWsReconnecting(false);
      setWsError('');
      reset();  // Server state must match UI state

      loadStateFromCurrentExperiment();

      socket.send(JSON.stringify({ action: 'list-tools' }));
      socket.send(JSON.stringify({ action: 'get-username' }));
    };

    socket.onmessage = (event: MessageEvent) => {
      if (wsRef.current !== socket) return; // Ignore messages from old sockets

      const data: WebSocketMessage = JSON.parse(event.data);

      if (data.type === 'node') {
        setTreeNodes(prev => [...prev, data.node!]);
      } else if (data.type === 'node_update') {
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

        /*
        // When a node update includes reaction info, preserve alternatives and templatesSearched:
        setTreeNodes(prev => prev.map(n => {
          if (n.id === data.node!.id) {
            // Preserve existing alternatives and templatesSearched if not explicitly updated
            if (restData.reaction && n.reaction) {
              if (!restData.reaction.alternatives) {
                restData.reaction.alternatives = n.reaction.alternatives;
              }
              if (restData.reaction.templatesSearched === undefined) {
                restData.reaction.templatesSearched = n.reaction.templatesSearched;
              }
            }
            return { ...n, ...restData };
          }
          return n;
        }));
        */
      } else if (data.type === 'node_delete') {
        setTreeNodes(prev => {
          const descendants = findAllDescendants(data.node!.id, prev);
          return prev.filter(n => !descendants.has(n.id) && n.id !== data.node!.id);
        });
        setEdges(prev => prev.filter(e =>
          e.fromNode !== data.node!.id && e.toNode !== data.node!.id
        ));
      } else if (data.type === 'subtree_update') {
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
        setEdges(prev => [...prev, data.edge!]);
      } else if (data.type === 'edge_update') {
        const { id, ...restData } = data.edge!;
        setEdges(prev => prev.map(e =>
          e.id === data.edge!.id ? { ...e, ...restData } : e
         ));
      } else if (data.type === 'subtree_delete') {
        let descendantsSet: Set<string>;
        setTreeNodes(prev => {
          descendantsSet = findAllDescendants(data.node!.id, prev);
          return prev.filter(n => !descendantsSet.has(n.id));
        });
        setEdges(prev => prev.filter(e =>
          !descendantsSet!.has(e.fromNode) && !descendantsSet!.has(e.toNode)
        ));
      } else if (data.type === 'complete') {
        setIsComputing(false);
        setIsComputingTemplates(false);
        unhighlightNodes();
        saveStateToExperiment();  // Keep experiment up to date
      } else if (data.type === 'response') {
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
      } else if (data.type === 'reaction-alternatives-response') {
        // Backend sends template alternatives for a specific node
        setTreeNodes(prev => prev.map(n =>
          n.id === data.node!.id && n.reaction
            ? {
                ...n,
                reaction: {
                  ...n.reaction,
                  alternatives: data.alternatives ?? [],
                  templatesSearched: true  // Mark that templates have been searched
                }
              }
            : n
        ));
        setIsComputingTemplates(false);
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
      }
    };
  };

  // Connect WebSocket on mount
  useEffect(() => {
    reconnectWS();

    return () => {
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    };
  }, []);

  const reset = (): void => {
    setTreeNodes([]);
    setEdges([]);
    setIsComputing(false);
    graphState.setOffset({ x: 50, y: 50 });
    graphState.setZoom(1);
    setContextMenu({node: null, isReaction: false, x:0, y:0});
    setCustomQueryModal(null);
    metricsDashboardState.setMetricsHistory([]);
    sidebarState.setMessages([]);
    setSaveDropdownOpen(false);
    sidebarState.setSourceFilterOpen(false);
    setWsTooltipPinned(false);
    // Clear and close reaction alternatives sidebar
    setReactionSidebarOpen(false);
    setSelectedReactionNode(null);
    setIsComputingTemplates(false);
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      wsRef.current.send(JSON.stringify({ action: 'reset' }));
    }
  };

  const unhighlightNodes = (): void => {
    setTreeNodes(prev => prev.map(n =>
      n.highlight === "yellow" ? { ...n, highlight: "normal" } : n
    ));
  };

  const stop = (): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      // alert('WebSocket not connected');
      return;
    }

    console.log('Sending stop command to server');
    setIsComputing(false);
    setIsComputingTemplates(false);
    wsRef.current.send(JSON.stringify({ action: 'stop' }));
    unhighlightNodes();
    saveStateToExperiment();
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
      systemPrompt, problemPrompt, propertyType, customPropertyName, customPropertyDesc, customPropertyAscending, projectData, projectSidebar]);



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
                   customPropertyDesc, customPropertyAscending, treeNodes, edges, graphState, metricsDashboardState, sidebarState, experimentContext };
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

    // Delete subtree locally
    setTreeNodes(prev => {
      const descendants = findAllDescendants(nodeId, prev);
      return prev.filter(n => !descendants.has(n.id) && n.id !== nodeId);
    });
    setEdges(prev => prev.filter(e =>
      e.fromNode !== nodeId && e.toNode !== nodeId
    ));

    // Tell backend that the alternative subtree has been chosen
    wsRef.current.send(JSON.stringify({
      action: "set-reaction-alternative",
      nodeId: nodeId,
      alternativeId: alt.id
    }));

    setReactionSidebarOpen(false);
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
  }, [selectedReactionNode?.id, selectedReactionNode?.reaction?.alternatives?.length]);


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
        xpos: customQueryModal?.x
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
          onReset={reset}
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
              <button onClick={loadContextFromFile} disabled={isComputing} className="btn btn-secondary btn-sm">
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
                </svg>
                Load
              </button>
              <div className="dropdown">
                <button onClick={(e) => { e.stopPropagation(); setSaveDropdownOpen(!saveDropdownOpen); }} onMouseDown={(e) => e.stopPropagation()} disabled={isComputing || treeNodes.length === 0 || !wsConnected} className="btn btn-secondary btn-sm">
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
                    <button onClick={saveTree} className="dropdown-item">Save Tree Only</button>
                    <button onClick={requestSaveContext} className="dropdown-item dropdown-divider">Save Full Context</button>
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
                            ⚠️
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
                            ⚠️
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
                    onClick={() => setShowToolSelectionModal(true)}
                    disabled={isComputing}
                    className="btn btn-tertiary mt-5"
                  >
                    Select Tools {selectedTools.length > 0 && `(${selectedTools.length})`}
                  </button>
                </div>

                <div className="input-row-actions">
                  <div className="relative group">
                    <div>
                      <label className="form-label">
                        Actions
                      </label>
                      <button onClick={() => {
                        if (treeNodes.length > 0) {
                          if (!window.confirm('Are you sure you want to rerun this computation? This will clear all previous progress.')) {
                            return;
                          }
                        }
                        runComputation();
                      }} disabled={!wsConnected || isComputing || !smiles} className="btn btn-primary">
                        {isComputing ? <><Loader2 className="w-5 h-5 animate-spin" />Computing</> : (treeNodes.length === 0 ? <><Play className="w-5 h-5" />Run</> : <><RefreshCw className="w-5 h-5" />Rerun</>)}
                      </button>
                      {(!wsConnected || isComputing || !smiles) && (
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

          {/* Lead Molecule Optimization */ (problemType === "optimization") && (
            <>
            <button onClick={() => {
              const nodeId = contextMenu.node!.id;
                // Delete all nodes from this point on
                setTreeNodes(prev => {
                  return prev.filter(n => n.x <= contextMenu.node!.x);
                });
                sendMessageToServer("optimize-from", {nodeId: nodeId, propertyType, customPropertyName, customPropertyDesc, customPropertyAscending, smiles: contextMenu.node!.smiles, xpos: contextMenu.node!.x});
                setIsComputing(true);
              }}  className="context-menu-item context-menu-divider">
              <StepForward className="w-4 h-4" />
              Refine search from here
            </button>
            <button onClick={() => {
                const nodeId = contextMenu.node!.id;
                // Delete all nodes from this point on
                setTreeNodes(prev => {
                  return prev.filter(n => n.x <= contextMenu.node!.x);
                });
                handleCustomQuery(contextMenu.node!, "optimize-from");
              }}  className="context-menu-item">
              <MessageSquareShare className="w-4 h-4" />
              Refine search (with prompt)
            </button>
            <button onClick={() => { createNewRetrosynthesisExperiment(contextMenu.node!.smiles); }}  className="context-menu-item">
              <FlaskConical className="w-4 h-4" />
              Plan synthesis pathway
            </button>
            </>
          )}

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

      {/* Use the MultiSelectToolModal component */}
      <MultiSelectToolModal
        isOpen={showToolSelectionModal}
        onClose={() => setShowToolSelectionModal(false)}
        availableToolsMap={availableToolsMap}
        selectedTools={selectedTools}
        onSelectionChange={setSelectedTools}
        onConfirm={handleToolSelectionConfirm}
        title="Select Tools to use for Task" // Optional, defaults to "Select Tools"
      />

    </div>
  );
};

export default ChemistryTool;
