/* eslint-disable react-hooks/exhaustive-deps */
import React, { useState, useRef, useEffect, useCallback } from 'react';
import { Loader2, FlaskConical, TestTubeDiagonal, Network, Play, RotateCcw, X, Send, RefreshCw, Sparkles, MessageCircleQuestion, StepForward, MessageSquareShare } from 'lucide-react';
import 'recharts';

import { WS_SERVER, VERSION } from './config';
import { TreeNode, Edge, ContextMenuState, SidebarMessage, Tool, WebSocketMessageToServer, WebSocketMessage, SelectableTool, Experiment, ProfileSettings } from './types';

import { loadRDKit } from './components/molecule';
import { ReasoningSidebar, useSidebarState } from './components/sidebar';
import { MoleculeGraph, useGraphState } from './components/graph';
import { MultiSelectToolModal } from './components/multi_select_tools';
import { ProjectSidebar, useProjectSidebar, useProjectManagement } from './components/project_sidebar';
import { ProfileButton } from './components/profile_button';

import { findAllDescendants, hasDescendants, isRootNode, relayoutTree } from './tree_utils';
import { copyToClipboard } from './utils';

import './animations.css';
import { MetricsDashboard, useMetricsDashboardState } from './components/metrics';
import { useProjectData } from './hooks/useProjectData';


const ChemistryTool: React.FC = () => {
  const [smiles, setSmiles] = useState<string>('');
  const [problemType, setProblemType] = useState<string>('retrosynthesis');
  const [propertyType, setPropertyType] = useState<string>('density');
  const [systemPrompt, setSystemPrompt] = useState<string>('');
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
  const [contextMenu, setContextMenu] = useState<ContextMenuState>({node: null, x: 0, y: 0});
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
    const saved = localStorage.getItem('profileSettings');
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
      apiKey: ''
    };
  };
  const [profileSettings, setProfileSettings] = useState(getInitialSettings());

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
    // await fetch('/api/save-selection', { method: 'POST', body: JSON.stringify(payload) });
  };

  // Callback function to send updated profile to backend
  const handleProfileUpdateConfirm = async (
    settings: ProfileSettings,
  ): Promise<void> => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    console.log(`Updated Profile Saved`);

    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      const message = {
        action: 'update-profile-settings',
        backend: settings.backend,
        customUrl: settings.customUrl,
        model: settings.model,
        apiKey: settings.apiKey
      };
      wsRef.current.send(JSON.stringify(message));
     }
  };

  useEffect(() => {
    const handleClickOutside = (): void => {
      setContextMenu({node: null, x:0, y:0});
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
      const experimentName = `Experiment 1`;

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
      const experimentCount = project ? project.experiments.length + 1 : 1;
      const experimentName = `Experiment ${experimentCount}`;

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
        }
        setTreeNodes(prev => prev.map(n =>
          n.id === data.node!.id ? { ...n, ...restData } : n
        ));
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
        saveStateToExperiment();  // Keep experiment up to date
      } else if (data.type === 'response') {
        addSidebarMessage(data.message!);
        console.log('Server response:', data.message);
      } else if (data.type === 'available-tools-response') {
        setAvailableTools(data.tools || []);
        availableToolsMap.splice(0, availableToolsMap.length);
        setSelectedTools([]);
        availableTools.forEach((server: Tool, index: number, array: Tool[]) => {
          console.log(`Tool Server Element at index ${index}: ${server.server}`);
          availableToolsMap.push({id: index, tool_server: server})
        });
      } else if (data.type === 'update-orchestrator-profile') {
        // Handle profile settings updates from server
        const newSettings: ProfileSettings = {
          backend: data.profileSettings!.backend,
          useCustomUrl: data.profileSettings!.useCustomUrl,
          customUrl: data.profileSettings!.customUrl,
          model: data.profileSettings!.model,
          // Don't take the use custom model field from the backend
          // Check the model against the list of models in copilot
          // useCustomModel: data.profileSettings.useCustomModel,
          apiKey: data.profileSettings!.apiKey,
        };
        setProfileSettings(newSettings);
        console.log('Updating the profile settings ', newSettings);
        localStorage.setItem('profileSettings', JSON.stringify(newSettings));
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
    setContextMenu({node: null, x:0, y:0});
    setCustomQueryModal(null);
    metricsDashboardState.setMetricsHistory([]);
    sidebarState.setMessages([]);
    setSaveDropdownOpen(false);
    sidebarState.setSourceFilterOpen(false);
    setWsTooltipPinned(false);
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      wsRef.current.send(JSON.stringify({ action: 'reset' }));
    }
  };

  const stop = (): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      // alert('WebSocket not connected');
      return;
    }

    console.log('Sending stop command to server');
    setIsComputing(false);
    wsRef.current.send(JSON.stringify({ action: 'stop' }));
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
    }
    setProblemType(problem_type);
  };

  const createNewRetrosynthesisExperiment = async (startingSmiles: string): Promise<void> => {
    // Close context menu
    setContextMenu({node: null, x: 0, y: 0});

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
          x: e.clientX,
          y: e.clientY
      });
  };

  const sendMessageToServer = (message: string, data?: Omit<WebSocketMessageToServer, 'action'>): void => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    const msg: WebSocketMessageToServer = {
      action: message,
      ...data
    };
    wsRef.current.send(JSON.stringify(msg));
    setContextMenu({node: null, x: 0, y: 0});
  };

  const handleCustomQuery = (node: TreeNode, queryType: string | null): void => {
    setCustomQueryModal(node);
    setCustomQueryText('');
    setCustomQueryType(queryType);
    setContextMenu({node: null, x: 0, y: 0});
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
        xpos: customQueryModal?.x
      };
    }

    const message: WebSocketMessageToServer = {
      action: customQueryType ?? (problemType === "optimization" ? "optimize-from" : "recompute-reaction"),
      nodeId: customQueryModal?.id,
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
    <div className="min-h-screen bg-gradient-to-br from-slate-900 via-purple-900 to-slate-900">
      <div className="flex min-h-screen">
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
        <div className="flex-1 min-w-0 p-8">
          <div className="w-full">
            <div className="absolute top-10 text-white">
              <svg version="1.1" id="Layer_1" height="60px" viewBox="0 0 40 40">
                <g>
                  <rect x="1.73" y="0.01" fill="#FFFFFF" width="34.19" height="34.19"/>
                  <path fill="#1E59AE" d="M35.92,0.01v17.53H18.95V0.01H35.92z M15.88,21.82c-1.12-0.07-1.72-0.78-1.79-2.1V0.01h-0.76v19.73
          c0.09,1.72,1,2.75,2.53,2.84h15.28l-4.83,5l-11.79,0c-3.04-0.36-6.22-2.93-6.14-6.98V0.01H7.6V20.6c-0.09,4.49,3.45,7.34,6.86,7.75
          h11.09l-4.59,4.75h-6.68C9.71,32.93,3.19,29.44,2.99,21.13V0.01H0.05v37.3h35.87V17.62l-4.05,4.19L15.88,21.82z"/>
                </g>
              </svg>
            </div>
            <div className="text-center mb-8">
              <div className="flex items-center justify-center gap-3 mb-2">
                <FlaskConical className="w-10 h-10 text-purple-400" />
                <h1 className="text-4xl font-bold text-white">FLASK Copilot</h1>
              </div>
              <p className="text-purple-300">Real-time molecular assistant</p>
            </div>

            <div className="flex justify-end gap-2 mb-4">
              <button onClick={loadContextFromFile} disabled={isComputing} className="px-4 py-2 bg-blue-500/30 text-white rounded-lg text-sm font-semibold hover:bg-blue-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
                </svg>
                Load
              </button>
              <div className="relative">
                <button onClick={(e) => { e.stopPropagation(); setSaveDropdownOpen(!saveDropdownOpen); }} onMouseDown={(e) => e.stopPropagation()} disabled={isComputing || treeNodes.length === 0 || !wsConnected} className="px-4 py-2 bg-blue-500/30 text-white rounded-lg text-sm font-semibold hover:bg-blue-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 7H5a2 2 0 00-2 2v9a2 2 0 002 2h14a2 2 0 002-2V9a2 2 0 00-2-2h-3m-1 4l-3 3m0 0l-3-3m3 3V4" />
                  </svg>
                  Save
                  <svg className={`w-3 h-3 transition-transform ${saveDropdownOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                  </svg>
                </button>
                {saveDropdownOpen && (
                  <div className="absolute top-full mt-2 left-0 bg-slate-800 border-2 border-purple-400 rounded-lg shadow-2xl py-2 min-w-48 z-[100]" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                    <button onClick={saveTree} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors">Save Tree Only</button>
                    <button onClick={requestSaveContext} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors border-t border-purple-400/30">Save Full Context</button>
                  </div>
                )}
              </div>
              <button onClick={() => setSidebarOpen(!sidebarOpen)} className="px-4 py-2 bg-purple-500/30 text-white rounded-lg text-sm font-semibold hover:bg-purple-500/50 transition-all flex items-center gap-2">
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z" />
                </svg>
                Reasoning
                {/* Badges for number of messages
                sidebarMessages.length > 0 && (
                  <span className="bg-pink-500 text-white text-xs rounded-full w-5 h-5 flex items-center justify-center font-bold">{sidebarMessages.length}</span>
                ) */}
              </button>
              <ProfileButton
                initialSettings={profileSettings}
                onSettingsChange={handleProfileUpdateConfirm}
                username={username}
              />
              {/* WebSocket Status Indicator */}
                <div
                  className="absolute top-10 group"
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
                    <div className={`w-4 h-4 rounded-full absolute ${
                      wsReconnecting ? 'bg-yellow-400 animate-ping' :
                      wsConnected ? 'bg-green-400' :
                      'bg-red-400 animate-pulse'
                    }`} />
                    <div className={`w-4 h-4 rounded-full ${
                      wsReconnecting ? 'bg-yellow-400' :
                      wsConnected ? 'bg-green-400' :
                      'bg-red-400 animate-ping'
                    } ${wsTooltipPinned ? 'ring-2 ring-white ring-offset-2 ring-offset-slate-900' : ''}`} />
                  </div>

                  <div className={`absolute right-0 top-8 transition-opacity z-50 ${
                    wsTooltipPinned ? 'opacity-100' : 'opacity-0 group-hover:opacity-100 pointer-events-none'
                  }`}>
                    <div
                      className="bg-slate-800 border-2 border-purple-400 rounded-lg px-3 py-2 text-sm shadow-xl"
                      style={{ minWidth: '220px', maxWidth: '300px' }}
                      onClick={(e) => e.stopPropagation()}
                      onMouseDown={(e) => e.stopPropagation()}
                    >
                      <div className="flex items-center justify-between mb-1">
                        <div className={`font-semibold ${
                          wsReconnecting ? 'text-yellow-400' :
                          wsConnected ? 'text-green-400' :
                          'text-red-400'
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
                            className="text-purple-400 hover:text-white transition-colors cursor-pointer"
                          >
                            <X className="w-3 h-3" />
                          </button>
                        )}
                      </div>
                      <div className="text-purple-200 text-xs mt-1">
                        {WS_SERVER}
                        {wsError && (
                          <div className="text-purple-300 text-xs mt-1">
                            {wsError}
                          </div>
                        )}
                        {wsConnected && availableTools.length > 0 && (
                          <div className="mt-3 pt-2 border-t border-purple-400/30">
                            <div className="text-purple-300 text-xs font-semibold mb-1.5">
                              Available Tool Servers ({availableTools.length})
                            </div>
                            <div className="space-y-1 max-h-60 overflow-y-auto pr-1">
                              {availableTools.map((tool, idx) => (
                                <div key={idx} className="text-xs bg-purple-900/30 rounded px-2 py-1">
                                  <div className="text-purple-100 font-medium">{tool.server || "server" as string}</div>
                                  {tool.names && (
                                    <div className="text-purple-300 mt-0.5 text-[10px] leading-tight">
                                      {tool.names.join(", ")}
                                    </div>
                                  )}
                                </div>
                              ))}
                            </div>
                          </div>
                        )}
                        {wsConnected && availableTools.length === 0 && (
                          <div className="mt-2 text-purple-300 text-xs italic">
                            No MCP tool servers detected!
                          </div>
                        )}
                      </div>
                      {!wsConnected && !wsReconnecting && (
                        <div className="mt-2">
                          <div className="text-purple-300 text-xs italic">
                            Backend server required for computation
                          </div>
                          {wsTooltipPinned && (
                            <button
                              onClick={(e) => {
                                e.stopPropagation();
                                reconnectWS();
                              }}
                              className="mt-2 w-full px-3 py-1.5 bg-purple-600 hover:bg-purple-500 text-white text-xs font-semibold rounded transition-colors flex items-center justify-center gap-1"
                            >
                              <RefreshCw className="w-3 h-3" />
                              Reconnect
                            </button>
                          )}
                        </div>
                      )}
                      {!wsTooltipPinned && (
                        <div className="text-purple-400 text-[10px] mt-2 italic text-center border-t border-purple-400/30 pt-1.5">
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

            <div className="bg-white/10 backdrop-blur-lg rounded-2xl shadow-2xl p-6 mb-6 border border-white/20">
              <div className="flex items-end gap-4 mb-4">
                <div className="flex-1">
                  <label className="block text-sm font-medium text-purple-200 mb-2">Starting Molecule (SMILES)</label>
                  <input type="text" value={smiles} onChange={(e) => setSmiles(e.target.value)} disabled={isComputing} placeholder="Enter SMILES notation" className="w-full px-4 py-3 bg-white/20 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none transition-colors font-mono text-lg text-white placeholder-purple-300/50 disabled:opacity-50" />
                </div>
                <div>
                  <label className="flex items-center gap-2 text-sm text-purple-200 mb-3 cursor-pointer">
                    <input type="checkbox" checked={autoZoom} onChange={(e) => setAutoZoom(e.target.checked)} className="w-4 h-4 rounded border-purple-400/50 bg-white/20 text-purple-600 focus:ring-purple-500 focus:ring-offset-0 disabled:opacity-50" />
                    Auto-zoom to fit
                  </label>
                </div>
              </div>

              <div className="flex items-center gap-4">
                <div className="flex items-end gap-2">
                  <div>
                    <label className="block text-sm font-medium text-purple-200 mb-2">
                      Problem Type
                        {problemType === "custom" && (!systemPrompt || !problemPrompt) && (
                          <span className="ml-2 text-amber-400 cursor-help relative group inline-block">
                            ⚠️
                          <div className="absolute bottom-full mb-2 left-1/2 transform -translate-x-1/2 opacity-0 group-hover:opacity-90 transition-opacity pointer-events-none">
                            <div className="bg-slate-800 border-2 border-purple-400 rounded-lg px-4 py-2 text-sm whitespace-nowrap shadow-xl text-purple-200">
                              Custom problem description not given
                            </div>
                          </div>
                          </span>
                        )}
                    </label>
                    <select value={problemType} onChange={(e) => {reset(); resetProblemType(e.target.value)}} disabled={isComputing} className="w-full px-4 py-2.5 bg-white/20 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none transition-colors text-white disabled:opacity-50 cursor-pointer text-sm">
                      <option value="retrosynthesis" className="bg-slate-800">Retrosynthesis</option>
                      <option value="optimization" className="bg-slate-800">Lead Molecule Optimization</option>
                      <option value="custom" className="bg-slate-800">Custom</option>
                    </select>
                  </div>
                  { /* Problem-specific UI */ }
                  {problemType === "optimization" &&
                    <div>
                      <label className="block text-sm font-medium text-purple-200 mb-2">
                        Property
                        {propertyType === "custom" && (!customPropertyName || !customPropertyDesc) && (
                          <span className="ml-2 text-amber-400 cursor-help relative group inline-block">
                            ⚠️
                          <div className="absolute bottom-full mb-2 left-1/2 transform -translate-x-1/2 opacity-0 group-hover:opacity-90 transition-opacity pointer-events-none">
                            <div className="bg-slate-800 border-2 border-purple-400 rounded-lg px-4 py-2 text-sm whitespace-nowrap shadow-xl text-purple-200">
                              Property name or description not given
                            </div>
                          </div>
                          </span>
                        )}
                      </label>
                      <select value={propertyType} onChange={(e) => {setPropertyType(e.target.value)}} disabled={isComputing} className="w-48 px-4 py-2.5 bg-white/20 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none transition-colors text-white disabled:opacity-50 cursor-pointer text-sm">
                        <option value="density" className="bg-slate-800">Crystalline Density</option>
                        <option value="hof" className="bg-slate-800">Heat of Formation</option>
                        <option value="bandgap" className="bg-slate-800">Band Gap</option>
                        <option value="custom" className="bg-slate-800">Other</option>
                      </select>
                    </div>
                  }
                  {problemType === "optimization" && propertyType === "custom" &&
                    <button onClick={() => setEditPropertyModal(true)} disabled={isComputing} className="px-3 py-2.5 bg-white/10 text-purple-200 rounded-lg text-sm font-medium hover:bg-white/20 transition-all disabled:opacity-50 disabled:cursor-not-allowed">
                      Property...
                    </button>
                  }
                  {problemType === "custom" &&
                    <button onClick={() => setEditPromptsModal(true)} disabled={isComputing || problemType !== "custom"} className="px-3 py-2.5 bg-white/10 text-purple-200 rounded-lg text-sm font-medium hover:bg-white/20 transition-all disabled:opacity-50 disabled:cursor-not-allowed">
                      Edit
                    </button>
                  }
                  <button
                    onClick={() => setShowToolSelectionModal(true)}
                    disabled={isComputing}
                    className="px-5 py-2.5 bg-white/10 text-purple-200 rounded-lg text-sm font-medium hover:bg-white/20 transition-all disabled:opacity-50 disabled:cursor-not-allowed"
                  >
                    Select Tools {selectedTools.length > 0 && `(${selectedTools.length})`}
                  </button>
                </div>

                <div className="flex gap-3 flex-1 flex-wrap justify-end">
                  <div className="relative group">
                    <button onClick={runComputation} disabled={!wsConnected || isComputing || !smiles} className="px-4 py-2 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all disabled:opacity-50 disabled:cursor-not-allowed disabled:transform-none flex items-center gap-2">
                      {isComputing ? <><Loader2 className="w-5 h-5 animate-spin" />Computing</> : <><Play className="w-5 h-5" />Run</>}
                    </button>
                    {(!wsConnected || isComputing || !smiles) && (
                      <div className="absolute bottom-full mb-2 left-1/2 transform -translate-x-1/2 opacity-0 group-hover:opacity-100 transition-opacity pointer-events-none">
                        <div className="bg-slate-800 border-2 border-purple-400 rounded-lg px-4 py-2 text-sm whitespace-nowrap shadow-xl">
                          <div className="text-purple-200">
                            {!wsConnected ? 'Backend server not connected' :
                            isComputing ? 'Computation already running' :
                            'Enter a SMILES string first'}
                          </div>
                          {!wsConnected && (
                            <div className="text-purple-300 text-xs mt-1">
                              Start the backend server at {WS_SERVER}
                            </div>
                          )}
                        </div>
                      </div>
                    )}
                  </div>
                  <button onClick={() => {
                    const [updatedNodes, updatedEdges] = relayoutTree(treeNodes, edges);
                    setTreeNodes(updatedNodes);
                    setEdges(updatedEdges);
                  }} disabled={isComputing || treeNodes.length === 0} className="px-4 py-2 bg-purple-500/30 text-white rounded-lg font-semibold hover:bg-purple-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                    <Sparkles className="w-5 h-5" />Relayout
                  </button>
                  <button onClick={reset} disabled={isComputing} className="px-4 py-2 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                    <RotateCcw className="w-5 h-5" />Reset
                  </button>
                  <button onClick={stop} disabled={!isComputing} className="px-4 py-2 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                    <X className="w-5 h-5" />Stop
                  </button>
                </div>
              </div>
            </div>

            <div className="bg-white/10 backdrop-blur-lg rounded-2xl shadow-2xl border border-white/20 overflow-hidden relative" style={{ height: '600px' }}>
              {treeNodes.length === 0 && !isComputing ? (
                <div className="flex flex-col items-center justify-center h-full text-purple-300">
                  <FlaskConical className="w-16 h-16 mb-4 opacity-50" />
                  <p className="text-center text-lg">
                    {wsConnected ?
                      `Click "Run" to start ${problemType === "optimization" ? "molecular discovery" : "the molecular computation tree"}` :
                      "Waiting for backend connection..."
                    }
                  </p>
                  <p className="text-sm text-purple-400 mt-2">
                    {autoZoom ? 'Auto-zoom will fit all molecules' : 'Drag to pan • Scroll to zoom'}
                  </p>
                  {!wsConnected && (
                    <div className="mt-4 bg-amber-500/20 border border-amber-400/50 rounded-lg p-4 max-w-md">
                      <div className="text-amber-200 text-sm text-center">
                        <strong>Backend Required:</strong> Start your Python backend server at <code className="bg-black/30 px-2 py-1 rounded">{WS_SERVER}</code> to enable molecular computations.
                      </div>
                    </div>
                  )}
                </div>
              ) : (
                <MoleculeGraph {...graphState} nodes={treeNodes} edges={edges} autoZoom={autoZoom} setAutoZoom={setAutoZoom} ctx={contextMenu} handleNodeClick={handleNodeClick} rdkitModule={rdkitModule} />
              )}
            </div>

            {isComputing && (
              <div className="mt-6 bg-purple-500/20 backdrop-blur-lg rounded-xl p-4 border border-purple-400/50 animate-pulse">
                <div className="flex items-center gap-3 text-purple-200">
                  <Loader2 className="w-5 h-5 animate-spin" />
                  <span className="font-medium">Streaming molecules... {treeNodes.length} nodes discovered</span>
                </div>
              </div>
            )}

            {!isComputing && treeNodes.length > 0 && (
              <div className="mt-6 bg-green-500/20 backdrop-blur-lg rounded-xl p-4 border border-green-400/50">
                <div className="flex items-center gap-3 text-green-200">
                  <span className="font-medium">
                    Computation complete! Generated {treeNodes.length} molecules
                  </span>
                </div>
              </div>
            )}

            {/* Metrics Dashboard */}
            {treeNodes.length > 0 && (
              <MetricsDashboard {...metricsDashboardState} treeNodes={treeNodes} />
            )}


          <div className="mt-8 pt-6 border-t border-purple-400/30 text-center text-purple-300 text-sm">
            <p>This work was performed under the auspices of the U.S. Department of Energy
            by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344
            (LLNL-CODE-2006345).</p>
            {VERSION && (
              <p>Server version: {VERSION}</p>
            )}
          </div>
          </div>
        </div>

        {sidebarOpen && (
          <ReasoningSidebar {...sidebarState} setSidebarOpen={setSidebarOpen} rdkitModule={rdkitModule} />
        )}
      </div>

      {contextMenu && contextMenu.node && (
        <div className="fixed z-50 bg-slate-800 border-2 border-purple-400 rounded-lg shadow-2xl py-2 min-w-48" style={{ left: `${contextMenu.x + 10}px`, top: `${contextMenu.y + 10}px` }} onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
          <div className="px-3 py-2 border-b border-purple-400/30">
            <div className="text-xs text-purple-300">Actions for</div>
            <div className="text-sm font-semibold text-white">{contextMenu.node.label}</div>
          </div>

          { (problemType === "optimization") && (
            <>
            <button onClick={() => {
              const nodeId = contextMenu.node!.id;
                // Delete all nodes from this point on
                setTreeNodes(prev => {
                  return prev.filter(n => n.x <= contextMenu.node!.x);
                });
                sendMessageToServer("optimize-from", {nodeId: nodeId, propertyType, customPropertyName, customPropertyDesc, customPropertyAscending, smiles: contextMenu.node!.smiles, xpos: contextMenu.node!.x});
                setIsComputing(true);
              }}  className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
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
              }}  className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <MessageSquareShare className="w-4 h-4" />
              Refine search (with prompt)
            </button>
            <button onClick={() => { createNewRetrosynthesisExperiment(contextMenu.node!.smiles); }}  className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <FlaskConical className="w-4 h-4" />
              Plan synthesis pathway
            </button>
            </>
          )}

          <button
            onClick={() => copyToClipboard(contextMenu.node!.smiles, 'smiles', setCopiedField)}
            className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2"
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

          { (problemType === "retrosynthesis" && !hasDescendants(contextMenu.node.id, treeNodes)) && (
            <button onClick={() => {
              sendMessageToServer("compute-reaction-from", {nodeId: contextMenu.node!.id});
            }} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <TestTubeDiagonal className="w-4 h-4" />
              How do I make this?
            </button>
            ) }

          { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) && (
            <>
            <button disabled={true} onClick={() => {sendMessageToServer("recompute-reaction", {nodeId: contextMenu.node!.id});}} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 disabled:text-gray-400 disabled:cursor-not-allowed disabled:hover:bg-transparent">
              <RefreshCw className="w-4 h-4" />Find Another Reaction
            </button>
            </>
          )}
          { (problemType === "retrosynthesis" && !isRootNode(contextMenu.node.id, treeNodes)) && (
            <>
            <button disabled={true} onClick={() => {sendMessageToServer("recompute-parent-reaction", {nodeId: contextMenu.node!.id});}} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 disabled:text-gray-400 disabled:cursor-not-allowed disabled:hover:bg-transparent">
              <Network className="w-4 h-4" />Substitute Molecule
            </button>
            </>
          )}

          { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) && (
          <button disabled={true} onClick={() => handleCustomQuery(contextMenu.node!, null)} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 border-t border-purple-400/30 disabled:text-gray-400 disabled:cursor-not-allowed disabled:hover:bg-transparent">
            <Send className="w-4 h-4" />
            { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) ? (<>Find Another Reaction with Custom Prompt...</>) : (<>Custom Query...</>) }
          </button>
          )}

          { (problemType === "retrosynthesis" && contextMenu.node) && (
          <button onClick={() => handleCustomQuery(contextMenu.node!, "query-retro-molecule")} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 border-t border-purple-400/30">
            <MessageCircleQuestion className="w-4 h-4" /> Ask about molecule...
          </button>
          )}
          { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) && (
          <button onClick={() => handleCustomQuery(contextMenu.node!, "query-retro-product")} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 border-t border-purple-400/30">
            <MessageCircleQuestion className="w-4 h-4" /> Ask about reaction product...
          </button>
          )}
          { (problemType === "retrosynthesis" && !isRootNode(contextMenu.node.id, treeNodes)) && (
          <button onClick={() => handleCustomQuery(contextMenu.node!, "query-retro-reactant")} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 border-t border-purple-400/30">
            <MessageCircleQuestion className="w-4 h-4" /> Ask about reactant...
          </button>
          )}
        </div>
      )}

      {customQueryModal && (
        <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-2xl w-full p-6">
            <div className="flex items-center justify-between mb-4">
              <div>
                <h2 className="text-xl font-bold text-white">Custom Query</h2>
                <p className="text-sm text-purple-300">for {customQueryModal.label}</p>
              </div>
              <button onClick={() => setCustomQueryModal(null)} className="text-purple-300 hover:text-white transition-colors">
                <X className="w-6 h-6" />
              </button>
            </div>

            <textarea
              value={customQueryText}
              onChange={(e) => setCustomQueryText(e.target.value)}
              placeholder="Enter your custom query here..."
              className="w-full h-40 px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 resize-none"
            />

            <div className="flex gap-3 mt-4">
              <button onClick={submitCustomQuery} disabled={!customQueryText.trim()} className="flex-1 px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all disabled:opacity-50 disabled:cursor-not-allowed disabled:transform-none flex items-center justify-center gap-2">
                <Send className="w-5 h-5" />
                Submit Query
              </button>

              <button onClick={() => setCustomQueryModal(null)} className="px-6 py-3 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all">
                Cancel
              </button>
            </div>
          </div>
        </div>
      )}

      {editPromptsModal && (
        <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-3xl w-full p-6">
            <div className="flex items-center justify-between mb-4">
              <div>
                <h2 className="text-xl font-bold text-white">Edit Prompts</h2>
                <p className="text-sm text-purple-300">Configure system and problem-specific prompts</p>
              </div>
              <button onClick={() => setEditPromptsModal(false)} className="text-purple-300 hover:text-white transition-colors">
                <X className="w-6 h-6" />
              </button>
            </div>

            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">System Prompt</label>
                <textarea value={systemPrompt} onChange={(e) => setSystemPrompt(e.target.value)} placeholder="Enter system-level instructions..." className="w-full h-32 px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 resize-none" />
              </div>

              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">Problem Prompt</label>
                <textarea value={problemPrompt} onChange={(e) => setProblemPrompt(e.target.value)} placeholder="Enter problem-specific instructions..." className="w-full h-32 px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 resize-none" />
              </div>
            </div>

            <div className="flex gap-3 mt-4">
              <button onClick={() => { savePrompts('', ''); }} className="px-4 py-3 bg-white/10 text-purple-200 rounded-lg font-medium hover:bg-white/20 transition-all flex items-center gap-2">
                <RotateCcw className="w-4 h-4" />
                Reset
              </button>
              <button onClick={() => {savePrompts(systemPrompt, problemPrompt); setProblemType('custom');}} className="flex-1 px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all flex items-center justify-center gap-2">
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
        <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-3xl w-full p-6">
            <div className="flex items-center justify-between mb-4">
              <div>
                <h2 className="text-xl font-bold text-white">Custom Property</h2>
                <p className="text-sm text-purple-300">Configure custom property type and units</p>
              </div>
              <button onClick={() => setEditPropertyModal(false)} className="text-purple-300 hover:text-white transition-colors">
                <X className="w-6 h-6" />
              </button>
            </div>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">Property Name</label>
                <input
                    type="text"
                    value={customPropertyName}
                    onChange={(e) => setCustomPropertyName(e.target.value)}
                    placeholder="Enter a name for the property"
                    className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50"
                  />
              </div>
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">Property Description</label>
                <textarea value={customPropertyDesc} onChange={(e) => setCustomPropertyDesc(e.target.value)} placeholder="Enter a description of the property and its units..." className="w-full h-32 px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 resize-none" />
              </div>
            </div>
            <div className="py-4 flex items-center justify-center gap-3">
              <span className="text-sm text-purple-200">Higher is better</span>
              <button
                onClick={() => setCustomPropertyAscending(!customPropertyAscending)}
                className={`relative w-14 h-7 rounded-full transition-colors ${customPropertyAscending ? 'bg-purple-400/30' : 'bg-purple-600'}`}
              >
                <div
                  className={`absolute top-1 left-1 w-5 h-5 bg-white rounded-full transition-transform ${customPropertyAscending ? 'translate-x-0' : 'translate-x-7'}`}
                />
              </button>
              <span className="text-sm text-purple-200">Lower is better</span>
            </div>
            <div className="flex gap-3 mt-4">
              <button onClick={() => { saveCustomProperty('', '', true); }} className="px-4 py-3 bg-white/10 text-purple-200 rounded-lg font-medium hover:bg-white/20 transition-all flex items-center gap-2">
                <RotateCcw className="w-4 h-4" />
                Reset
              </button>
              <button onClick={() => {saveCustomProperty(customPropertyName, customPropertyDesc, customPropertyAscending);}} className="flex-1 px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all flex items-center justify-center gap-2">
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
