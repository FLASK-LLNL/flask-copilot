/* eslint-disable react-hooks/exhaustive-deps */
import React, { useState, useRef, useEffect, useMemo, useCallback } from 'react';
import { Loader2, FlaskConical, TestTubeDiagonal, Network, Play, RotateCcw, Move, X, Send, RefreshCw, Sparkles } from 'lucide-react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

import { WS_SERVER } from './config';
import { TreeNode, Edge, ContextMenuState, MetricHistoryItem, VisibleMetrics, SidebarMessage, VisibleSources, Tool, WebSocketMessageToServer, WebSocketMessage, MetricDefinitions } from './types';

import { MoleculeSVG, loadRDKit } from './components/molecule';
import { MarkdownText } from './components/markdown';
import { MoleculeGraph, useGraphState } from './components/graph';

import { findAllDescendants, hasDescendants, isRootNode, relayoutTree } from './tree_utils';

import './animations.css';

const ChemistryTool: React.FC = () => {
  const [smiles, setSmiles] = useState<string>('O=C\\C1=C(\\C=C/CC1(C)C)C');
  const [problemType, setProblemType] = useState<string>('retrosynthesis');
  const [problemName, setProblemName] = useState<string>('retro-safranal');
  const [systemPrompt, setSystemPrompt] = useState<string>('');
  const [problemPrompt, setProblemPrompt] = useState<string>('');
  const [promptsModified, setPromptsModified] = useState<boolean>(false);
  const [editPromptsModal, setEditPromptsModal] = useState<boolean>(false);
  const [isComputing, setIsComputing] = useState<boolean>(false);
  const [autoZoom, setAutoZoom] = useState<boolean>(true);
  const [treeNodes, setTreeNodes] = useState<TreeNode[]>([]);
  const [edges, setEdges] = useState<Edge[]>([]);
  const [contextMenu, setContextMenu] = useState<ContextMenuState>({node: null, x: 0, y: 0});
  const [customQueryModal, setCustomQueryModal] = useState<TreeNode | null>(null);
  const [customQueryText, setCustomQueryText] = useState<string>('');
  const [metricsHistory, setMetricsHistory] = useState<MetricHistoryItem[]>([]);
  const [websocket, setWebsocket] = useState<WebSocket | null>(null);
  const [wsConnected, setWsConnected] = useState<boolean>(false);
  const [saveDropdownOpen, setSaveDropdownOpen] = useState<boolean>(false);
  const [wsError, setWsError] = useState<string>('');
  const [wsReconnecting, setWsReconnecting] = useState<boolean>(false);
  const [visibleMetrics, setVisibleMetrics] = useState<VisibleMetrics>({
    cost: true,
    bandgap: false,
    density: false,
    yield: false,
  });
  const rdkitModule = loadRDKit();
  const [sidebarOpen, setSidebarOpen] = useState<boolean>(false);
  const [sidebarMessages, setSidebarMessages] = useState<SidebarMessage[]>([]);
  const [copiedField, setCopiedField] = useState<string | null>(null);
  const [visibleSources, setVisibleSources] = useState<VisibleSources>({
    'System': true,
    'Reasoning': true,
    'Logger (Error)': true,
    'Logger (Warning)': true,
    'Logger (Info)': false,
    'Logger (Debug)': false
  });
  const [sourceFilterOpen, setSourceFilterOpen] = useState<boolean>(false);
  const [availableTools, setAvailableTools] = useState<Tool[]>([]);
  const [wsTooltipPinned, setWsTooltipPinned] = useState<boolean>(false);

  const wsRef = useRef<WebSocket | null>(null);
  const sidebarRef = useRef<HTMLDivElement>(null);

  const graphState = useGraphState();
  
  const copyToClipboard = async (text: string, fieldName: string): Promise<void> => {
    try {
      await navigator.clipboard.writeText(text);
      setCopiedField(fieldName);
      setTimeout(() => setCopiedField(null), 2000);
    } catch (err) {
      console.error('Failed to copy:', err);
    }
  };

  useEffect(() => {
    const handleClickOutside = (): void => {
      setContextMenu({node: null, x:0, y:0});
      setSaveDropdownOpen(false);
      setSourceFilterOpen(false);
      setWsTooltipPinned(false);
      setCopiedField(null);
    };
    if (contextMenu || saveDropdownOpen || sourceFilterOpen || wsTooltipPinned) {
      window.addEventListener('mousedown', handleClickOutside);
      return () => window.removeEventListener('mousedown', handleClickOutside);
    }
  }, [contextMenu, saveDropdownOpen, sourceFilterOpen, wsTooltipPinned]);
  
  const runComputation = async (): Promise<void> => {
    setSidebarOpen(true);
    setIsComputing(true);
    setTreeNodes([]);
    setEdges([]);
    graphState.setOffset({ x: 50, y: 50 });
    graphState.setZoom(1);

    const message: WebSocketMessageToServer = {
      action: 'compute',
      smiles: smiles,
      problemType: problemType
    };

    websocket?.send(JSON.stringify(message));
  };

  const reconnectWS = (): void => {
    if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
      wsRef.current.close();
    }
    
    setWsReconnecting(true);
    
    const socket = new WebSocket(WS_SERVER);
    wsRef.current = socket;
    
    socket.onopen = () => {
      console.log('WebSocket connected');
      setWebsocket(socket);
      setWsConnected(true);
      setWsReconnecting(false);
      setWsError('');
      reset();  // Server state must match UI state

      socket.send(JSON.stringify({ action: 'list-tools' }));
    };
    
    socket.onmessage = (event: MessageEvent) => {
      const data: WebSocketMessage = JSON.parse(event.data);
      
      if (data.type === 'node') {
        setTreeNodes(prev => [...prev, data.node!]);
      } else if (data.type === 'node_update') {
        const { id, ...restData } = data.node!;

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
      } else if (data.type === 'response') {
        addSidebarMessage(data.message!);
        console.log('Server response:', data.message);
      } else if (data.type === 'available-tools-response') {
        setAvailableTools(data.tools || []);
      } else if (data.type === 'error') {
        console.error(data.message);
        alert("Server error: " + data.message);
      }
    };
    
    socket.onerror = (error: Event) => {
      console.error('WebSocket error:', error);
      setWsReconnecting(false);
      setIsComputing(false);
      setWsError((error as any).message || 'Connection failed');
      setAvailableTools([]);
    };

    socket.onclose = () => {
      console.log('WebSocket closed');
      wsRef.current = null;
      setWebsocket(null);
      setWsConnected(false);
      setIsComputing(false);
      setWsReconnecting(false);
      setAvailableTools([]);
    };
  };

  // Connect WebSocket on mount
  useEffect(() => {
    reconnectWS();
  }, []);

  const reset = (): void => {
    setTreeNodes([]);
    setEdges([]);
    setIsComputing(false);
    graphState.setOffset({ x: 50, y: 50 });
    graphState.setZoom(1);
    setContextMenu({node: null, x:0, y:0});
    setCustomQueryModal(null);
    setMetricsHistory([]);
    setSidebarMessages([]);
    setSaveDropdownOpen(false);
    setSourceFilterOpen(false);
    setWsTooltipPinned(false);
    if (websocket && websocket.readyState === WebSocket.OPEN) {
      websocket.send(JSON.stringify({ action: 'reset' }));
    }
  };

  const stop = (): void => {
    if (!websocket || websocket.readyState !== WebSocket.OPEN) {
      // alert('WebSocket not connected');
      return;
    }

    console.log('Sending stop command to server');
    setIsComputing(false);
    websocket.send(JSON.stringify({ action: 'stop' }));
  };

  // TODO: Improve saveTree and saveFullContext by handing the information to/from the backend
  const saveTree = (): void => {
    const data = { version: '1.0', type: 'tree', timestamp: new Date().toISOString(), smiles, problemType, systemPrompt, problemPrompt, nodes: treeNodes, edges };
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `molecule-tree-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
    setSaveDropdownOpen(false);
  };

  const saveFullContext = (): void => {
    const data = { version: '1.0', type: 'full-context', timestamp: new Date().toISOString(), smiles, problemType, systemPrompt, problemPrompt, nodes: treeNodes, edges, metricsHistory, visibleMetrics, sidebarMessages, visibleSources };
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `molecule-context-${Date.now()}.json`;
    a.click();
    URL.revokeObjectURL(url);
    setSaveDropdownOpen(false);
  };

  const loadContext = (): void => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.json';
    input.onchange = (e: Event) => {
      const file = (e.target as HTMLInputElement).files?.[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = (e: ProgressEvent<FileReader>) => {
        try {
          const data = JSON.parse(e.target?.result as string);
          if (data.type === 'tree' || data.type === 'full-context') {
            setSmiles(data.smiles || 'CCO');
            setProblemType(data.problemType || 'retrosynthesis');
            setSystemPrompt(data.systemPrompt || '');
            setProblemPrompt(data.problemPrompt || '');
            setPromptsModified(!!(data.systemPrompt || data.problemPrompt));
            setTreeNodes(data.nodes || []);
            setEdges(data.edges || []);
            if (data.type === 'full-context') {
              setMetricsHistory(data.metricsHistory || []);
              setVisibleMetrics(data.visibleMetrics || { cost: true, bandgap: false, yield: false, density: false });
              graphState.setZoom(data.zoom || 1);
              graphState.setOffset(data.offset || { x: 50, y: 50 });
              setSidebarMessages(data.sidebarMessages || []);
              setVisibleSources(data.visibleSources || {
                'System': true,
                'Reasoning': true,
                'Logger (Error)': true,
                'Logger (Warning)': true,
                'Logger (Info)': false,
                'Logger (Debug)': false });
            }
          } else {
            alert('Invalid file format');
          }
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
    setPromptsModified(!!(newSystemPrompt || newProblemPrompt));
    setEditPromptsModal(false);
  };

  const resetProblemType = (problem_name: string): void => {
    setSystemPrompt('');
    setProblemPrompt('');
    setPromptsModified(false);
    setProblemName(problem_name);
    if (problem_name === 'retro-safranal') {
      setSmiles('O=C\\C1=C(\\C=C/CC1(C)C)C')
      setVisibleMetrics({cost: false, bandgap: false, yield: false, density: false});
      setProblemType('retrosynthesis');
    } else if (problem_name === 'retro-nirmatrelvir') {
      setSmiles('CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C')
      setVisibleMetrics({cost: false, bandgap: false, yield: false, density: false});
      setProblemType('retrosynthesis');
    } else if (problem_name === 'optimization-bandgap') {
      setSmiles('C1(N2C3=CC=CC=C3N(C4=CC=CC=C4)C5=C2C=CC=C5)=CC=CC=C1');
      setVisibleMetrics({cost: false, bandgap: true, yield: false, density: false});
      setProblemType('optimization');
    } else if (problem_name === 'optimization-density') {
      setSmiles('CCC(=O)CCc1cccc(CCC(=O)CC)c1ONC(C)=O');
      setVisibleMetrics({cost: false, bandgap: false, yield: false, density: true});
      setProblemType('optimization');
    }
  };


  // Metric definitions for extensibility
  const metricDefinitions: MetricDefinitions = {
    cost: { 
      label: 'Reaction Cost ($)', 
      color: '#EC4899',
      calculate: (nodes: TreeNode[]) => nodes.reduce((sum, node) => sum + (node.cost || 0), 0)
    },
    bandgap: { 
      label: 'Band Gap (eV)', 
      color: '#F59E0B',
      calculate: (nodes: TreeNode[]) => nodes[nodes.length-1]?.bandgap || 0
    },
    density: { 
      label: 'Molecular Density (g/cm³)', 
      color: '#10B981',
      calculate: (nodes: TreeNode[]) => nodes[nodes.length-1]?.density || 0
    },
    /*yield: { 
      label: 'Yield (%)', 
      color: '#10B981',
      calculate: (nodes) => nodes.length > 0 ? (nodes.reduce((sum, node) => sum + (node.yield || Math.random() * 100), 0) / nodes.length) : 0
    },*/
  };

  // Calculate all metrics
  const calculateMetrics = (nodes: TreeNode[]): MetricHistoryItem => {
    const metrics: MetricHistoryItem = { nodeCount: nodes.length, step: 0 };
    Object.keys(metricDefinitions).forEach(key => {
      metrics[key] = metricDefinitions[key].calculate(nodes);
    });
    return metrics;
  };

  // Update metrics history whenever tree changes
  useEffect(() => {
    if (treeNodes.length > 0) {
      let metrics = calculateMetrics(treeNodes);
      setMetricsHistory(prev => {
        metrics.step = prev.length;
        return [...prev, { ...metrics }];
      });
    }
  }, [treeNodes.length]);

  // Auto-scroll sidebar to bottom when new messages arrive
  useEffect(() => {
    if (sidebarRef.current && sidebarMessages.length > 0) {
      // Delay scroll to allow molecules to render
      setTimeout(() => {
        if (sidebarRef.current) {
          sidebarRef.current.scrollTo({
            top: sidebarRef.current.scrollHeight,
            behavior: 'smooth'
          });
        }
      }, 100);
    }
  }, [sidebarMessages]);
  
  const handleNodeClick = (e: React.MouseEvent<HTMLDivElement>, node: TreeNode): void => {
      e.stopPropagation();
      if (isComputing) return; // Don't open menu while computing
      setContextMenu({
          node,
          x: e.clientX,
          y: e.clientY
      });
  };

  const sendMessageToServer = (message: string, nodeId: string): void => {
    if (!websocket || websocket.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    const msg: WebSocketMessageToServer = {
      action: message,
      nodeId: nodeId
    };
    websocket.send(JSON.stringify(msg));
    setContextMenu({node: null, x: 0, y: 0});
  };

  const handleCustomQuery = (node: TreeNode): void => {
    if (problemType === "optimization") {
      // Invalidate all nodes below this one
      setTreeNodes(prev => {
        return prev.filter(n => n.y <= node.y);
      });
    }
    setCustomQueryModal(node);
    setCustomQueryText('');
    setContextMenu({node: null, x: 0, y: 0});
  };

  const submitCustomQuery = (): void => {
    if (!websocket || websocket.readyState !== WebSocket.OPEN) {
      alert('WebSocket not connected');
      return;
    }
    console.log(`Custom query for ${customQueryModal?.label}: ${customQueryText}`);
    
    const message: WebSocketMessageToServer = {
      action: problemType === "optimization" ? "optimize-from" : "recompute-reaction",
      nodeId: customQueryModal?.id,
      query: customQueryText
    };
    websocket.send(JSON.stringify(message));

    setCustomQueryModal(null);
    setCustomQueryText('');
    setIsComputing(true); // If expecting new nodes
  };

  const addSidebarMessage = (message: SidebarMessage): void => {
    message.id = Date.now();
    message.timestamp = new Date().toISOString();
    if (!message.source) {
      message.source = "Backend";
    }

    setSidebarMessages(prev => [...prev, message]);
    setSidebarOpen(true);
    
    setVisibleSources(prev => {
      if (!(message.source in prev)) {
        return { ...prev, [message.source]: true };
      }
      return prev;
    });
  };

  // Memoize the metrics charts to prevent re-render on mouse move
  const metricsCharts = useMemo(() => {
    if (metricsHistory.length === 0) return null;
    
    return (
      <div className={`grid gap-6 ${Object.values(visibleMetrics).filter(Boolean).length > 1 ? 'grid-cols-2' : 'grid-cols-1'}`}>
        {Object.keys(metricDefinitions).filter(key => visibleMetrics[key as keyof VisibleMetrics]).map(metricKey => {
          const metric = metricDefinitions[metricKey];
          return (
            <div key={metricKey} className="bg-white/5 rounded-xl p-4">
              <h4 className="text-sm font-semibold text-purple-200 mb-2">{metric.label}</h4>
              <ResponsiveContainer width="100%" height={200}>
                <LineChart data={metricsHistory}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#8B5CF6" opacity={0.1} />
                  <XAxis 
                    dataKey="step" 
                    stroke="#A78BFA"
                    tick={{ fontSize: 12 }}
                  />
                  <YAxis 
                    stroke="#A78BFA"
                    tick={{ fontSize: 12 }}
                  />
                  <Tooltip 
                    contentStyle={{ 
                      backgroundColor: '#1E1B4B', 
                      border: '2px solid #8B5CF6',
                      borderRadius: '8px',
                      color: '#E9D5FF'
                    }}
                    // NOTE: Returning `as any` since no other type information worked
                    formatter={(value: number) => value.toFixed(2) as any}
                  />
                  <Line 
                    type="monotone" 
                    dataKey={metricKey}
                    stroke={metric.color}
                    strokeWidth={3}
                    dot={{ fill: metric.color, r: 3 }}
                    activeDot={{ r: 5 }}
                  />
                </LineChart>
              </ResponsiveContainer>
              <div className="mt-2 text-center">
                <div className="text-2xl font-bold text-white">
                  {metricsHistory[metricsHistory.length - 1][metricKey].toFixed(2)}
                </div>
                <div className="text-xs text-purple-300">Current Value</div>
              </div>
            </div>
          );
        })}
      </div>
    );
  }, [metricsHistory, visibleMetrics]);

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-900 via-purple-900 to-slate-900">
     <div className="flex min-h-screen">
     <div className={`p-8 ${sidebarOpen ? 'flex-1' : 'w-full'}`}>
      <div className="max-w-7xl mx-auto">
        <div className="text-center mb-8">
          <div className="flex items-center justify-center gap-3 mb-2">
            <FlaskConical className="w-10 h-10 text-purple-400" />
            <h1 className="text-4xl font-bold text-white">FLASK Copilot</h1>
          </div>
          <p className="text-purple-300">Real-time molecular assistant</p>
        </div>

        <div className="flex justify-end gap-2 mb-4">
          <button onClick={loadContext} disabled={isComputing} className="px-4 py-2 bg-blue-500/30 text-white rounded-lg text-sm font-semibold hover:bg-blue-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12" />
            </svg>
            Load
          </button>
          <div className="relative">
            <button onClick={(e) => { e.stopPropagation(); setSaveDropdownOpen(!saveDropdownOpen); }} onMouseDown={(e) => e.stopPropagation()} disabled={isComputing || treeNodes.length === 0} className="px-4 py-2 bg-blue-500/30 text-white rounded-lg text-sm font-semibold hover:bg-blue-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
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
                <button onClick={saveFullContext} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors border-t border-purple-400/30">Save Full Context</button>
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
                          Available Tools ({availableTools.length})
                        </div>
                        <div className="space-y-1 max-h-60 overflow-y-auto pr-1">
                          {availableTools.map((tool, idx) => (
                            <div key={idx} className="text-xs bg-purple-900/30 rounded px-2 py-1">
                              <div className="text-purple-100 font-medium">{tool.name || tool as string}</div>
                              {tool.description && (
                                <div className="text-purple-300 mt-0.5 text-[10px] leading-tight">
                                  {tool.description}
                                </div>
                              )}
                            </div>
                          ))}
                        </div>
                      </div>
                    )}
                    {wsConnected && availableTools.length === 0 && (
                      <div className="mt-2 text-purple-300 text-xs italic">
                        Loading tools...
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
              <div className="w-100">
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Problem Type
                  {promptsModified && <span className="ml-2 text-xs text-amber-400">●</span>}
                </label>
                <select value={problemName} onChange={(e) => {reset(); resetProblemType(e.target.value)}} disabled={isComputing} className="w-full px-4 py-2.5 bg-white/20 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none transition-colors text-white disabled:opacity-50 cursor-pointer text-sm">
                  <option value="retro-safranal" className="bg-slate-800">Synthesizing Safranal</option>
                  <option value="retro-nirmatrelvir" className="bg-slate-800">Synthesizing Nirmatrelvir</option>
                  <option value="optimization-bandgap" className="bg-slate-800">Optimizing OLED Molecule for Band Gap</option>
                  <option value="optimization-density" className="bg-slate-800">Molecular Discovery for Crystalline Density</option>
                  {/*<option value="custom" className="bg-slate-800">Custom</option>*/}
                </select>
              </div>
              <button onClick={() => setEditPromptsModal(true)} disabled={true} className="px-3 py-2.5 bg-white/10 text-purple-200 rounded-lg text-sm font-medium hover:bg-white/20 transition-all disabled:opacity-50 disabled:cursor-not-allowed">
                Edit
              </button>
            </div>
            
            <div className="flex gap-3 flex-1 justify-end">
              <div className="relative group">
                <button onClick={runComputation} disabled={!wsConnected || isComputing || !smiles} className="px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all disabled:opacity-50 disabled:cursor-not-allowed disabled:transform-none flex items-center gap-2">
                  {isComputing ? <><Loader2 className="w-5 h-5 animate-spin" />Computing</> : <><Play className="w-5 h-5" />Run</>}
                </button>
                {(!wsConnected || isComputing || !smiles) && (
                  <div className="absolute bottom-full mb-2 left-1/2 transform -translate-x-1/2 opacity-0 group-hover:opacity-100 transition-opacity pointer-events-none">
                    <div className="bg-slate-800 border-2 border-purple-400 rounded-lg px-3 py-2 text-sm whitespace-nowrap shadow-xl">
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
              }} disabled={isComputing || treeNodes.length === 0} className="px-6 py-3 bg-purple-500/30 text-white rounded-lg font-semibold hover:bg-purple-500/50 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                <Sparkles className="w-5 h-5" />Relayout
              </button>
              <button onClick={reset} disabled={isComputing} className="px-6 py-3 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
                <RotateCcw className="w-5 h-5" />Reset
              </button>
              <button onClick={stop} disabled={!isComputing} className="px-6 py-3 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2">
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
          <div className="mt-6 bg-white/10 backdrop-blur-lg rounded-2xl shadow-2xl p-6 border border-white/20">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-xl font-semibold text-white">Optimization Metrics</h3>
              
              {/* Metric Toggles */}
              <div className="flex gap-2">
                {Object.keys(metricDefinitions).map(key => (
                  <button
                    key={key}
                    onClick={() => setVisibleMetrics(prev => ({ ...prev, [key]: !prev[key as keyof VisibleMetrics] }))}
                    className={`px-3 py-1 rounded-lg text-sm font-medium transition-all ${
                      visibleMetrics[key as keyof VisibleMetrics]
                        ? 'bg-purple-500 text-white'
                        : 'bg-white/10 text-purple-300 hover:bg-white/20'
                    }`}
                  >
                    {metricDefinitions[key].label}
                  </button>
                ))}
              </div>
            </div>

            {metricsHistory.length > 0 ? (
              <>
                {/* Plots Grid - Memoized to prevent re-render on mouse move */}
                {metricsCharts}

                {/* Node Count Line */}
                {Object.values(visibleMetrics).filter(Boolean).length > 0 && (
                  <div className="mt-4 text-sm text-purple-300 text-center">
                    Total Molecules: {treeNodes.length}
                  </div>
                )}
              </>
            ) : (
              <div className="text-purple-300 text-center py-8">
                Waiting for metrics data...
              </div>
            )}
          </div>
        )}

      <div className="absolute top-10 left-10 text-white">
        <svg version="1.1" id="Layer_1" height="60px" viewBox="0 0 40 40">
          <g>
	          <rect x="1.73" y="0.01" fill="#FFFFFF" width="34.19" height="34.19"/>
	          <path fill="#1E59AE" d="M35.92,0.01v17.53H18.95V0.01H35.92z M15.88,21.82c-1.12-0.07-1.72-0.78-1.79-2.1V0.01h-0.76v19.73
		c0.09,1.72,1,2.75,2.53,2.84h15.28l-4.83,5l-11.79,0c-3.04-0.36-6.22-2.93-6.14-6.98V0.01H7.6V20.6c-0.09,4.49,3.45,7.34,6.86,7.75
		h11.09l-4.59,4.75h-6.68C9.71,32.93,3.19,29.44,2.99,21.13V0.01H0.05v37.3h35.87V17.62l-4.05,4.19L15.88,21.82z"/>
          </g>
        </svg>
      </div>
      <div className="mt-8 pt-6 border-t border-purple-400/30 text-center text-purple-300 text-sm">
        <p>This work was performed under the auspices of the U.S. Department of Energy
        by Lawrence Livermore National Laboratory (LLNL) under Contract DE-AC52-07NA27344
        (LLNL-CODE-2006345).</p>
      </div>
      </div>
    </div>

    {sidebarOpen && (
      <div className="w-96 bg-slate-900 border-l-2 border-purple-400 shadow-2xl flex-shrink-0 sticky top-0 h-screen overflow-hidden animate-slideInSidebar">
        <div className="flex flex-col h-full">
          <div className="flex items-center justify-between p-4 border-b border-purple-400/30">
            <h3 className="text-lg font-semibold text-white">Reasoning</h3>
            <div className="flex items-center gap-2">
              <div className="relative">
                <button 
                  onClick={(e) => { e.stopPropagation(); setSourceFilterOpen(!sourceFilterOpen); }}
                  onMouseDown={(e) => e.stopPropagation()}
                  className="px-3 py-1 bg-purple-500/30 text-white rounded-lg text-xs font-semibold hover:bg-purple-500/50 transition-all flex items-center gap-2"
                >
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
                  </svg>
                  Filter
                  <svg className={`w-3 h-3 transition-transform ${sourceFilterOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                  </svg>
                </button>
                {sourceFilterOpen && (
                  <div className="absolute top-full mt-2 right-0 bg-slate-800 border-2 border-purple-400 rounded-lg shadow-2xl py-2 min-w-48 z-[100]" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                    <div className="px-3 py-2 border-b border-purple-400/30 text-xs font-semibold text-purple-300">
                      Message Sources
                    </div>
                    {Object.keys(visibleSources).map(source => (
                      <label key={source} className="flex items-center gap-2 px-3 py-2 hover:bg-purple-600/30 cursor-pointer transition-colors">
                        <input
                          type="checkbox"
                          checked={visibleSources[source]}
                          onChange={() => setVisibleSources(prev => ({ ...prev, [source]: !prev[source] }))}
                          className="w-4 h-4 rounded border-purple-400/50 bg-white/20 text-purple-600 focus:ring-purple-500 focus:ring-offset-0"
                        />
                        <span className="text-sm text-white">{source}</span>
                        <span className="ml-auto text-xs text-purple-400">
                          ({sidebarMessages.filter(m => m.source === source).length})
                        </span>
                      </label>
                    ))}
                    <div className="px-3 py-2 border-t border-purple-400/30 flex gap-2">
                      <button
                        onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: true }), {}))}
                        className="flex-1 px-2 py-1 bg-purple-500/30 text-white rounded text-xs font-semibold hover:bg-purple-500/50 transition-all"
                      >
                        All
                      </button>
                      <button
                        onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: false }), {}))}
                        className="flex-1 px-2 py-1 bg-purple-500/30 text-white rounded text-xs font-semibold hover:bg-purple-500/50 transition-all"
                      >
                        None
                      </button>
                    </div>
                  </div>
                )}
              </div>
              <button onClick={() => setSidebarOpen(false)} className="text-purple-300 hover:text-white transition-colors">
                <X className="w-5 h-5" />
              </button>
            </div>
          </div>

          <div className="flex-1 overflow-y-auto p-4 space-y-3" ref={sidebarRef} style={{ overflowX: 'hidden' }}>
            {sidebarMessages.filter(msg => visibleSources[msg.source]).length === 0 ? (
              <div className="text-center text-purple-400 mt-8">
                <svg className="w-12 h-12 mx-auto mb-2 opacity-50" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z" />
                </svg>
                <p className="text-sm">
                  {sidebarMessages.length === 0 ? 'No messages yet' : 'No messages match the selected filters'}
                </p>
              </div>
            ) : (
              sidebarMessages.filter(msg => visibleSources[msg.source]).map((msg, idx) => (
                <div key={msg.id} className="bg-white/5 rounded-lg p-4 border border-purple-400/30 animate-slideIn opacity-0" style={{ animationDelay: `${idx * 50}ms`, animationFillMode: 'forwards' }}>
                  <div className="flex items-center justify-between mb-2">
                    <div className="text-xs text-purple-400">
                      {new Date(msg.timestamp).toLocaleTimeString()}
                    </div>
                    <div className="px-2 py-0.5 bg-purple-500/30 rounded text-xs font-semibold text-purple-200">
                      {msg.source}
                    </div>
                  </div>
                  <div className="text-sm text-purple-100">
                    <MarkdownText text={msg.message} />
                  </div>
                  {msg.smiles && (
                    <div className="mt-3 bg-white/50 rounded-lg p-2 flex justify-center">
                      <MoleculeSVG smiles={msg.smiles} height={120} rdkitModule={rdkitModule} />
                    </div>
                  )}
                </div>
              ))
            )}
          </div>
        </div>
      </div>
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
                return prev.filter(n => n.y <= contextMenu.node!.y);
              });
              sendMessageToServer("optimize-from", nodeId);
            }}  className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
            <RotateCcw className="w-4 h-4" />
            Restart from here
            </button>
            <button 
              onClick={() => copyToClipboard(contextMenu.node!.smiles, 'smiles')} 
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
            </>
          )}

          { (problemType === "retrosynthesis" && !hasDescendants(contextMenu.node.id, treeNodes)) && (
            <button onClick={() => {
              sendMessageToServer("compute-reaction-from", contextMenu.node!.id);
            }} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <TestTubeDiagonal className="w-4 h-4" />
              How do I make this?
            </button>
            ) }

          { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) && (
            <>
            <button onClick={() => {sendMessageToServer("recompute-reaction", contextMenu.node!.id);}} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <RefreshCw className="w-4 h-4" />Find Another Reaction
            </button>
            </>
          )}
          { (problemType === "retrosynthesis" && !isRootNode(contextMenu.node.id, treeNodes)) && (
            <>
            <button onClick={() => {sendMessageToServer("recompute-parent-reaction", contextMenu.node!.id);}} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2">
              <Network className="w-4 h-4" />Substitute Molecule
            </button>
            </>
          )}

          { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) && (
          <button onClick={() => handleCustomQuery(contextMenu.node!)} className="w-full px-4 py-2 text-left text-sm text-white hover:bg-purple-600/50 transition-colors flex items-center gap-2 border-t border-purple-400/30">
            <Send className="w-4 h-4" />
            { (problemType === "retrosynthesis" && hasDescendants(contextMenu.node.id, treeNodes)) ? (<>Find Another Reaction with Custom Prompt...</>) : (<>Custom Query...</>) }
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
              <button onClick={() => { savePrompts('', ''); setProblemType('retrosynthesis'); }} className="px-4 py-3 bg-white/10 text-purple-200 rounded-lg font-medium hover:bg-white/20 transition-all flex items-center gap-2">
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
    </div>
  );
};

export default ChemistryTool;
