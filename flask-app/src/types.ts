// TypeScript interfaces and types
import { RDKitModule } from '@rdkit/rdkit';
import { NODE_STYLES } from "./constants";
import { Dispatch, SetStateAction } from 'react';

export interface TreeNode {
  id: string;
  smiles: string;
  label: string;
  hoverInfo: string;
  level: number;
  parentId: string | null;
  cost?: number;
  bandgap?: number;
  density?: number;
  yield?: number;
  x: number;
  y: number;
  highlight?: keyof typeof NODE_STYLES;
}

export interface Edge {
  id: string;
  fromNode: string;
  toNode: string;
  reactionType?: string;
  status?: 'complete' | 'computing';
  label?: string;
}


export interface SidebarMessage {
  id: number;
  timestamp: string;
  message: string;
  smiles: string | null;
  source: string;
}

export interface Tool {
  server?: string;
  names?: string[];
  description?: string;
}

export interface SelectableTool {
  id: number;
  tool_server: Tool;
}

export interface ToolMap {
  selectedIds?: number[];
  selectedTools?: SelectableTool[];
}

export interface ProfileSettings {
  backend: string;
  useCustomUrl: boolean;
  customUrl?: string;
  model: string;
  apiKey: string;
}

export interface WebSocketMessageToServer {
  action?: string;
  smiles?: string;
  problemType?: string;
  nodeId?: string;
  query?: string;
  experimentContext?: string;
  enabledTools?: ToolMap;
}

// Messages received from backend
export interface WebSocketMessage {
  type: string;

  node?: TreeNode;
  edge?: Edge;
  message?: SidebarMessage;
  tools?: Tool[];
  experimentContext?: string;

  withNode?: boolean;
}


export interface MetricDefinition {
  label: string;
  color: string;
  calculate: (nodes: TreeNode[]) => number;
}

export interface MetricDefinitions {
  [key: string]: MetricDefinition;
}

export interface VisibleMetrics {
  cost: boolean;
  bandgap: boolean;
  density: boolean;
  yield: boolean;
}

export interface VisibleSources {
  [key: string]: boolean;
}

export interface ContextMenuState {
  node: TreeNode | null;
  x: number;
  y: number;
}

export interface Position {
  x: number;
  y: number;
}

export interface MetricHistoryItem {
  step: number;
  nodeCount: number;
  [key: string]: number;
}

export interface MoleculeSVGProps {
  smiles: string;
  height?: number;
  rdkitModule: RDKitModule | null;
}

export interface MarkdownTextProps {
  text: string;
}

export interface MoleculeGraphState {
  offset: Position;
  setOffset: Dispatch<SetStateAction<Position>>;
  zoom: number;
  setZoom: Dispatch<SetStateAction<number>>;
}

export interface MoleculeGraphProps extends MoleculeGraphState {
    nodes: TreeNode[];
    edges: Edge[];
    ctx: ContextMenuState;
    autoZoom: boolean;
    setAutoZoom: Dispatch<SetStateAction<boolean>>;
    handleNodeClick: (e: React.MouseEvent<HTMLDivElement>, node: TreeNode) => void;
    rdkitModule: RDKitModule | null;
}

export interface SidebarState {
    messages: SidebarMessage[];
    setMessages: Dispatch<SetStateAction<SidebarMessage[]>>;
    sourceFilterOpen: boolean;
    setSourceFilterOpen: Dispatch<SetStateAction<boolean>>;
    visibleSources: VisibleSources;
    setVisibleSources: Dispatch<SetStateAction<VisibleSources>>;
}

export interface SidebarProps extends SidebarState {
    // General state from app
    setSidebarOpen: Dispatch<SetStateAction<boolean>>;
    rdkitModule: RDKitModule | null;
}

export interface MetricsDashboardState {
    metricsHistory: MetricHistoryItem[];
    setMetricsHistory: Dispatch<SetStateAction<MetricHistoryItem[]>>;
    visibleMetrics: VisibleMetrics;
    setVisibleMetrics: Dispatch<SetStateAction<VisibleMetrics>>;
}

export interface MetricsDashboardProps extends MetricsDashboardState {
    // General state from app
    treeNodes: TreeNode[];
}

// Experiment types
export interface Experiment {
  id: string;
  name: string;
  createdAt: string;
  lastModified: string;
  isRunning?: boolean;  // Track if experiment is currently computing

  // System state
  smiles?: string;
  problemType?: string;
  problemName?: string;
  systemPrompt?: string;
  problemPrompt?: string;
  treeNodes?: TreeNode[];
  edges?: Edge[];
  metricsHistory?: MetricHistoryItem[];
  visibleMetrics?: VisibleMetrics;
  graphState?: MoleculeGraphState;
  autoZoom?: boolean;
  sidebarState?: SidebarState;

  // Experiment state
  experimentContext?: string;
}

export interface Project {
  id: string;
  name: string;
  createdAt: string;
  lastModified: string;
  experiments: Experiment[];
}

export interface ProjectSelection {
  projectId: string | null;
  experimentId: string | null;
}
