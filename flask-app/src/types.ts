// TypeScript interfaces and types
import { RDKitModule } from '@rdkit/rdkit';
import { NODE_STYLES } from "./constants";
import { Dispatch, SetStateAction } from 'react';

// Class that only exists for UI preview purposes
export interface PathwayStep {
  smiles: string[];
  label: string[];
  parents: number[];
}

export interface ReactionAlternative {
  id: string;
  name: string;
  type: 'exact' | 'template' | 'ai';
  status: 'active' | 'available' | 'computing';
  hoverInfo: string;
  disabled?: boolean;
  disabledReason?: string;
  pathway: PathwayStep[];  // For UI preview
}

export interface Reaction {
  id: string;
  label?: string;
  hoverInfo: string;
  highlight: keyof typeof NODE_STYLES;
  alternatives?: ReactionAlternative[];
  templatesSearched: boolean;  // Whether to show the "Search Templates" button
  mappedReaction?: RdkitjsReactionPayload;
}

export interface TreeNode {
  id: string;
  smiles: string;
  label: string;
  hoverInfo: string;
  level: number;
  parentId: string | null;
  cost?: number | null;
  bandgap?: number | null;
  density?: number | null;
  yield?: number | null;
  sascore?: number | null;
  purchasable?: boolean | null;
  x: number;
  y: number;
  highlight?: keyof typeof NODE_STYLES;
  reaction?: Reaction;
}

export interface Edge {
  id: string;
  fromNode: string;
  toNode: string;
  status?: 'complete' | 'computing';
  label?: string;
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

export interface RunSettings {
  moleculeName: MoleculeNameFormat;
  promptDebugging: boolean;
}

import type {
  OrchestratorSettings,
  SidebarMessage,
  SidebarState,
  SidebarProps,
  VisibleSources,
  MarkdownTextProps,
  ToolServer,
  MoleculeNameFormat
} from 'lcc-conductor';

// Optimization customization options
export interface OptimizationCustomization {
  enableConstraints?: boolean;
  molecularSimilarity?: number;
  diversityPenalty?: number;
  explorationRate?: number;
  additionalConstraints?: string[];  // Array of constraint types
}

export interface ConstraintOption {
  value: string;
  label: string;
  description: string;
}

export interface WebSocketMessageToServer {
  action?: string;
  smiles?: string;
  problemType?: string;
  nodeId?: string;
  query?: string;
  experimentContext?: string;
  enabledTools?: ToolMap;

  // Lead molecule optimization
  propertyType?: string;
  customPropertyName?: string;
  customPropertyDesc?: string;
  customPropertyAscending?: boolean;
  xpos?: number;

  // Settings
  runSettings?: RunSettings;

  // Optimization customization
  customization?: OptimizationCustomization;

  // Custom problem
  systemPrompt?: string;
  userPrompt?: string;

  // Retrosynthesis
  alternativeId?: string;

  // Prompt debugging
  prompt?: string;
  metadata?: any;
}

// Messages received from backend
export interface WebSocketMessage {
  type: string;

  node?: TreeNode;
  edge?: Edge;
  message?: SidebarMessage;
  tools?: Tool[];
  experimentContext?: string;
  orchestratorSettings?: OrchestratorSettings;

  withNode?: boolean;
  username?: string;

  alternatives?: ReactionAlternative[];

  // Prompt debugging
  prompt?: string;
  metadata?: any;
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
  sascore: boolean;
  density: boolean;
  yield: boolean;
}

export interface ContextMenuState {
  node: TreeNode | null;
  isReaction: boolean;
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
  highlightAtomIdxs?: number[];
  highlightRgb?: RGB;
  highlightAlpha?: number;
}

// rdkit.js payload types for reaction highlighting
export type RGB = [number, number, number];

export interface RdkitjsMolPayload {
  smiles: string;
  highlight_atom_idxs: number[];
  highlight_atom_mapnums: number[];
}

export interface RdkitjsReactionPayload {
  reactants: RdkitjsMolPayload[];
  products: RdkitjsMolPayload[];
  main_product_index: number;
  highlight_rgb: RGB;
  highlight_alpha: number;
  reactant_mcs_smarts?: (string | null)[];
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
    handleReactionClick: (e: React.MouseEvent<HTMLDivElement>, node: TreeNode) => void;
    handleReactionCardClick: (node: TreeNode) => void;
    selectedReactionNodeId?: string;
    reactionSidebarOpen: boolean;
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
  propertyType?: string;
  customPropertyName?: string;
  customPropertyDesc?: string;
  customPropertyAscending?: boolean;
  customization?: OptimizationCustomization;
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
