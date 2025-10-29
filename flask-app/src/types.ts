// TypeScript interfaces and types
import { RDKitModule } from '@rdkit/rdkit';
import { NODE_STYLES } from "./constants";

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
  name?: string;
  description?: string;
}

export interface WebSocketMessageToServer {
  action?: string;
  smiles?: string;
  problemType?: string;
  nodeId?: string;
  query?: string;
}

// Messages received from backend
export interface WebSocketMessage {
  type: string;
  
  node?: TreeNode;
  edge?: Edge;
  message?: SidebarMessage;
  tools?: Tool[];

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
  setOffset: (position: Position) => void;
  zoom: number;
  setZoom: (n: number) => void;
}

export interface MoleculeGraphProps extends MoleculeGraphState {
    nodes: TreeNode[];
    edges: Edge[];
    ctx: ContextMenuState;
    autoZoom: boolean;
    setAutoZoom: (v: boolean) => void;
    handleNodeClick: (e: React.MouseEvent<HTMLDivElement>, node: TreeNode) => void;
    rdkitModule: RDKitModule | null;
}
