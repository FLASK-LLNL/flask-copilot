import { TreeNode, Edge, Position } from './types';
import { BOX_GAP, BOX_WIDTH } from './constants';

export function findAllDescendants(nodeId: string, nodes: TreeNode[]): Set<string> {
    const descendants = new Set<string>();
    const findDescendants = (id: string): void => {
      nodes.forEach(n => {
        if (n.parentId === id && !descendants.has(n.id)) {
          descendants.add(n.id);
          findDescendants(n.id);
        }
      });
    };
    findDescendants(nodeId);
    return descendants;
}

export function hasDescendants(nodeId: string, nodes: TreeNode[]): boolean {
    return nodes.some(n => n.parentId === nodeId);
}

export function isRootNode(nodeId: string, nodes: TreeNode[]): boolean {
    const node = nodes.find(n => n.id === nodeId);
    return !node?.parentId;
}

  // Relayouts the molecule graph for better visibility (assumes tree)
export function relayoutTree(treeNodes: TreeNode[], edges: Edge[]): [TreeNode[], Edge[]] {
    if (treeNodes.length === 0) return [treeNodes, edges];

    // Build parent-to-children map
    const childrenMap = new Map<string, TreeNode[]>();
    treeNodes.forEach(node => {
      if (node.parentId) {
        if (!childrenMap.has(node.parentId)) {
          childrenMap.set(node.parentId, []);
        }
        childrenMap.get(node.parentId)!.push(node);
      }
    });

    // Find root node(s)
    const roots = treeNodes.filter(n => n.parentId === null || !treeNodes.find(t => t.id === n.parentId));

    const levelGap = BOX_WIDTH + BOX_GAP;
    const nodeSpacing = 150;
    const newPositions = new Map<string, Position>();

    // First pass: calculate subtree sizes (number of leaf descendants)
    const subtreeSizes = new Map<string, number>();
    const calculateSubtreeSize = (nodeId: string): number => {
      const children = childrenMap.get(nodeId) || [];
      if (children.length === 0) {
        subtreeSizes.set(nodeId, 1);
        return 1;
      }
      const size = children.reduce((sum, child) => sum + calculateSubtreeSize(child.id), 0);
      subtreeSizes.set(nodeId, size);
      return size;
    };

    roots.forEach(root => calculateSubtreeSize(root.id));

    // Second pass: assign positions based on subtree sizes
    const assignPositions = (nodeId: string, level: number, startY: number): number => {
      const node = treeNodes.find(n => n.id === nodeId);
      if (!node) return startY;

      const children = childrenMap.get(nodeId) || [];

      if (children.length === 0) {
        // Leaf node - place at next available Y
        newPositions.set(nodeId, {
          x: 100 + level * levelGap,
          y: startY
        });
        return startY + nodeSpacing;
      }

      // Internal node - place children first, then center parent
      let currentY = startY;
      const childPositions: number[] = [];

      children.forEach(child => {
        currentY = assignPositions(child.id, level + 1, currentY);
        childPositions.push(newPositions.get(child.id)!.y);
      });

      // Center parent among its children
      const avgChildY = childPositions.reduce((sum, y) => sum + y, 0) / childPositions.length;
      newPositions.set(nodeId, {
        x: 100 + level * levelGap,
        y: avgChildY
      });

      return currentY;
    };

    let currentY = 100;
    roots.forEach(root => {
      currentY = assignPositions(root.id, 0, currentY);
    });

    // Update nodes with new positions
    const updatedNodes = treeNodes.map(node => {
      const pos = newPositions.get(node.id);
      if (pos) {
        return { ...node, x: pos.x, y: pos.y };
      }
      return node;
    });

    // Create node map for edges
    const nodeMap: { [key: string]: TreeNode } = {};
    updatedNodes.forEach(n => {
      nodeMap[n.id] = n;
    });

    // Update all edges with new positions
    const updatedEdges = edges.map(e => ({
      ...e,
    }));

    return [updatedNodes, updatedEdges];
}
