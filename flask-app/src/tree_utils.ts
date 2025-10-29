import { useCallback } from 'react';
import { TreeNode } from './types';


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
