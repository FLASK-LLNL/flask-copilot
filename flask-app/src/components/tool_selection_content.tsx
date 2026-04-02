import React from 'react';
import { AlertTriangle } from 'lucide-react';
import { SelectableTool } from '../types';

interface ToolSelectionContentProps {
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onSelectionChange: (selectedIds: number[]) => void;
}

interface MCPToolGroup {
  key: string;
  server: string;
  executionScope: 'backend' | 'local';
  items: SelectableTool[];
}

const selectableToolName = (item: SelectableTool): string =>
  (
    item.tool_name ||
    item.tool_server.names?.[0] ||
    item.tool_server.identifier ||
    item.tool_server.server ||
    ''
  ).trim();

const groupMcpTools = (items: SelectableTool[]): MCPToolGroup[] => {
  const groups = new Map<string, MCPToolGroup>();

  items.forEach((item) => {
    const executionScope = item.tool_server.executionScope === 'local' ? 'local' : 'backend';
    const server = item.tool_server.server || item.tool_server.identifier || 'MCP Server';
    const key = `${executionScope}:${server}`;
    const existing = groups.get(key);

    if (existing) {
      existing.items.push(item);
      return;
    }

    groups.set(key, {
      key,
      server,
      executionScope,
      items: [item],
    });
  });

  return Array.from(groups.values());
};

export const ToolSelectionContent: React.FC<ToolSelectionContentProps> = ({
  availableToolsMap,
  selectedTools,
  onSelectionChange,
}) => {
  const selectedToolIds = React.useMemo(() => new Set(selectedTools), [selectedTools]);
  const availableToolById = React.useMemo(
    () => new Map(availableToolsMap.map((tool) => [tool.id, tool])),
    [availableToolsMap]
  );
  const backendMcpGroups = groupMcpTools(
    availableToolsMap.filter(
      (tool) => tool.tool_server.kind !== 'builtin' && tool.tool_server.executionScope !== 'local'
    )
  );
  const localMcpGroups = groupMcpTools(
    availableToolsMap.filter(
      (tool) => tool.tool_server.kind !== 'builtin' && tool.tool_server.executionScope === 'local'
    )
  );
  const builtinTools = availableToolsMap.filter((tool) => tool.tool_server.kind === 'builtin');

  const duplicateAvailableToolNames = React.useMemo(() => {
    const duplicateMap = new Map<string, Set<string>>();
    const nameCounts = new Map<string, number>();

    availableToolsMap.forEach((item) => {
      const duplicateKey = selectableToolName(item);
      if (!duplicateKey) {
        return;
      }
      nameCounts.set(duplicateKey, (nameCounts.get(duplicateKey) || 0) + 1);
    });

    availableToolsMap.forEach((item) => {
      const duplicateKey = selectableToolName(item);
      if (!duplicateKey || (nameCounts.get(duplicateKey) || 0) < 2) {
        return;
      }
      const serverLabel = item.tool_server.server || item.tool_server.identifier || 'MCP Server';
      const scopeLabel = item.tool_server.executionScope === 'local' ? 'local' : 'backend';
      const servers = duplicateMap.get(duplicateKey) || new Set<string>();
      servers.add(`${scopeLabel}:${serverLabel}`);
      duplicateMap.set(duplicateKey, servers);
    });

    return Array.from(duplicateMap.entries())
      .filter(([, servers]) => servers.size > 1)
      .map(([name, servers]) => ({
        name,
        servers: Array.from(servers).sort(),
      }))
      .sort((left, right) => left.name.localeCompare(right.name));
  }, [availableToolsMap]);

  const maxSelectableCount = React.useMemo(() => {
    const uniqueNames = new Set<string>();
    availableToolsMap.forEach((item) => {
      const name = selectableToolName(item);
      if (name) {
        uniqueNames.add(name);
      }
    });
    return uniqueNames.size;
  }, [availableToolsMap]);

  const applyExclusiveSelection = React.useCallback(
    (nextSelectedIds: number[]) => {
      const uniqueSelection: number[] = [];
      const seenNames = new Set<string>();

      nextSelectedIds.forEach((toolId) => {
        const item = availableToolById.get(toolId);
        if (!item) {
          return;
        }

        const name = selectableToolName(item);
        if (!name || seenNames.has(name)) {
          return;
        }

        seenNames.add(name);
        uniqueSelection.push(toolId);
      });

      onSelectionChange(uniqueSelection);
    },
    [availableToolById, onSelectionChange]
  );

  const toggleToolSelection = (toolId: number): void => {
    if (selectedToolIds.has(toolId)) {
      onSelectionChange(selectedTools.filter((id: number) => id !== toolId));
      return;
    }

    const item = availableToolById.get(toolId);
    const name = item ? selectableToolName(item) : '';
    const baseSelection = selectedTools.filter((id) => {
      const selectedItem = availableToolById.get(id);
      return !selectedItem || selectableToolName(selectedItem) !== name;
    });
    applyExclusiveSelection([...baseSelection, toolId]);
  };

  const setServerSelection = (items: SelectableTool[], selected: boolean): void => {
    const itemIds = items.map((item) => item.id);

    if (selected) {
      const itemNames = new Set(items.map((item) => selectableToolName(item)).filter(Boolean));
      const baseSelection = selectedTools.filter((id) => {
        const selectedItem = availableToolById.get(id);
        return !selectedItem || !itemNames.has(selectableToolName(selectedItem));
      });
      applyExclusiveSelection([...baseSelection, ...itemIds]);
      return;
    }

    onSelectionChange(selectedTools.filter((id) => !itemIds.includes(id)));
  };

  const handleSelectAll = (): void => {
    applyExclusiveSelection(availableToolsMap.map((tool) => tool.id));
  };

  const handleClearAll = (): void => {
    onSelectionChange([]);
  };

  const renderBuiltinSection = () => {
    if (builtinTools.length === 0) {
      return (
        <div className="bg-secondary/20 rounded-lg p-4">
          <div className="text-sm font-medium text-secondary">No built-in tools available</div>
          <div className="text-xs text-tertiary mt-1">
            This backend has not registered any direct Python tools.
          </div>
        </div>
      );
    }

    return (
      <div className="space-y-3">
        <div>
          <div className="text-sm font-medium text-primary">Built-in Tools</div>
          <div className="text-xs text-secondary mt-1">
            Backend Python functions exposed directly to the agent without MCP.
          </div>
        </div>
        <div className="tool-list space-y-2 custom-scrollbar">
          {builtinTools.map((item) => (
            <label key={item.id} className="tool-list-item">
              <input
                type="checkbox"
                checked={selectedToolIds.has(item.id)}
                onChange={() => toggleToolSelection(item.id)}
                className="form-checkbox"
              />
              <div className="flex-1">
                <div className="tool-list-item-label">
                  {item.tool_server.server || item.tool_name || 'Built-in Tool'}
                </div>
                {item.tool_name && item.tool_name !== item.tool_server.server && (
                  <div className="text-xs text-tertiary mt-0.5">{item.tool_name}</div>
                )}
                {item.tool_description && (
                  <div className="text-xs text-secondary mt-1">{item.tool_description}</div>
                )}
              </div>
            </label>
          ))}
        </div>
      </div>
    );
  };

  const renderMcpSection = (
    title: string,
    subtitle: string,
    groups: MCPToolGroup[],
    emptyTitle: string,
    emptySubtitle: string
  ) => {
    if (groups.length === 0) {
      return (
        <div className="bg-secondary/20 rounded-lg p-4">
          <div className="text-sm font-medium text-secondary">{emptyTitle}</div>
          <div className="text-xs text-tertiary mt-1">{emptySubtitle}</div>
        </div>
      );
    }

    return (
      <div className="space-y-3">
        <div>
          <div className="text-sm font-medium text-primary">{title}</div>
          <div className="text-xs text-secondary mt-1">{subtitle}</div>
        </div>
        <div className="space-y-3">
          {groups.map((group) => {
            const selectedCount = group.items.filter((item) => selectedToolIds.has(item.id)).length;
            const scopeBadge = group.executionScope === 'local' ? 'MCP local' : 'MCP';
            const scopeDescription =
              group.executionScope === 'local'
                ? 'Reached from this browser and proxied over the websocket.'
                : 'Reached directly by the backend process.';

            return (
              <div key={group.key} className="glass-panel space-y-3">
                <div className="flex items-start justify-between gap-3">
                  <div>
                    <div className="text-sm font-medium text-primary">{group.server}</div>
                    <div className="text-xs text-secondary mt-1">{scopeDescription}</div>
                  </div>
                  <div className="text-right">
                    <div className="text-[10px] tracking-wide text-tertiary">{scopeBadge}</div>
                    <div className="text-xs text-secondary mt-1">
                      {selectedCount} / {group.items.length} selected
                    </div>
                  </div>
                </div>
                <div className="flex gap-2">
                  <button
                    onClick={() => setServerSelection(group.items, true)}
                    disabled={selectedCount === group.items.length}
                    className="btn btn-tertiary btn-sm"
                  >
                    Select Server
                  </button>
                  <button
                    onClick={() => setServerSelection(group.items, false)}
                    disabled={selectedCount === 0}
                    className="btn btn-tertiary btn-sm"
                  >
                    Clear Server
                  </button>
                </div>
                <div className="space-y-2">
                  {group.items.map((item) => (
                    <label key={item.id} className="tool-list-item">
                      <input
                        type="checkbox"
                        checked={selectedToolIds.has(item.id)}
                        onChange={() => toggleToolSelection(item.id)}
                        className="form-checkbox"
                      />
                      <div className="flex-1">
                        <div className="tool-list-item-label">
                          {item.tool_name || 'Unnamed tool'}
                        </div>
                        {item.tool_description && (
                          <div className="text-xs text-secondary mt-1">{item.tool_description}</div>
                        )}
                      </div>
                    </label>
                  ))}
                </div>
              </div>
            );
          })}
        </div>
      </div>
    );
  };

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <p className="text-sm text-secondary">
          Select which MCP and built-in tools to use for this workflow
        </p>
        <div className="flex gap-2">
          <button
            onClick={handleSelectAll}
            disabled={selectedTools.length === maxSelectableCount}
            className="btn btn-tertiary btn-sm"
          >
            Select All
          </button>
          <button
            onClick={handleClearAll}
            disabled={selectedTools.length === 0}
            className="btn btn-tertiary btn-sm"
          >
            Clear All
          </button>
        </div>
      </div>

      {duplicateAvailableToolNames.length > 0 && (
        <div className="glass-panel" style={{ border: '1px solid rgba(251, 191, 36, 0.45)' }}>
          <div className="flex items-start gap-2">
            <AlertTriangle className="w-4 h-4 mt-0.5 text-yellow-300 flex-shrink-0" />
            <div>
              <div className="text-sm font-medium text-primary">
                Duplicate tool names require manual choice
              </div>
              <div className="text-xs text-secondary mt-1">
                Agent Framework cannot run with duplicate tool names. These names are not
                auto-selected, and selecting one will replace any currently selected tool with the
                same name.
              </div>
              <div className="text-xs text-tertiary mt-2">
                {duplicateAvailableToolNames
                  .map((entry) => `${entry.name} (${entry.servers.join(', ')})`)
                  .join('; ')}
              </div>
            </div>
          </div>
        </div>
      )}

      {availableToolsMap.length === 0 ? (
        <div className="text-center py-8 text-tertiary">
          <p className="text-sm">No tools available</p>
          <p className="text-xs mt-1">
            Connect MCP servers or enable backend built-ins to use tools
          </p>
        </div>
      ) : (
        <div className="space-y-5" style={{ maxHeight: '400px', overflowY: 'auto' }}>
          {renderBuiltinSection()}
          {renderMcpSection(
            'Backend MCP Tool Servers',
            'External MCP servers the backend can reach directly. Select only the tools you want to expose from each server.',
            backendMcpGroups,
            'No backend MCP tool servers available',
            'Connect backend-accessible MCP servers in settings to enable them.'
          )}
          {renderMcpSection(
            'Local MCP Tool Servers',
            'MCP servers reached from this browser and proxied to the orchestrator over the websocket. Select only the tools you want to expose from each server.',
            localMcpGroups,
            'No local MCP tool servers available',
            'Add local MCP servers in settings to expose browser-reachable endpoints such as 127.0.0.1.'
          )}
        </div>
      )}

      <div className="bg-secondary/30 rounded-lg p-3 text-sm text-secondary">
        <strong>
          Selected: {selectedTools.length} of {availableToolsMap.length} tools
        </strong>
        <p className="text-xs text-tertiary mt-1">
          Selected tools will be available to the AI agent during workflow execution
        </p>
      </div>
    </div>
  );
};
