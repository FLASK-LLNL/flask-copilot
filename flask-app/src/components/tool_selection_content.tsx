import React from 'react';
import { SelectableTool } from '../types';

interface ToolSelectionContentProps {
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onSelectionChange: (selectedIds: number[]) => void;
}

export const ToolSelectionContent: React.FC<ToolSelectionContentProps> = ({
  availableToolsMap,
  selectedTools,
  onSelectionChange,
}) => {
  const mcpTools = availableToolsMap.filter((tool) => tool.tool_server.kind !== 'builtin');
  const builtinTools = availableToolsMap.filter((tool) => tool.tool_server.kind === 'builtin');

  const toggleToolSelection = (toolId: number): void => {
    const newSelection = selectedTools.includes(toolId)
      ? selectedTools.filter((id: number) => id !== toolId)
      : [...selectedTools, toolId];
    onSelectionChange(newSelection);
  };

  const handleSelectAll = (): void => {
    onSelectionChange(availableToolsMap.map((tool) => tool.id));
  };

  const handleClearAll = (): void => {
    onSelectionChange([]);
  };

  const renderToolSection = (
    title: string,
    subtitle: string,
    items: SelectableTool[],
    emptyTitle: string,
    emptySubtitle: string
  ) => {
    if (items.length === 0) {
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
        <div className="tool-list space-y-2 custom-scrollbar">
          {items.map((item: SelectableTool) => (
            <label key={item.id} className="tool-list-item">
              <input
                type="checkbox"
                checked={selectedTools.includes(item.id)}
                onChange={() => toggleToolSelection(item.id)}
                className="form-checkbox"
              />
              <div className="flex-1">
                <span className="tool-list-item-label">{item.tool_server.server}</span>
                {item.tool_server.names && item.tool_server.names.length > 0 && (
                  <div className="text-xs text-tertiary mt-0.5">
                    {item.tool_server.names.join(', ')}
                  </div>
                )}
                {item.tool_server.description && (
                  <div className="text-xs text-secondary mt-1">{item.tool_server.description}</div>
                )}
              </div>
            </label>
          ))}
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
            disabled={selectedTools.length === availableToolsMap.length}
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

      {availableToolsMap.length === 0 ? (
        <div className="text-center py-8 text-tertiary">
          <p className="text-sm">No tools available</p>
          <p className="text-xs mt-1">
            Connect MCP servers or enable backend built-ins to use tools
          </p>
        </div>
      ) : (
        <div className="space-y-5" style={{ maxHeight: '400px', overflowY: 'auto' }}>
          {renderToolSection(
            'Built-in Tools',
            'Backend Python functions exposed directly to the agent without MCP.',
            builtinTools,
            'No built-in tools available',
            'This backend has not registered any direct Python tools.'
          )}
          {renderToolSection(
            'MCP Tool Servers',
            'External tool servers connected through MCP.',
            mcpTools,
            'No MCP tool servers available',
            'Connect MCP servers in settings to enable external tools.'
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
