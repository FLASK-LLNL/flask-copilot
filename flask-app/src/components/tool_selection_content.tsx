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

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <p className="text-sm text-secondary">Select which tool servers to use for this workflow</p>
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
          <p className="text-sm">No tool servers available</p>
          <p className="text-xs mt-1">Connect MCP servers in settings to enable tools</p>
        </div>
      ) : (
        <div
          className="tool-list space-y-2 custom-scrollbar"
          style={{ maxHeight: '400px', overflowY: 'auto' }}
        >
          {availableToolsMap.map((item: SelectableTool) => (
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
