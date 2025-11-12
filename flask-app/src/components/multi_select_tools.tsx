import React from 'react';

import { Tool } from './types';

export interface SelectableTool {
  id: number;
  tool_server: Tool;
}

interface MultiSelectToolModalProps {
  isOpen: boolean;
  onClose: () => void;
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onSelectionChange: (selectedIds: number[]) => void;
  onConfirm?: (selectedIds: number[], selectedItemsData: SelectableItem[]) => void | Promise<void>;
  title?: string;
}

export const MultiSelectToolModal: React.FC<MultiSelectToolModalProps> = ({
  isOpen,
  onClose,
  availableToolsMap,
  selectedTools,
  onSelectionChange,
  onConfirm,
  title = "Select Tools"
}) => {
  if (!isOpen) return null;

  const toggleToolSelection = (toolId: number): void => {
    const newSelection = selectedTools.includes(toolId)
      ? selectedTools.filter((id: number) => id !== toolId)
      : [...selectedTools, toolId];
    onSelectionChange(newSelection);
  };

  const handleDone = async (): Promise<void> => {
    // Get the full item data for selected tools
    const selectedToolsData = availableToolsMap.filter(tool =>
      selectedTools.includes(tool.id)
    );

    // Call the onConfirm callback if provided
    if (onConfirm) {
      await onConfirm(selectedTools, selectedToolsData);
    }

    onClose();
  };

  const handleClearAll = (): void => {
    onSelectionChange([]);
  };

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50">
      <div className="bg-slate-800 rounded-lg p-6 max-w-md w-full mx-4 border-2 border-purple-400/50">
        <h3 className="text-xl font-semibold text-purple-200 mb-4">{title}</h3>

        <div className="max-h-96 overflow-y-auto space-y-2 mb-4">
          {availableToolsMap.map((item: SelectableTool) => (
            <label
              key={item.id}
              className="flex items-center gap-3 p-3 bg-white/10 rounded-lg hover:bg-white/20 cursor-pointer transition-colors"
            >
              <input
                type="checkbox"
                checked={selectedTools.includes(item.id)}
                onChange={() => toggleToolSelection(item.id)}
                className="w-4 h-4 accent-purple-500"
              />
              <span className="text-white">{item.tool_server.server}: [{item.tool_server.names.join(", ")}]</span>
            </label>
          ))}
        </div>

        <div className="flex gap-3">
          <button
            onClick={handleDone}
            className="flex-1 px-4 py-2 bg-purple-500 hover:bg-purple-600 text-white rounded-lg font-medium transition-colors"
          >
            Done
          </button>
          <button
            onClick={handleClearAll}
            className="px-4 py-2 bg-white/10 hover:bg-white/20 text-purple-200 rounded-lg font-medium transition-colors"
          >
            Clear All
          </button>
        </div>
      </div>
    </div>
  );
};
