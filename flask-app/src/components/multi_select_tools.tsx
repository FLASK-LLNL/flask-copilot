import React from 'react';

import { Tool, SelectableTool } from '../types';


interface MultiSelectToolModalProps {
  isOpen: boolean;
  onClose: () => void;
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onSelectionChange: (selectedIds: number[]) => void;
  onConfirm?: (selectedIds: number[], selectedItemsData: SelectableTool[]) => void | Promise<void>;
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
    <div className="modal-overlay">
      <div className="modal-content modal-content-sm">
        <h3 className="modal-title mb-4">{title}</h3>

        <div className="tool-list space-y-2 custom-scrollbar">
          {availableToolsMap.map((item: SelectableTool) => (
            <label
              key={item.id}
              className="tool-list-item"
            >
              <input
                type="checkbox"
                checked={selectedTools.includes(item.id)}
                onChange={() => toggleToolSelection(item.id)}
                className="form-checkbox"
              />
              <span className="tool-list-item-label">{item.tool_server.server}: [{item.tool_server.names!.join(", ")}]</span>
            </label>
          ))}
        </div>

        <div className="flex gap-3 mt-4">
          <button
            onClick={handleDone}
            className="btn btn-secondary flex-1"
          >
            Done
          </button>
          <button
            onClick={handleClearAll}
            className="btn btn-tertiary"
          >
            Clear All
          </button>
        </div>
      </div>
    </div>
  );
};
