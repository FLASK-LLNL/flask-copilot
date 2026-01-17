import React, { useState } from 'react';
import { X, Wrench, Network } from 'lucide-react';
import { OptimizationCustomization, SelectableTool } from '../types';
import { ToolSelectionContent } from './tool_selection_content';
import { OptimizationCustomizationContent } from './optimization_customization_content';

interface CombinedCustomizationModalProps {
  isOpen: boolean;
  onClose: () => void;

  // Tool selection props
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onToolSelectionChange: (selectedIds: number[]) => void;
  onToolConfirm?: (selectedIds: number[], selectedItemsData: SelectableTool[]) => void | Promise<void>;

  // Optimization customization props
  initialCustomization: OptimizationCustomization;
  onCustomizationSave: (customization: OptimizationCustomization) => void;

  // Show optimization tab only for optimization problem type
  showOptimizationTab?: boolean;
}

type TabType = 'tools' | 'optimization';

export const CombinedCustomizationModal: React.FC<CombinedCustomizationModalProps> = ({
  isOpen,
  onClose,
  availableToolsMap,
  selectedTools,
  onToolSelectionChange,
  onToolConfirm,
  initialCustomization,
  onCustomizationSave,
  showOptimizationTab = true,
}) => {
  const [activeTab, setActiveTab] = useState<TabType>('tools');
  const [pendingCustomization, setPendingCustomization] = useState<OptimizationCustomization>(initialCustomization);

  if (!isOpen) return null;

  const handleToolConfirm = async () => {
    const selectedToolsData = availableToolsMap.filter(tool =>
      selectedTools.includes(tool.id)
    );

    if (onToolConfirm) {
      await onToolConfirm(selectedTools, selectedToolsData);
    }
  };

  const handleCustomizationSave = () => {
    onCustomizationSave(pendingCustomization);
  };

  const handleApplyAndClose = async () => {
    // Apply tool selection
    await handleToolConfirm();

    // Apply customization if tab is shown
    if (showOptimizationTab) {
      handleCustomizationSave();
    }

    onClose();
  };

  return (
    <div className="modal-overlay">
      <div className="modal-content modal-content-lg">
        <div className="modal-header">
          <div>
            <h2 className="modal-title">Customize Workflow</h2>
            <p className="modal-subtitle">Configure tools and optimization parameters</p>
          </div>
          <button onClick={onClose} className="btn-icon">
            <X className="w-6 h-6" />
          </button>
        </div>

        {/* Tab Navigation */}
        <div className="flex border-b border-secondary">
          <button
            onClick={() => setActiveTab('tools')}
            className={`flex items-center gap-2 px-6 py-3 font-medium transition-colors border-b-2 ${
              activeTab === 'tools'
                ? 'border-primary text-primary'
                : 'border-transparent text-secondary hover:text-primary'
            }`}
          >
            <Network className="w-4 h-4" />
            Tool Selection
            {selectedTools.length > 0 && (
              <span className="ml-1 px-2 py-0.5 text-xs rounded-full bg-primary/20 text-primary">
                {selectedTools.length}
              </span>
            )}
          </button>

          {showOptimizationTab && (
            <button
              onClick={() => setActiveTab('optimization')}
              className={`flex items-center gap-2 px-6 py-3 font-medium transition-colors border-b-2 ${
                activeTab === 'optimization'
                  ? 'border-primary text-primary'
                  : 'border-transparent text-secondary hover:text-primary'
              }`}
            >
              <Wrench className="w-4 h-4" />
              Optimization
              {pendingCustomization.enableConstraints && (
                <span className="ml-1 px-2 py-0.5 text-xs rounded-full bg-primary/20 text-primary">
                  ON
                </span>
              )}
            </button>
          )}
        </div>

        {/* Tab Content */}
        <div className="modal-body" style={{ minHeight: '400px' }}>
          {activeTab === 'tools' && (
            <ToolSelectionContent
              availableToolsMap={availableToolsMap}
              selectedTools={selectedTools}
              onSelectionChange={onToolSelectionChange}
            />
          )}

          {activeTab === 'optimization' && showOptimizationTab && (
            <OptimizationCustomizationContent
              customization={pendingCustomization}
              onCustomizationChange={setPendingCustomization}
              showResetButton={true}
            />
          )}
        </div>

        {/* Footer with Apply/Cancel */}
        <div className="modal-footer">
          <button
            onClick={onClose}
            className="btn btn-tertiary"
          >
            Cancel
          </button>
          <button
            onClick={handleApplyAndClose}
            className="btn btn-primary flex-1"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
            </svg>
            Apply Changes
          </button>
        </div>
      </div>
    </div>
  );
};
