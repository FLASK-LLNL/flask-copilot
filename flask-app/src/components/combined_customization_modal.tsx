import React, { useState } from 'react';
import { BookOpen, FileText, X, Wrench, Network, FlaskConical } from 'lucide-react';
import { AttachmentUpload } from 'lc-conductor';
import type { AgentAttachment } from 'lc-conductor';
import { OptimizationCustomization, PdfReferenceMetadata, SelectableTool } from '../types';
import { ToolSelectionContent } from './tool_selection_content';
import { OptimizationCustomizationContent } from './optimization_customization_content';
import { MoleculePropertiesContent } from './molecule_properties';

interface CombinedCustomizationModalProps {
  isOpen: boolean;
  onClose: () => void;

  // Tool selection props
  availableToolsMap: SelectableTool[];
  selectedTools: number[];
  onToolSelectionChange: (selectedIds: number[]) => void;
  onToolConfirm?: (
    selectedIds: number[],
    selectedItemsData: SelectableTool[]
  ) => void | Promise<void>;

  // Optimization customization props
  initialCustomization: OptimizationCustomization;
  onCustomizationSave: (customization: OptimizationCustomization) => void;

  // Molecule properties props
  initialMoleculeName?: string;
  onMoleculeNameSave?: (moleculeName: string) => void;

  // Reference document props
  referenceDocument?: PdfReferenceMetadata | null;
  referenceUploadDisabled?: boolean;
  onReferenceDocumentSave?: (reference: AgentAttachment | null | undefined) => void | Promise<void>;

  // Show optimization tab only for optimization problem type
  showOptimizationTab?: boolean;
}

type TabType = 'tools' | 'optimization' | 'molecule' | 'references';

const sizeLabel = (sizeBytes: number): string => {
  if (sizeBytes < 1024 * 1024) {
    return `${Math.max(1, Math.round(sizeBytes / 1024))} KB`;
  }
  return `${(sizeBytes / (1024 * 1024)).toFixed(1)} MB`;
};

export const CombinedCustomizationModal: React.FC<CombinedCustomizationModalProps> = (props) => {
  if (!props.isOpen) return null;
  return <CombinedCustomizationModalContent {...props} />;
};

const CombinedCustomizationModalContent: React.FC<CombinedCustomizationModalProps> = ({
  onClose,
  availableToolsMap,
  selectedTools,
  onToolSelectionChange,
  onToolConfirm,
  initialCustomization,
  onCustomizationSave,
  initialMoleculeName = 'brand',
  onMoleculeNameSave,
  referenceDocument,
  referenceUploadDisabled = false,
  onReferenceDocumentSave,
  showOptimizationTab = true,
}) => {
  const [activeTab, setActiveTab] = useState<TabType>('tools');
  const [pendingCustomization, setPendingCustomization] =
    useState<OptimizationCustomization>(initialCustomization);
  const [pendingMoleculeName, setPendingMoleculeName] = useState<string>(initialMoleculeName);
  const [pendingReference, setPendingReference] = useState<AgentAttachment[]>([]);
  const [referenceRemoved, setReferenceRemoved] = useState(false);
  const [pendingReferenceUploaded, setPendingReferenceUploaded] = useState(false);
  const hasLiveReference =
    referenceDocument?.status === 'uploading' || referenceDocument?.status === 'available';

  const handleToolConfirm = async () => {
    const selectedToolsData = availableToolsMap.filter((tool) => selectedTools.includes(tool.id));

    if (onToolConfirm) {
      await onToolConfirm(selectedTools, selectedToolsData);
    }
  };

  const handleCustomizationSave = () => {
    onCustomizationSave(pendingCustomization);
  };

  const handleMoleculeNameSave = () => {
    if (onMoleculeNameSave) {
      onMoleculeNameSave(pendingMoleculeName);
    }
  };

  const handleReferenceDocumentSave = async () => {
    if (!onReferenceDocumentSave) {
      return;
    }
    if (pendingReference.length > 0 && !pendingReferenceUploaded) {
      await onReferenceDocumentSave(pendingReference[0]);
      return;
    }
    if (referenceRemoved) {
      await onReferenceDocumentSave(null);
      return;
    }
    await onReferenceDocumentSave(undefined);
  };

  const handleApplyAndClose = async () => {
    // Apply tool selection
    await handleToolConfirm();

    // Apply customization if tab is shown
    if (showOptimizationTab) {
      handleCustomizationSave();
    }

    // Apply molecule name preference
    handleMoleculeNameSave();

    await handleReferenceDocumentSave();

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

          <button
            onClick={() => setActiveTab('molecule')}
            className={`flex items-center gap-2 px-6 py-3 font-medium transition-colors border-b-2 ${
              activeTab === 'molecule'
                ? 'border-primary text-primary'
                : 'border-transparent text-secondary hover:text-primary'
            }`}
          >
            <FlaskConical className="w-4 h-4" />
            Molecule Display
          </button>

          <button
            onClick={() => setActiveTab('references')}
            className={`flex items-center gap-2 px-6 py-3 font-medium transition-colors border-b-2 ${
              activeTab === 'references'
                ? 'border-primary text-primary'
                : 'border-transparent text-secondary hover:text-primary'
            }`}
          >
            <BookOpen className="w-4 h-4" />
            References
            {referenceDocument && (
              <span className="ml-1 px-2 py-0.5 text-xs rounded-full bg-primary/20 text-primary">
                {referenceDocument.status === 'available' ? 'PDF' : 'Reupload'}
              </span>
            )}
          </button>
        </div>

        {/* Tab Content */}
        <div className="modal-body">
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

          {activeTab === 'molecule' && (
            <MoleculePropertiesContent
              moleculeName={pendingMoleculeName}
              onMoleculeNameChange={setPendingMoleculeName}
            />
          )}

          {activeTab === 'references' && (
            <div className="space-y-4">
              {referenceDocument && pendingReference.length === 0 && !referenceRemoved && (
                <div className="attachment-file-card">
                  <div className="attachment-file-thumbnail attachment-file-icon">
                    <FileText className="w-5 h-5" />
                  </div>
                  <div className="attachment-file-meta">
                    <div className="attachment-file-name">
                      {referenceDocument.title || referenceDocument.name}
                    </div>
                    <div className="attachment-file-detail">
                      {referenceDocument.name} - {sizeLabel(referenceDocument.sizeBytes)}
                      {referenceDocument.status === 'available' &&
                        ` - ${referenceDocument.pageCount || '?'} pages`}
                      {referenceDocument.status === 'uploading' && ' - uploading'}
                      {referenceDocument.status === 'missing' && ' - reupload required'}
                      {referenceDocument.status === 'error' && ' - error'}
                    </div>
                  </div>
                  <button
                    type="button"
                    className="attachment-file-remove"
                    onClick={() => setReferenceRemoved(true)}
                    aria-label={`Remove ${referenceDocument.name}`}
                  >
                    <X className="w-4 h-4" />
                  </button>
                </div>
              )}

              {(!referenceDocument ||
                referenceDocument.status === 'missing' ||
                referenceDocument.status === 'error' ||
                referenceRemoved ||
                pendingReference.length > 0) && (
                <AttachmentUpload
                  value={pendingReference}
                  onChange={(next) => {
                    const nextReference = next.slice(-1);
                    setPendingReference(nextReference);
                    setReferenceRemoved(false);
                    if (nextReference[0] && onReferenceDocumentSave) {
                      setPendingReferenceUploaded(true);
                      void onReferenceDocumentSave(nextReference[0]);
                    }
                  }}
                  accept="application/pdf,.pdf"
                  acceptedMimeTypes={['application/pdf']}
                  maxFiles={1}
                  maxSizeBytes={100 * 1024 * 1024}
                  label="PDF reference"
                  emptyLabel="Attach PDF"
                  invalidTypeMessage="Only PDF files can be attached."
                  disabled={referenceUploadDisabled}
                />
              )}

              {referenceUploadDisabled && (
                <div className="alert alert-info">
                  <div className="text-sm text-secondary">
                    PDF upload requires an active websocket connection.
                  </div>
                </div>
              )}

              {referenceDocument?.status === 'missing' && pendingReference.length === 0 && (
                <div className="alert alert-info">
                  <div className="text-sm text-secondary">
                    This PDF was used in a previous session and must be reuploaded before
                    consult_with_document can access it.
                  </div>
                </div>
              )}

              {referenceDocument?.status === 'error' && (
                <div className="alert alert-info">
                  <div className="text-sm text-secondary">
                    {referenceDocument.error || 'The PDF could not be loaded.'}
                  </div>
                </div>
              )}
            </div>
          )}
        </div>

        {/* Footer with Apply/Cancel */}
        <div className="modal-footer">
          <button onClick={onClose} className="btn btn-tertiary">
            Cancel
          </button>
          <button onClick={handleApplyAndClose} className="btn btn-primary flex-1">
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M5 13l4 4L19 7"
              />
            </svg>
            Apply Changes
          </button>
        </div>
      </div>
    </div>
  );
};
