import React, { useState } from 'react';
import { X, RotateCcw } from 'lucide-react';
import { OptimizationCustomization } from '../types';

interface CustomizationModalProps {
  isOpen: boolean;
  onClose: () => void;
  initialCustomization: OptimizationCustomization;
  onSave: (customization: OptimizationCustomization) => void;
}

const DEFAULT_CUSTOMIZATION: OptimizationCustomization = {
  molecularSimilarity: 0.7,
  diversityPenalty: 0.0,
  explorationRate: 0.5,
  enableConstraints: false,
};

export const CustomizationModal: React.FC<CustomizationModalProps> = ({
  isOpen,
  onClose,
  initialCustomization,
  onSave,
}) => {
  const [customization, setCustomization] = useState<OptimizationCustomization>(initialCustomization);

  if (!isOpen) return null;

  const handleSave = () => {
    onSave(customization);
  };

  const handleReset = () => {
    setCustomization(DEFAULT_CUSTOMIZATION);
  };

  return (
    <div className="modal-overlay">
      <div className="modal-content modal-content-lg">
        <div className="modal-header">
          <div>
            <h2 className="modal-title">Optimization Customization</h2>
            <p className="modal-subtitle">Configure advanced optimization parameters</p>
          </div>
          <button onClick={onClose} className="btn-icon">
            <X className="w-6 h-6" />
          </button>
        </div>

        <div className="modal-body space-y-6">
          {/* Molecular Similarity */}
          <div className="form-group">
            <label className="form-label">
              Molecular Similarity Threshold
              <span className="text-sm text-tertiary ml-2">
                (0.0 - 1.0, higher = more similar to parent)
              </span>
            </label>
            <div className="flex items-center gap-4">
              <input
                type="range"
                min="0"
                max="1"
                step="0.05"
                value={customization.molecularSimilarity}
                onChange={(e) => setCustomization({
                  ...customization,
                  molecularSimilarity: parseFloat(e.target.value)
                })}
                className="flex-1"
              />
              <input
                type="number"
                min="0"
                max="1"
                step="0.05"
                value={customization.molecularSimilarity}
                onChange={(e) => setCustomization({
                  ...customization,
                  molecularSimilarity: parseFloat(e.target.value)
                })}
                className="form-input w-20"
              />
            </div>
            <p className="text-xs text-tertiary mt-1">
              Controls how similar generated molecules should be to the lead molecule. Higher values promote conservative modifications.
            </p>
          </div>

          {/* Diversity Penalty */}
          <div className="form-group">
            <label className="form-label">
              Diversity Penalty
              <span className="text-sm text-tertiary ml-2">
                (0.0 - 1.0, higher = more penalty for similar molecules)
              </span>
            </label>
            <div className="flex items-center gap-4">
              <input
                type="range"
                min="0"
                max="1"
                step="0.05"
                value={customization.diversityPenalty}
                onChange={(e) => setCustomization({
                  ...customization,
                  diversityPenalty: parseFloat(e.target.value)
                })}
                className="flex-1"
              />
              <input
                type="number"
                min="0"
                max="1"
                step="0.05"
                value={customization.diversityPenalty}
                onChange={(e) => setCustomization({
                  ...customization,
                  diversityPenalty: parseFloat(e.target.value)
                })}
                className="form-input w-20"
              />
            </div>
            <p className="text-xs text-tertiary mt-1">
              Penalizes generation of molecules too similar to previously generated ones. Higher values encourage chemical space exploration.
            </p>
          </div>

          {/* Exploration Rate */}
          <div className="form-group">
            <label className="form-label">
              Exploration Rate
              <span className="text-sm text-tertiary ml-2">
                (0.0 - 1.0, higher = more exploration)
              </span>
            </label>
            <div className="flex items-center gap-4">
              <input
                type="range"
                min="0"
                max="1"
                step="0.05"
                value={customization.explorationRate}
                onChange={(e) => setCustomization({
                  ...customization,
                  explorationRate: parseFloat(e.target.value)
                })}
                className="flex-1"
              />
              <input
                type="number"
                min="0"
                max="1"
                step="0.05"
                value={customization.explorationRate}
                onChange={(e) => setCustomization({
                  ...customization,
                  explorationRate: parseFloat(e.target.value)
                })}
                className="form-input w-20"
              />
            </div>
            <p className="text-xs text-tertiary mt-1">
              Balance between exploiting good candidates (low) vs exploring new regions of chemical space (high).
            </p>
          </div>

          {/* Chemical Constraints */}
          <div className="form-group">
            <label className="form-label flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={customization.enableConstraints}
                onChange={(e) => setCustomization({
                  ...customization,
                  enableConstraints: e.target.checked
                })}
                className="form-checkbox"
              />
              Enable Chemical Constraints
            </label>
            <p className="text-xs text-tertiary mt-1">
              Apply drug-likeness rules and synthetic accessibility filters during optimization.
            </p>
          </div>
        </div>

        <div className="modal-footer">
          <button
            onClick={handleReset}
            className="btn btn-tertiary"
          >
            <RotateCcw className="w-4 h-4" />
            Reset to Defaults
          </button>
          <button
            onClick={handleSave}
            className="btn btn-primary flex-1"
          >
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
            </svg>
            Apply Customization
          </button>
        </div>
      </div>
    </div>
  );
};
