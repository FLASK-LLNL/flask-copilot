import React from 'react';
import { RotateCcw } from 'lucide-react';
import { OptimizationCustomization, ConstraintOption } from '../types';

interface OptimizationCustomizationContentProps {
  customization: OptimizationCustomization;
  onCustomizationChange: (customization: OptimizationCustomization) => void;
  showResetButton?: boolean;
}

export const DEFAULT_CUSTOMIZATION: OptimizationCustomization = {
  enableConstraints: false,
  molecularSimilarity: 0.7,
  diversityPenalty: 0.0,
  explorationRate: 0.5,
  additionalConstraints: [],
};

export const CONSTRAINT_OPTIONS: ConstraintOption[] = [
  {
    value: 'drug-likeness',
    label: 'Drug-likeness (Lipinski\'s Rule)',
    description: 'Ensure molecules satisfy Lipinski\'s Rule of Five for oral bioavailability'
  },
  {
    value: 'synthesizability',
    label: 'Synthetic Accessibility',
    description: 'Prioritize molecules with high synthetic accessibility scores'
  },
  {
    value: 'lead-likeness',
    label: 'Lead-likeness',
    description: 'Apply lead-likeness criteria for early drug discovery'
  },
  {
    value: 'pan-assay-interference',
    label: 'PAINS Filter',
    description: 'Filter out Pan-Assay Interference Compounds (PAINS)'
  },
  {
    value: 'toxicity-rules',
    label: 'Toxicity Rules',
    description: 'Apply structural alerts for potential toxicity'
  },
  {
    value: 'reactive-groups',
    label: 'Reactive Group Filter',
    description: 'Avoid molecules with highly reactive functional groups'
  }
];

export const OptimizationCustomizationContent: React.FC<OptimizationCustomizationContentProps> = ({
  customization,
  onCustomizationChange,
  showResetButton = true,
}) => {
  const handleReset = () => {
    onCustomizationChange(DEFAULT_CUSTOMIZATION);
  };

  const toggleConstraint = (value: string) => {
    const constraints = customization.additionalConstraints || [];
    const newConstraints = constraints.includes(value)
      ? constraints.filter(c => c !== value)
      : [...constraints, value];

    onCustomizationChange({
      ...customization,
      additionalConstraints: newConstraints
    });
  };

  const isCustomizationEnabled = customization.enableConstraints ?? false;

  return (
    <div className="space-y-6">
      {/* Enable/Disable Customization */}
      <div className="form-group">
        <label className="form-label flex items-center gap-2 cursor-pointer">
          <input
            type="checkbox"
            checked={isCustomizationEnabled}
            onChange={(e) => onCustomizationChange({
              ...customization,
              enableConstraints: e.target.checked
            })}
            className="form-checkbox"
          />
          <span className="font-semibold">Enable Custom Optimization Strategy</span>
        </label>
        <p className="text-xs text-tertiary mt-1">
          When disabled, default optimization strategy will be used. Enable to customize parameters below.
        </p>
      </div>

      {/* Divider */}
      {isCustomizationEnabled && (
        <div className="border-t border-secondary pt-4 space-y-6">
          {showResetButton && (
            <div className="flex justify-end">
              <button
                onClick={handleReset}
                className="btn btn-tertiary btn-sm"
              >
                <RotateCcw className="w-4 h-4" />
                Reset to Defaults
              </button>
            </div>
          )}

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
                onChange={(e) => onCustomizationChange({
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
                onChange={(e) => onCustomizationChange({
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
                onChange={(e) => onCustomizationChange({
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
                onChange={(e) => onCustomizationChange({
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
                onChange={(e) => onCustomizationChange({
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
                onChange={(e) => onCustomizationChange({
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

          {/* Additional Constraints */}
          <div className="form-group">
            <label className="form-label mb-3">
              Additional Constraints
              <span className="text-sm text-tertiary ml-2">
                (select all that apply)
              </span>
            </label>
            <div className="space-y-3 max-h-64 overflow-y-auto custom-scrollbar pr-2">
              {CONSTRAINT_OPTIONS.map((option) => (
                <label
                  key={option.value}
                  className="flex items-start gap-3 p-3 rounded-lg border border-secondary hover:border-primary transition-colors cursor-pointer"
                >
                  <input
                    type="checkbox"
                    checked={(customization.additionalConstraints || []).includes(option.value)}
                    onChange={() => toggleConstraint(option.value)}
                    className="form-checkbox mt-1"
                  />
                  <div className="flex-1">
                    <div className="text-sm font-medium text-primary">
                      {option.label}
                    </div>
                    <div className="text-xs text-tertiary mt-1">
                      {option.description}
                    </div>
                  </div>
                </label>
              ))}
            </div>
            <p className="text-xs text-tertiary mt-2">
              Selected constraints will guide molecule generation to meet specific chemical criteria.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};
