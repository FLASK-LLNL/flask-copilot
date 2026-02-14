import React from 'react';

export const MOLECULE_NAME_OPTIONS = [
  { value: 'brand', label: 'Brand/Common Name' },
  { value: 'iupac', label: 'IUPAC Name' },
  { value: 'formula', label: 'Chemical Formula' },
  { value: 'smiles', label: 'SMILES' }
];

interface MoleculePropertiesContentProps {
  moleculeName: string;
  onMoleculeNameChange: (name: string) => void;
}

export const MoleculePropertiesContent: React.FC<MoleculePropertiesContentProps> = ({
  moleculeName,
  onMoleculeNameChange,
}) => {
  return (
    <div className="space-y-4">
      <div>
        <h3 className="heading-3">Display Preferences</h3>
        <p className="helper-text">
          Configure how molecular information is displayed throughout the application
        </p>
      </div>

      {/* Molecule Name Preference */}
      <div className="form-group">
        <label className="form-label">
          Preferred Molecule Name Format
        </label>
        <select
          value={moleculeName || 'brand'}
          onChange={(e) => onMoleculeNameChange(e.target.value)}
          className="form-select"
        >
          {MOLECULE_NAME_OPTIONS.map(option => (
            <option key={option.value} value={option.value}>
              {option.label}
            </option>
          ))}
        </select>
        <p className="helper-text">
          Choose how molecule names are displayed throughout the application
        </p>
      </div>

      {/* Format Examples */}
      <div className="glass-panel">
        <h4 className="text-sm font-semibold mb-3 text-primary">Format Examples</h4>
        <div className="space-y-2 text-sm">
          <div className="flex items-start gap-3">
            <span className="text-mono emphasized-text min-w-[120px]">Brand Name:</span>
            <span className="text-secondary">Aspirin, Tylenol, Advil</span>
          </div>
          <div className="flex items-start gap-3">
            <span className="text-mono emphasized-text min-w-[120px]">IUPAC Name:</span>
            <span className="text-secondary">2-acetoxybenzoic acid</span>
          </div>
          <div className="flex items-start gap-3">
            <span className="text-mono emphasized-text min-w-[120px]">Formula:</span>
            <span className="text-secondary">C₉H₈O₄</span>
          </div>
          <div className="flex items-start gap-3">
            <span className="text-mono emphasized-text min-w-[120px]">SMILES:</span>
            <span className="text-secondary">CC(=O)Oc1ccccc1C(=O)O</span>
          </div>
        </div>
      </div>
    </div>
  );
};
