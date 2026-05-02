import React from 'react';
import { RsaSettingsPanel } from 'lc-conductor';

interface RetrosynthesisCustomizationContentProps {
  useAiBased: boolean;
  onUseAiBasedChange: (value: boolean) => void;
  useRsa: boolean;
  onUseRsaChange: (value: boolean) => void;
  rsaMode: 'standalone' | 'rag';
  onRsaModeChange: (mode: 'standalone' | 'rag') => void;
  rsaN: number;
  onRsaNChange: (value: number) => void;
  rsaK: number;
  onRsaKChange: (value: number) => void;
  rsaT: number;
  onRsaTChange: (value: number) => void;
}

export const RetrosynthesisCustomizationContent: React.FC<
  RetrosynthesisCustomizationContentProps
> = ({
  useAiBased,
  onUseAiBasedChange,
  useRsa,
  onUseRsaChange,
  rsaMode,
  onRsaModeChange,
  rsaN,
  onRsaNChange,
  rsaK,
  onRsaKChange,
  rsaT,
  onRsaTChange,
}) => {
  return (
    <div className="space-y-6">
      {/* Use AI-based Approach */}
      <div className="glass-panel">
        <div className="flex items-start gap-4">
          <input
            type="checkbox"
            checked={useAiBased}
            onChange={(e) => onUseAiBasedChange(e.target.checked)}
            className="form-checkbox mt-1"
            id="use-ai-based"
          />
          <div className="flex-1">
            <label htmlFor="use-ai-based" className="form-label-block cursor-pointer">
              Use AI-based Approach
            </label>
            <p className="text-sm text-tertiary mt-1">
              Use AI-based retrosynthesis instead of template-based methods. AI models can suggest
              novel reaction pathways beyond known templates.
            </p>
          </div>
        </div>
      </div>

      {/* Enable RSA Mode - Available independently */}
      <RsaSettingsPanel
        useRsa={useRsa}
        onUseRsaChange={onUseRsaChange}
        rsaMode={rsaMode}
        onRsaModeChange={onRsaModeChange}
        rsaN={rsaN}
        onRsaNChange={onRsaNChange}
        rsaK={rsaK}
        onRsaKChange={onRsaKChange}
        rsaT={rsaT}
        onRsaTChange={onRsaTChange}
      />

      {/* Information Panel */}
      <div className="glass-panel bg-primary/5 border border-primary/20">
        <div className="text-sm text-secondary">
          <p className="font-semibold mb-2">About Retrosynthesis Approaches:</p>
          <ul className="list-disc list-inside space-y-1 text-xs text-tertiary">
            <li>
              <strong>AI-based:</strong> Uses neural networks to suggest novel reactions and
              pathways not found in databases
            </li>
            <li>
              <strong>Template-based (disabled):</strong> Uses known reaction patterns from chemical
              databases
            </li>
            <li>RSA can be used with either approach for enhanced optimization</li>
            <li>Template-based is more conservative but limited to known chemistry</li>
          </ul>
        </div>
      </div>
    </div>
  );
};
