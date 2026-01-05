import React from 'react';
import { X, Loader2, FlaskConical, BookOpen, Check, ChevronDown, ChevronRight, AlertCircle } from 'lucide-react';
import { ReactionAlternative } from '../types';
import { RDKitModule } from '@rdkit/rdkit';

interface ReactionAlternativesSidebarProps {
  isOpen: boolean;
  onClose: () => void;
  productMolecule: string;
  productSmiles: string;
  alternatives: ReactionAlternative[];
  onSelectAlternative: (alt: ReactionAlternative) => void;
  onComputeTemplates: () => void;
  onComputeFlaskAI: () => void;
  isComputingTemplates: boolean;
  isComputing: boolean;
  templatesSearched: boolean;
  rdkitModule: RDKitModule | null;
}

// Mini molecule preview component
const MiniMoleculeSVG: React.FC<{ smiles: string; rdkitModule: RDKitModule | null }> = ({ smiles, rdkitModule }) => {
  const [svg, setSvg] = React.useState<string>('');
  const renderedRef = React.useRef<string>('');

  React.useEffect(() => {
    if (!rdkitModule || renderedRef.current === smiles || !smiles) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (mol && mol.is_valid()) {
        // Configure drawing options for dark mode
        const drawOpts = JSON.stringify({
          width: 45,
          height: 45,
          backgroundColour: [1, 1, 1, 0.0], // Transparent [R, G, B, A]
          // Lighter colors for dark backgrounds
          bondLineWidth: 1,
        });

        const svgStr = mol.get_svg_with_highlights(drawOpts);
        //const svgStr = mol.get_svg(45, 45);
        setSvg(svgStr);
        renderedRef.current = smiles;
        mol.delete();
      }
    } catch (e) {
      console.error('Error rendering mini molecule:', e);
    }
  }, [smiles, rdkitModule]);

  if (!svg) {
    return (
      <div className="w-9 h-9 rounded bg-surface flex items-center justify-center flex-shrink-0">
        <div className="w-3 h-3 rounded-full bg-primary/30"></div>
      </div>
    );
  }

  return (
    <div
      className="w-9 h-9 rounded bg-white/70 flex items-center justify-center flex-shrink-0"
      dangerouslySetInnerHTML={{ __html: svg }}
    />
  );
};

// Multi-step pathway preview with support for multiple reactants per step
const PathwayPreview: React.FC<{
  alternative: ReactionAlternative;
  rdkitModule: RDKitModule | null;
  isActive: boolean;
  isDisabled: boolean;
}> = ({ alternative, rdkitModule, isActive, isDisabled }) => {
  if (!alternative.pathway || alternative.pathway.length === 0) {
    return (
      <div className="flex items-center gap-2 py-2 text-xs text-tertiary italic">
        No preview available
      </div>
    );
  }

  const totalSteps = alternative.pathway.length - 1;

  return (
    <div className={`flex flex-col gap-2 py-2 ${isDisabled ? 'opacity-50' : ''}`}>
      {/* Step count indicator */}
      <div className="text-xs text-tertiary">
        {totalSteps} step{totalSteps !== 1 ? 's' : ''}
      </div>

      {/* Pathway visualization */}
      <div className="flex items-center gap-1 overflow-x-auto custom-scrollbar pb-2">
        {alternative.pathway.map((step, stepIdx) => (
          <React.Fragment key={stepIdx}>
            {/* Group of molecules at this step */}
            <div className="flex flex-col gap-1">
              <div className={`flex items-center gap-1 ${step.smiles.length > 1 ? 'flex-wrap max-w-[120px]' : ''}`}>
                {step.smiles.map((smiles, molIdx) => (
                  <React.Fragment key={molIdx}>
                    <div className="flex flex-col items-center gap-0.5 flex-shrink-0">
                      <MiniMoleculeSVG smiles={smiles} rdkitModule={rdkitModule} />
                      <div className="text-[8px] text-tertiary max-w-[45px] truncate text-center" title={step.label[molIdx]}>
                        {step.label[molIdx]}
                      </div>
                    </div>
                    {molIdx < step.smiles.length - 1 && (
                      <div className="text-[10px] text-muted">+</div>
                    )}
                  </React.Fragment>
                ))}
              </div>
            </div>

            {/* Arrow to next step */}
            {stepIdx < alternative.pathway.length - 1 && (
              <ChevronRight className="w-3 h-3 text-muted flex-shrink-0 mx-0.5" />
            )}
          </React.Fragment>
        ))}
      </div>
    </div>
  );
};

// Collapsible section component
const CollapsibleSection: React.FC<{
  title: string;
  count: number;
  dotColor: string;
  subtitle?: string;
  defaultExpanded?: boolean;
  children: React.ReactNode;
}> = ({ title, count, dotColor, subtitle, defaultExpanded = false, children }) => {
  const [isExpanded, setIsExpanded] = React.useState(defaultExpanded);

  return (
    <div>
      <button
        onClick={() => setIsExpanded(!isExpanded)}
        className="w-full flex items-center justify-between p-2 hover:bg-surface-hover transition-colors rounded"
      >
        <div className="flex items-center gap-2">
          <div className={`w-2 h-2 rounded-full ${dotColor}`}></div>
          <span className="section-label">{title} ({count})</span>
        </div>
        {isExpanded ? (
          <ChevronDown className="w-4 h-4 text-muted" />
        ) : (
          <ChevronRight className="w-4 h-4 text-muted" />
        )}
      </button>
      {subtitle && (
        <div className="text-xs text-tertiary ml-6 mb-2 italic">{subtitle}</div>
      )}
      {isExpanded && (
        <div className="space-y-2 ml-2 mt-2">
          {children}
        </div>
      )}
    </div>
  );
};

export const ReactionAlternativesSidebar: React.FC<ReactionAlternativesSidebarProps> = ({
  isOpen,
  onClose,
  productMolecule,
  productSmiles,
  alternatives,
  onSelectAlternative,
  onComputeTemplates,
  onComputeFlaskAI,
  isComputingTemplates,
  isComputing,
  templatesSearched,
  rdkitModule
}) => {
  const [hoveredDisabled, setHoveredDisabled] = React.useState<string | null>(null);

  if (!isOpen) return null;

  const exactMatches = alternatives.filter(a => a.type === 'exact');
  const templateMatches = alternatives.filter(a => a.type === 'template');

  const hasExactMatches = exactMatches.length > 0;
  const hasTemplateMatches = templateMatches.length > 0;

  // Show template button only if templates haven't been searched yet
  const showTemplateButton = !templatesSearched;

  const AlternativeCard: React.FC<{ alt: ReactionAlternative }> = ({ alt }) => {
    const isActive = alt.status === 'active';
    const isDisabled = alt.disabled || false;
    const canClick = !isActive && !alt.disabled && alt.status !== 'computing';

    return (
      <div className="relative">
        <button
          onClick={() => canClick && onSelectAlternative(alt)}
          disabled={!canClick}
          onMouseEnter={() => isDisabled && alt.disabledReason && setHoveredDisabled(alt.id)}
          onMouseLeave={() => setHoveredDisabled(null)}
          className={`w-full glass-panel transition-all text-left relative ${
            isActive ? 'border-2 border-primary bg-primary/10' : ''
          } ${
            isDisabled ? 'opacity-60 cursor-not-allowed' : 'hover:bg-surface-hover'
          } ${
            alt.status === 'computing' || isActive ? 'cursor-default' : ''
          }`}
        >
          <div className="flex items-start justify-between mb-1">
            <div className="flex items-center gap-2 flex-1 min-w-0">
              {isDisabled ? (
                <AlertCircle className="w-3.5 h-3.5 text-warning flex-shrink-0" />
              ) : (
                <BookOpen className="w-3.5 h-3.5 text-muted flex-shrink-0" />
              )}
              <span className={`text-sm font-medium truncate ${isDisabled ? 'text-muted' : 'text-primary'}`}>
                {alt.name}
              </span>
            </div>
            {isActive && !isDisabled && (
              <Check className="w-4 h-4 text-success flex-shrink-0" />
            )}
            {alt.status === 'computing' && (
              <Loader2 className="w-4 h-4 animate-spin text-warning flex-shrink-0" />
            )}
          </div>
          <PathwayPreview
            alternative={alt}
            rdkitModule={rdkitModule}
            isActive={isActive}
            isDisabled={isDisabled}
          />
        </button>

        {/* Disabled reason tooltip */}
        {hoveredDisabled === alt.id && alt.disabledReason && (
          <div className="absolute left-0 right-0 top-full mt-1 z-50 pointer-events-none">
            <div className="bg-warning/20 border border-warning/50 rounded px-3 py-2 text-xs text-warning backdrop-blur-sm">
              <div className="flex items-start gap-2">
                <AlertCircle className="w-3 h-3 flex-shrink-0 mt-0.5" />
                <span>{alt.disabledReason}</span>
              </div>
            </div>
          </div>
        )}
      </div>
    );
  };

  return (
    <div
      className="absolute top-0 right-0 bottom-0 z-50 animate-slideInSidebar"
      style={{ width: '420px' }}
      onClick={(e) => e.stopPropagation()}
    >
      <div className="h-full flex flex-col sidebar-inner backdrop-blur-lg border-l-2 border-primary shadow-2xl">
        <div className="flex items-center justify-between p-4 border-b border-secondary">
          <div className="flex-1 min-w-0">
            <h3 className="text-lg font-semibold text-primary truncate">
              Synthesis Pathways
            </h3>
            <p className="text-sm text-secondary truncate" title={productMolecule}>
              for {productMolecule}
            </p>
          </div>
          <button onClick={onClose} className="btn-icon ml-2">
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto custom-scrollbar p-4 space-y-4">

          {/* AI Computation Button - Top Priority */}
          <div className="pb-2 border-b border-secondary">
            <button
              onClick={() => {
                onComputeFlaskAI();
                onClose();
              }}
              className="btn btn-primary w-full"
              disabled={isComputing}
            >
              <FlaskConical className="w-4 h-4" />
              Compute New Path with AI
            </button>
            <div className="text-xs text-tertiary mt-2 italic">
              Discovers novel pathways using machine learning
            </div>
          </div>

          {/* No Exact Matches Message */}
          {!hasExactMatches && templatesSearched && (
            <div className="glass-panel bg-warning/10 border border-warning/30 p-3">
              <div className="flex items-start gap-2">
                <AlertCircle className="w-4 h-4 text-warning flex-shrink-0 mt-0.5" />
                <div className="text-sm text-warning">
                  No exact matches found in database
                </div>
              </div>
            </div>
          )}

          {/* Exact Matches Section */}
          {hasExactMatches && (
            <CollapsibleSection
              title="EXACT MATCHES"
              count={exactMatches.length}
              dotColor="bg-success"
              subtitle="Known literature reactions"
              defaultExpanded={true}
            >
              {exactMatches.map(alt => (
                <AlternativeCard key={alt.id} alt={alt} />
              ))}
            </CollapsibleSection>
          )}

          {/* Template-Based Matches Section */}
          {(hasTemplateMatches || isComputingTemplates) && (
            <CollapsibleSection
              title="TEMPLATE-BASED MATCHES"
              count={templateMatches.length}
              dotColor="bg-info"
              subtitle="Adapted from similar reactions"
              defaultExpanded={!hasExactMatches}
            >
              {isComputingTemplates && templateMatches.length === 0 && (
                <div className="text-sm text-tertiary italic p-4 text-center glass-panel flex items-center justify-center gap-2">
                  <Loader2 className="w-4 h-4 animate-spin" />
                  Searching templates...
                </div>
              )}
              {templateMatches.map(alt => (
                <AlternativeCard key={alt.id} alt={alt} />
              ))}
            </CollapsibleSection>
          )}

          {/* No alternatives yet */}
          {!hasExactMatches && !hasTemplateMatches && !isComputingTemplates && templatesSearched && (
            <div className="text-sm text-tertiary italic p-4 text-center glass-panel">
              No alternative pathways found
            </div>
          )}

          {/* Find Template-Based Matches Button */}
          {showTemplateButton && (
            <div className="pt-2">
              <button
                onClick={onComputeTemplates}
                disabled={isComputingTemplates || isComputing}
                className="btn btn-secondary w-full"
              >
                {isComputingTemplates ? (
                  <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    Searching Templates...
                  </>
                ) : (
                  <>
                    <BookOpen className="w-4 h-4" />
                    {hasExactMatches
                      ? 'Find Template-Based Matches'
                      : 'Find Template-Based Pathways'}
                  </>
                )}
              </button>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
