import { useState, useEffect } from "react";
import { MOLECULE_WIDTH } from "../constants";
import { MoleculeSVGProps } from "../types";
import { RDKitModule } from '@rdkit/rdkit';

export const MoleculeSVG: React.FC<MoleculeSVGProps> = ({ smiles, height = 80, rdkitModule = null }) => {
  const [svg, setSvg] = useState<string | null>(null);
  const [error, setError] = useState<boolean>(false);

  useEffect(() => {
    if (rdkitModule && smiles) {
      try {
        const mol = rdkitModule.get_mol(smiles);
        if (mol && mol.is_valid()) {
          // Configure drawing options for dark mode
          const drawOpts = JSON.stringify({
            width: MOLECULE_WIDTH,
            height: height,
            backgroundColour: [1, 1, 1, 0.0], // Transparent [R, G, B, A]
            // Lighter colors for dark backgrounds
            bondLineWidth: 1,
          });

          const molSvg = mol.get_svg_with_highlights(drawOpts);

          // Use this to change bond colors to "dark mode" (maybe not a good idea)
          const modifiedSvg = molSvg
              //.replace(/fill:"#FFFFFF"/g, "fill='transparent'")
              //.replace(/stroke:#000000/g, "stroke:#E5E7EB")
              //.replace(/fill:#000000/g, "fill:#E5E7EB");
          setSvg(modifiedSvg);
          mol.delete();
          return;
        }
      } catch (e) {
        console.error('RDKit rendering error:', e);
        setError(true);
      }
    }
    setError(true);
  }, [smiles, rdkitModule, height]);


  // Fallback to dummy visualization
  if (error || !rdkitModule) {
    const colors = ['#3B82F6', '#8B5CF6', '#EC4899', '#10B981', '#F59E0B'];
    const hash = smiles.split('').reduce((acc, char) => acc + char.charCodeAt(0), 0);
    const color = colors[hash % colors.length];
    const nodes = 5 + (hash % 3);
    
    const points: [number, number][] = [];
    for (let i = 0; i < nodes; i++) {
      const angle = (i * 2 * Math.PI) / nodes;
      const r = 25;
      points.push([40 + r * Math.cos(angle), 40 + r * Math.sin(angle)]);
    }
    
    return (
      <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', width: '100%', height: height }}>
        <svg height={height} viewBox="0 0 80 80" style={{ display: 'block', margin: '0 auto' }}>
          <defs>
            <linearGradient id={`grad-${hash}`} x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" style={{ stopColor: color, stopOpacity: 0.8 }} />
              <stop offset="100%" style={{ stopColor: color, stopOpacity: 0.4 }} />
            </linearGradient>
          </defs>
          {points.map((point, i) => {
            const nextPoint = points[(i + 1) % points.length];
            return (
              <line
                key={i}
                x1={point[0]}
                y1={point[1]}
                x2={nextPoint[0]}
                y2={nextPoint[1]}
                stroke={color}
                strokeWidth="2"
              />
            );
          })}
          {points.map((point, i) => (
            <circle
              key={i}
              cx={point[0]}
              cy={point[1]}
              r="6"
              fill={`url(#grad-${hash})`}
              stroke={color}
              strokeWidth="2"
            />
          ))}
        </svg>
      </div>
    );
  }

  return <div dangerouslySetInnerHTML={{ __html: svg || '' }} />;
};

export const loadRDKit = (): RDKitModule | null => {
    const [rdkitModule, setRdkitModule] = useState<RDKitModule | null>(null);

    // Load RDKit.js on mount
    useEffect(() => {
        const script = document.createElement('script');
        script.src = '/rdkit/RDKit_minimal.js';
        script.async = true;
        
        script.onload = () => {
        window.initRDKitModule().then((RDKit: RDKitModule) => {
            console.log('RDKit loaded successfully!');
            setRdkitModule(RDKit);
        }).catch((error: any) => {
            console.error('RDKit initialization failed:', error);
        });
        };
        
        script.onerror = () => {
        console.error('Failed to load RDKit script');
        };
        
        document.body.appendChild(script);
        
        return () => {
        if (document.body.contains(script)) {
            document.body.removeChild(script);
        }
        };
    }, []);

    return rdkitModule;
};
