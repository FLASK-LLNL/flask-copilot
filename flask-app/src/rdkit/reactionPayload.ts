export type RGB = [number, number, number]

export type RdkitjsMolPayload = {
  molblock: string
  smiles: string
  highlight_atom_idxs: number[]
  highlight_atom_mapnums: number[]
}

export type RdkitjsReactionPayload = {
  version: 'rdkitjs-reaction-payload/v1'
  reactants: RdkitjsMolPayload[]
  products: RdkitjsMolPayload[]
  main_product_index: number
  highlight_rgb: RGB
  highlight_alpha: number
  reactant_mcs_smarts?: (string | null)[]
}

export type RdkitModule = {
  get_mol: (input: string) => RdkitMol
}

export type RdkitMol = {
  get_svg_with_highlights: (details: string) => string
  delete: () => void
}

function clamp01(x: number): number {
  if (x < 0) return 0
  if (x > 1) return 1
  return x
}

function blendToWhite(rgb: RGB, alpha: number): RGB {
  const a = clamp01(alpha)
  const [r, g, b] = rgb
  const rr = Math.round(a * r + (1 - a) * 255)
  const gg = Math.round(a * g + (1 - a) * 255)
  const bb = Math.round(a * b + (1 - a) * 255)
  return [rr, gg, bb]
}

export function renderMolWithHighlights(
  RDKit: RdkitModule,
  molblock: string,
  highlightAtomIdxs: number[],
  opts?: { highlightRgb?: RGB; highlightAlpha?: number; width?: number; height?: number },
): string {
  const width = opts?.width ?? 120
  const height = opts?.height ?? 90
  const rgb = opts?.highlightRgb ?? [220, 0, 0]
  const alpha = opts?.highlightAlpha ?? 0.45
  const softRgb = blendToWhite(rgb, alpha)
  const norm = (c: number) => Math.max(0, Math.min(1, c / 255))
  const softRgbNorm: [number, number, number] = [norm(softRgb[0]), norm(softRgb[1]), norm(softRgb[2])]

  const mol = RDKit.get_mol(molblock)
  try {
    const baseOpts = {
      width,
      height,
      backgroundColour: [1, 1, 1, 0.0],
      bondLineWidth: 1,
    }
    const base = mol.get_svg_with_highlights(JSON.stringify(baseOpts))

    // RDKit_minimal has had multiple highlight JSON schemas over time.
    // Prefer the variant that actually changes output from the base rendering.
    const candidates: any[] = [
      { ...baseOpts, atoms: highlightAtomIdxs },
      { ...baseOpts, atoms: highlightAtomIdxs, highlightColour: softRgbNorm },
      { ...baseOpts, highlightAtomIds: highlightAtomIdxs },
      { ...baseOpts, highlightAtomIds: highlightAtomIdxs, highlightColour: softRgbNorm },
      { ...baseOpts, highlightAtoms: highlightAtomIdxs },
    ]

    let lastErr: unknown = null
    for (const details of candidates) {
      try {
        const svg = mol.get_svg_with_highlights(JSON.stringify(details))
        if (svg && svg !== base) return svg
      } catch (e) {
        lastErr = e
      }
    }

    // Fallback to base render if none produced highlights.
    if (lastErr) {
      // Some builds accept keys but ignore them; keep base rather than throwing.
      return base
    }
    return base
  } finally {
    mol.delete()
  }
}

// Note: reaction-level rendering moved into the main graph hover highlighting.
