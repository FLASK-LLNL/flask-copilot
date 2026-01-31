export type RGB = [number, number, number]

export type RdkitjsMolPayload = {
  molblock: string
  smiles: string
  highlight_atom_idxs: number[]
  highlight_atom_mapnums: number[]
}

export type RdkitjsReactionPayload = {
  version: "rdkitjs-reaction-payload/v1"
  reactants: RdkitjsMolPayload[]
  products: RdkitjsMolPayload[]
  main_product_index: number
  highlight_rgb: RGB
  highlight_alpha: number
}

// Minimal rdkit.js surface needed for rendering.
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
  const width = opts?.width ?? 320
  const height = opts?.height ?? 240
  const rgb = opts?.highlightRgb ?? [220, 0, 0]
  const alpha = opts?.highlightAlpha ?? 0.45
  const softRgb = blendToWhite(rgb, alpha)

  const mol = RDKit.get_mol(molblock)
  try {
    const atomHighlights: Record<string, number[]> = {}
    for (const idx of highlightAtomIdxs) atomHighlights[String(idx)] = softRgb

    const details = {
      width,
      height,
      // rdkit.js expects "atoms" map from idx -> [r,g,b]
      // and "highlightAtomIds" array.
      atoms: atomHighlights,
      highlightAtomIds: highlightAtomIdxs,
    }

    return mol.get_svg_with_highlights(JSON.stringify(details))
  } finally {
    mol.delete()
  }
}

export function renderReactionPayload(
  RDKit: RdkitModule,
  payload: RdkitjsReactionPayload,
  opts?: { width?: number; height?: number },
): { reactantsSvg: string[]; productsSvg: string[] } {
  const width = opts?.width ?? 320
  const height = opts?.height ?? 240

  const reactantsSvg = payload.reactants.map((m) =>
    renderMolWithHighlights(RDKit, m.molblock, m.highlight_atom_idxs, {
      highlightRgb: payload.highlight_rgb,
      highlightAlpha: payload.highlight_alpha,
      width,
      height,
    }),
  )

  const productsSvg = payload.products.map((m) =>
    renderMolWithHighlights(RDKit, m.molblock, m.highlight_atom_idxs, {
      highlightRgb: payload.highlight_rgb,
      highlightAlpha: payload.highlight_alpha,
      width,
      height,
    }),
  )

  return { reactantsSvg, productsSvg }
}
