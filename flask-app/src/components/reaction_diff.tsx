import React from 'react'
import { RDKitModule } from '@rdkit/rdkit'
import { RdkitjsReactionPayload } from '../types'
import { renderReactionPayload } from '../rdkit/reactionPayload'

type Props = {
  payload: RdkitjsReactionPayload
  rdkitModule: RDKitModule | null
  width?: number
  height?: number
}

export const ReactionDiff: React.FC<Props> = ({ payload, rdkitModule, width = 120, height = 90 }) => {
  if (!payload) return null
  if (!rdkitModule) return null
  let reactants: string[] = []
  let products: string[] = []
  try {
    const { reactantsSvg, productsSvg } = renderReactionPayload(rdkitModule as any, payload, { width, height })
    reactants = reactantsSvg
    products = productsSvg
  } catch (e) {
    // If rendering fails, show nothing
    return null
  }

  return (
    <div className="space-y-2">
      <div className="flex items-center gap-2 flex-wrap">
        {reactants.map((svg, i) => (
          <div key={`r-${i}`} className="bg-white/70 rounded" dangerouslySetInnerHTML={{ __html: svg }} />
        ))}
      </div>
      <div className="text-xs text-tertiary text-center">→</div>
      <div className="flex items-center gap-2 flex-wrap">
        {products.map((svg, i) => (
          <div key={`p-${i}`} className="bg-white/70 rounded" dangerouslySetInnerHTML={{ __html: svg }} />
        ))}
      </div>
    </div>
  )
}
