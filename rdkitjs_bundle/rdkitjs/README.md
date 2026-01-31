# rdkit.js Payload Notes

This folder contains reference code + payload shapes for rendering reaction atom-change
highlights in a TypeScript frontend using `rdkit.js`.

Best-practice approach
- Compute atom-change highlights on the server (Python/RDKit).
- Send a stable, explicit serialization for each molecule (MolBlock) plus the atom
  indices to highlight.
- Render client-side with `rdkit.js` by loading the MolBlock and calling
  `get_svg_with_highlights`.

Why MolBlock + atom indices
- A Python RDKit `Mol` cannot be sent directly to `rdkit.js`.
- MolBlock preserves atom ordering and is directly loadable by rdkit.js.
- Atom indices are the most direct input for rdkit.js highlighting.

Also included
- The payload includes per-atom map numbers (1..N) for debugging and future-proofing.

Files
- `rdkitjs/reactionPayload.ts`: TypeScript types + a renderer helper.
- `services/rdkitjs_payload.py`: Python builder for the payload.
