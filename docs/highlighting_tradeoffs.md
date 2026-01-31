# Atom-Diff Highlighting On Graph: SMILES vs MolBlock

We moved changed-atom highlighting from the "Synthesis Pathways" sidebar into the main graph view.
Highlighting is triggered by hovering a reaction badge and draws highlights directly on the connected molecule cards.

This note captures the key rendering decision and its trade-offs.

## Decision

Render hover-highlights from **SMILES**, not from backend-provided **MolBlock**.

In practice:
- Normal molecule cards are rendered from `smiles`.
- Highlighted hover state is also rendered from `smiles`.

## Why Not MolBlock?

MolBlock-based rendering aligns atom indices between backend diff computation and frontend highlighting.
However, MolBlock depictions frequently differ in 2D orientation/rotation from SMILES depictions.

Observed UI impact:
- On hover, molecule pictures "jump" (rotation/layout changes), which is visually jarring and harms usability.

## Trade-offs

### SMILES-based highlights (current)

Pros:
- Visual continuity: highlighted/non-highlighted depictions look consistent (minimal "jump" on hover).
- Matches the existing rendering path used throughout the graph.

Cons:
- **Atom index mismatch risk**: backend computes `highlight_atom_idxs` against its RDKit molecule ordering, but the frontend recreates molecules from SMILES in rdkit.js.
- When atom ordering differs, highlights may land on the wrong atoms.

### MolBlock-based highlights

Pros:
- Correctness: atom indices from the backend match the rendered molecule exactly.
- Less sensitive to SMILES canonicalization and atom reordering.

Cons:
- Visual discontinuity: depictions can rotate/relayout when switching between SMILES and MolBlock rendering.
- This was noticeable enough to reject for the main graph hover interaction.

## Recommended Path If Correctness Becomes Priority

If we need both visual continuity and correctness, we should change the payload contract rather than toggling render inputs:

Option A (preferred):
- Backend sends a stable mapping identifier (atom-map numbers) and the frontend computes rdkit.js atom indices from the rendered molecule.

Option B:
- Backend computes highlights using the same canonical atom ordering as the frontend SMILES (hard to guarantee across toolchains).

## Current Status

- Sidebar diff panel removed.
- Graph reaction hover highlights enabled.
- Highlights rendered from SMILES for consistency, accepting the known index-mismatch risk.
