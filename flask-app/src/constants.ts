export const MOLECULE_WIDTH = 250;
export const BOX_WIDTH = 10 + MOLECULE_WIDTH + 10;
export const BOX_GAP = 160;
export const BOX_HEIGHT = 100;

export const NODE_STYLES = {
  "normal": 'node-normal',
  "red": 'node-error',
  "yellow": 'node-computing',
} as const;

export const DEFAULT_CUSTOM_SYSTEM_PROMPT = 'Respond with concise and clear answers, as well as JSONs for each SMILES string you return as `{"smiles": "<SMILES STRING GOES HERE>"}`.';