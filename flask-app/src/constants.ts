export const MOLECULE_WIDTH = 250;
export const BOX_WIDTH = 10 + MOLECULE_WIDTH + 10;
export const BOX_GAP = 160;
export const BOX_HEIGHT = 100;

export const NODE_STYLES = {
  "normal": 'from-purple-50/80 to-pink-300/80 border-purple-400/50 hover:border-purple-300',
  "red": 'from-amber-500/40 to-red-500/100 border-red-400 ring-4 ring-red-400/50',
  "yellow": 'from-amber-500/40 to-yellow-500/40 border-amber-400 ring-4 ring-amber-400/50 animate-pulse',
} as const;
