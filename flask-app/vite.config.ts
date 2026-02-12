import { defineConfig } from 'vite';
import tailwindcss from '@tailwindcss/vite'
import react from '@vitejs/plugin-react';
import { resolve } from 'path'
import { copyFileSync, mkdirSync, writeFileSync } from 'fs'

// Plugin to copy RDKit files during build
function copyRDKitFiles() {
  return {
    name: 'copy-rdkit-files',
    buildStart() {
      try {
        // Create dist directory in public folder
        mkdirSync(resolve(__dirname, 'public/rdkit'), { recursive: true });

        // Copy RDKit files
        copyFileSync(
          resolve(__dirname, 'node_modules/@rdkit/rdkit/dist/RDKit_minimal.js'),
          resolve(__dirname, 'public/rdkit/RDKit_minimal.js')
        );
        copyFileSync(
          resolve(__dirname, 'node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm'),
          resolve(__dirname, 'public/rdkit/RDKit_minimal.wasm')
        );
        console.log('✓ RDKit files copied to assets');
      } catch (error) {
        console.error('Failed to copy RDKit files:', error);
      }
    }
  }
}

export default defineConfig({
  plugins: [
    react(),
    tailwindcss(),
    copyRDKitFiles()
  ],
  define: {
    'window.APP_CONFIG.WS_SERVER': JSON.stringify(process.env.WS_SERVER || 'ws://localhost:8001/ws'),
    'window.APP_CONFIG.VERSION': JSON.stringify(process.env.SERVER_VERSION || '')
  },
  resolve: {
    dedupe: ['react', 'react-dom']
  },
  optimizeDeps: {
    // Force pre-bundling in dev AND specify for build
    include: ['react', 'react-dom', 'lc-conductor', 'scheduler'],
    force: true  // Force re-optimization
  },
  build: {
    rollupOptions: {
      // Completely disable CommonJS detection for scheduler
      commonjsOptions: {
        exclude: ['scheduler'],  // Don't let commonjs plugin touch it
      },
      // Explicitly tell Rollup: DO NOT externalize scheduler
      external: (id) => {
        if (id.includes('scheduler')) {
          console.log('🔍 Checking scheduler:', id, '-> BUNDLE IT!');
          return false;
        }
        return false;
      }
    }
  }
});
