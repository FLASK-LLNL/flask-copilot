export interface AppConfig {
  // WebSocket
  WS_SERVER: string;
  VERSION: string;
};

declare global {
  interface Window {
    APP_CONFIG?: Partial<AppConfig>;
    VERSION?: Partial<AppConfig>;
  }
}

const DEFAULT_CONFIG: AppConfig = {
  WS_SERVER: 'ws://localhost:8001/ws',
  VERSION: ''
};

let config: AppConfig | null = null;

export const getConfig = (): AppConfig => {
    if (!config) {
        config = {
        ...DEFAULT_CONFIG,
        ...(window.APP_CONFIG || {})
        };
    }
    return config;
};

function wsToHttp(wsUrl: string): string {
  // Replace ws:// with http:// and wss:// with https://
  let httpUrl = wsUrl.replace(/^wss:\/\//, 'https://').replace(/^ws:\/\//, 'http://');

  // Find the last slash and keep everything up to and including it
  const lastSlashIndex = httpUrl.lastIndexOf('/');

  // Ensure we're not cutting the protocol's slashes
  if (lastSlashIndex > httpUrl.indexOf('://') + 2) {
    httpUrl = httpUrl.substring(0, lastSlashIndex);
  }

  return httpUrl;
}

export const WS_SERVER = getConfig().WS_SERVER;
export const VERSION = getConfig().VERSION;
export const HTTP_SERVER = wsToHttp(WS_SERVER);
