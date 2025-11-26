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

export const WS_SERVER = getConfig().WS_SERVER;
export const VERSION = getConfig().VERSION;
