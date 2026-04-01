import { Client } from '@modelcontextprotocol/sdk/client/index.js';
import { StreamableHTTPClientTransport } from '@modelcontextprotocol/sdk/client/streamableHttp.js';

type LocalMcpClient = {
  listTools: () => Promise<{ tools?: Array<Record<string, any>> }>;
  callTool: (args: {
    name: string;
    arguments?: Record<string, any>;
  }) => Promise<Record<string, any>>;
  close?: () => Promise<void> | void;
};

type LocalMcpTransport = {
  close?: () => Promise<void> | void;
};

type LocalMcpClientEntry = {
  client: LocalMcpClient;
  transport: LocalMcpTransport;
};

export type LocalMcpToolDefinition = {
  name: string;
  description?: string;
  inputSchema?: Record<string, unknown>;
};

const clientCache = new Map<string, Promise<LocalMcpClientEntry>>();

export const normalizeMcpUrl = (url: string): string => {
  const trimmed = url.trim();
  if (!trimmed) {
    return trimmed;
  }
  return trimmed.endsWith('/mcp') ? trimmed : `${trimmed.replace(/\/+$/, '')}/mcp`;
};

const closeClientEntry = async (entry: LocalMcpClientEntry | undefined): Promise<void> => {
  if (!entry) {
    return;
  }

  await entry.client.close?.();
  await entry.transport.close?.();
};

const createClientEntry = async (serverUrl: string): Promise<LocalMcpClientEntry> => {
  const client = new Client(
    {
      name: 'flask-copilot-browser-proxy',
      version: '0.1.0',
    },
    {
      capabilities: {},
    }
  ) as LocalMcpClient;

  const transport = new StreamableHTTPClientTransport(new URL(serverUrl)) as LocalMcpTransport;

  await (client as any).connect(transport);
  return { client, transport };
};

const getClientEntry = async (serverUrl: string): Promise<LocalMcpClientEntry> => {
  const normalizedUrl = normalizeMcpUrl(serverUrl);
  const cached = clientCache.get(normalizedUrl);
  if (cached) {
    return cached;
  }

  const created = createClientEntry(normalizedUrl);
  clientCache.set(normalizedUrl, created);

  try {
    return await created;
  } catch (error) {
    clientCache.delete(normalizedUrl);
    throw error;
  }
};

const resetClientEntry = async (serverUrl: string): Promise<void> => {
  const normalizedUrl = normalizeMcpUrl(serverUrl);
  const cached = clientCache.get(normalizedUrl);
  clientCache.delete(normalizedUrl);
  if (!cached) {
    return;
  }

  try {
    await closeClientEntry(await cached);
  } catch {
    // Best-effort cleanup only.
  }
};

const withClientRetry = async <T>(
  serverUrl: string,
  callback: (entry: LocalMcpClientEntry) => Promise<T>
): Promise<T> => {
  const normalizedUrl = normalizeMcpUrl(serverUrl);

  for (let attempt = 0; attempt < 2; attempt += 1) {
    const entry = await getClientEntry(normalizedUrl);
    try {
      return await callback(entry);
    } catch (error) {
      await resetClientEntry(normalizedUrl);
      if (attempt === 1) {
        throw error;
      }
    }
  }

  throw new Error('Unable to establish an MCP session');
};

export const listLocalMcpTools = async (serverUrl: string): Promise<LocalMcpToolDefinition[]> =>
  withClientRetry(serverUrl, async ({ client }) => {
    const result = await client.listTools();
    return (result.tools || []).map((tool) => ({
      name: String(tool.name),
      description: typeof tool.description === 'string' ? tool.description : undefined,
      inputSchema:
        tool.inputSchema && typeof tool.inputSchema === 'object'
          ? (tool.inputSchema as Record<string, unknown>)
          : undefined,
    }));
  });

export const callLocalMcpTool = async (
  serverUrl: string,
  toolName: string,
  argumentsPayload: Record<string, any>
): Promise<Record<string, any>> =>
  withClientRetry(serverUrl, async ({ client }) =>
    client.callTool({
      name: toolName,
      arguments: argumentsPayload,
    })
  );
