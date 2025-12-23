import React from 'react';
import { User, Plus, Trash2, Edit2, Loader2, Settings, Wrench } from 'lucide-react';
import { ProfileSettings, ToolServer } from '../types';
import { HTTP_SERVER } from '../config';

interface ProfileButtonProps {
  onClick?: () => void;
  onSettingsChange?: (settings: ProfileSettings) => void;
  onServerAdded?: () => void;
  onServerRemoved?: () => void;
  initialSettings?: Partial<ProfileSettings>;
  username?: string;
  className?: string;
}

export const BACKEND_OPTIONS = [
  {
    value: 'openai',
    label: 'OpenAI',
    defaultUrl: 'https://api.openai.com/v1',
    models: ['gpt-5', 'gpt-5-mini', 'gpt-5-nano']
  },
  {
    value: 'livai',
    label: 'LivAI',
    defaultUrl: '',
    models: ['gpt-5', 'gpt-5-mini', 'gpt-5-nano', 'claude-sonnet-3.7']
  },
  {
    value: 'llamame',
    label: 'LLamaMe',
    defaultUrl: '',
    models: ['openai/gpt-oss-120b', 'meta-llama/Llama-3.3-70B-Instruct']
  },
  {
    value: 'alcf',
    label: 'ALCF Sophia',
    defaultUrl: '',
    models: ['openai/gpt-oss-120b', 'openai/gpt-oss-20b', 'meta-llama/Llama-4-Scout-17B-16E-Instruct']
  },
  {
    value: 'gemini',
    label: 'Google Gemini',
    defaultUrl: 'https://generativelanguage.googleapis.com/v1',
    models: ['gemini-2.0-flash-exp', 'gemini-1.5-pro', 'gemini-1.5-flash', 'gemini-1.0-pro']
  },
  {
    value: 'ollama',
    label: 'Ollama',
    defaultUrl: '',
    models: ['gpt-oss:latest', 'gpt-oss-120b', 'gpt-oss-20b']
  },
  {
    value: 'vllm',
    label: 'vLLM',
    defaultUrl: '',
    models: ['gpt-oss-120b', 'gpt-oss-20b']
  },
  {
    value: 'huggingface',
    label: 'HuggingFace Local',
    defaultUrl: '',
    models: ['']
  },
  {
    value: 'custom',
    label: 'Custom URL',
    defaultUrl: 'http://localhost:8000',
    models: ['']
  },
];

export const MOLECULE_NAME_OPTIONS = [
  { value: 'brand', label: 'Brand/Common Name' },
  { value: 'iupac', label: 'IUPAC Name' },
  { value: 'formula', label: 'Chemical Formula' },
  { value: 'smiles', label: 'SMILES' }
];

export const ProfileButton: React.FC<ProfileButtonProps> = ({
  onClick,
  onSettingsChange,
  onServerAdded,
  onServerRemoved,
  initialSettings,
  username,
  className = ''
}) => {
  const [isModalOpen, setIsModalOpen] = React.useState(false);
  const [activeTab, setActiveTab] = React.useState<'orchestrator' | 'tools'>('orchestrator');
  // Cache for storing backend-specific settings
  const [backendCache, setBackendCache] = React.useState<Record<string, {
    customUrl: string;
    model: string;
    useCustomModel: boolean;
  }>>({});

  // Default settings
  const defaultSettings: ProfileSettings = {
    backend: 'openai',
    backendLabel: 'OpenAI',
    useCustomUrl: false,
    customUrl: '',
    model: 'gpt-5-nano',
    useCustomModel: false,
    apiKey: '',
    moleculeName: 'brand',
    toolServers: [],
    ...initialSettings
  };

  const [settings, setSettings] = React.useState<ProfileSettings>(defaultSettings);
  const [tempSettings, setTempSettings] = React.useState<ProfileSettings>(settings);

  // Tool Servers state
  const [addingServer, setAddingServer] = React.useState(false);
  const [newServerUrl, setNewServerUrl] = React.useState('');
  const [editingServer, setEditingServer] = React.useState<string | null>(null);
  const [editServerUrl, setEditServerUrl] = React.useState('');
  const [connectivityStatus, setConnectivityStatus] = React.useState<Record<string, {
    status: 'checking' | 'connected' | 'disconnected';
    tools?: Array<{ name: string; description?: string }>;
    error?: string;
  }>>({});
  const [hoveredServer, setHoveredServer] = React.useState<string | null>(null);
  const [pinnedServer, setPinnedServer] = React.useState<string | null>(null);

  // Store active connections for cleanup
  const activeConnectionsRef = React.useRef<Map<string, AbortController>>(new Map());
  const tooltipRef = React.useRef<HTMLDivElement>(null);

  // Update settings when initialSettings prop changes
  React.useEffect(() => {
    if (initialSettings) {
      const backendOption = BACKEND_OPTIONS.find(opt => opt.value === initialSettings.backend);

      // Check if the model is in the predefined list for this backend
      const modelInList = backendOption?.models?.includes(initialSettings.model || '');

      const updatedSettings = {
        ...settings,
        ...initialSettings,
        backendLabel: backendOption!.label,
        // If model is not in the predefined list, automatically set useCustomModel to true
        // Unless useCustomModel is explicitly provided in initialSettings
        useCustomModel: initialSettings.useCustomModel !== undefined
          ? initialSettings.useCustomModel
          : !modelInList
      };
      setSettings(updatedSettings);
      setTempSettings(updatedSettings);
    }
  }, [initialSettings]);

  // Check connectivity for all tool servers when modal opens
  React.useEffect(() => {
    if (isModalOpen && tempSettings.toolServers && tempSettings.toolServers.length > 0) {
      tempSettings.toolServers.forEach(server => {
        checkMCPServerConnectivity(server.url);
      });
    }

    // Cleanup: abort all connections when modal closes
    return () => {
      activeConnectionsRef.current.forEach(controller => controller.abort());
      activeConnectionsRef.current.clear();
    };
  }, [isModalOpen]);

  // Handle click outside to close pinned tooltip
  React.useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (tooltipRef.current && !tooltipRef.current.contains(event.target as Node)) {
        // Check if click is on a connectivity indicator
        const target = event.target as HTMLElement;
        if (!target.closest('.connectivity-indicator')) {
          setPinnedServer(null);
        }
      }
    };

    if (pinnedServer) {
      document.addEventListener('mousedown', handleClickOutside);
      return () => document.removeEventListener('mousedown', handleClickOutside);
    }
  }, [pinnedServer]);

  const checkMCPServerConnectivity = async (url: string) => {
    // Cancel any existing connection for this URL
    const existingController = activeConnectionsRef.current.get(url);
    if (existingController) {
      existingController.abort();
    }

    // Create new abort controller
    const abortController = new AbortController();
    activeConnectionsRef.current.set(url, abortController);

    setConnectivityStatus(prev => ({
      ...prev,
      [url]: { status: 'checking' }
    }));

    try {
      // Call backend to validate the MCP server
      const response = await fetch(HTTP_SERVER + '/check-mcp-servers', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          urls: [url]
        }),
        signal: abortController.signal
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const data = await response.json();
      const result = data.results[url];

      if (result.status === 'connected') {
        setConnectivityStatus(prev => ({
          ...prev,
          [url]: {
            status: 'connected',
            tools: result.tools?.length > 0 ? result.tools : undefined
          }
        }));
      } else {
        setConnectivityStatus(prev => ({
          ...prev,
          [url]: {
            status: 'disconnected',
            error: result.error || 'Connection failed'
          }
        }));
      }

      activeConnectionsRef.current.delete(url);
    } catch (error: any) {
      if (error.name === 'AbortError') {
        // Request was cancelled, don't update status
        return;
      }

      setConnectivityStatus(prev => ({
        ...prev,
        [url]: {
          status: 'disconnected',
          error: error.message || 'Connection failed'
        }
      }));
      activeConnectionsRef.current.delete(url);
    }
  };

  const handleOpenModal = () => {
    setTempSettings(settings);
    setIsModalOpen(true);
    setActiveTab('orchestrator');
    onClick?.();
  };

  const handleSave = () => {
    setSettings(tempSettings);
    setIsModalOpen(false);
    setPinnedServer(null);
    console.log('Settings saved:', tempSettings);

    // Call the callback with the saved settings
    if (onSettingsChange) {
      onSettingsChange(tempSettings);
    }
  };

  const handleCancel = () => {
    setTempSettings(settings);
    setIsModalOpen(false);
    setAddingServer(false);
    setEditingServer(null);
    setPinnedServer(null);
    setActiveTab('orchestrator');
  };

  const handleBackendChange = (newBackend: string) => {
    const newBackendOption = BACKEND_OPTIONS.find(opt => opt.value === newBackend);

    // Cache the current settings before switching backends
    const updatedCache = {
      ...backendCache,
      [tempSettings.backend]: {
        customUrl: tempSettings.customUrl || '',
        model: tempSettings.model || '',
        useCustomModel: tempSettings.useCustomModel || false
      }
    };
    setBackendCache(updatedCache);

    // Restore cached URL and model for new backend, or use defaults
    const cached = updatedCache[newBackend];
    const urlToUse = tempSettings.useCustomUrl
      ? (cached?.customUrl || newBackendOption?.defaultUrl || '')
      : (newBackendOption?.defaultUrl || '');
    const modelToUse = cached?.model
      ? cached?.model
      : (cached?.useCustomModel
        ? (tempSettings.model || '')
        : (newBackendOption?.models[0] || ''));
    const useCustomModelToUse = cached?.useCustomModel || false;

    setTempSettings({
      ...tempSettings,
      backend: newBackend,
      customUrl: urlToUse,
      model: modelToUse,
      useCustomModel: useCustomModelToUse,
      backendLabel: newBackendOption!.label
    });
  };

  const handleCustomUrlToggle = (enabled: boolean) => {
    const backendOption = BACKEND_OPTIONS.find(opt => opt.value === tempSettings.backend);
    const cached = backendCache[tempSettings.backend];

    setTempSettings({
      ...tempSettings,
      useCustomUrl: enabled,
      // When enabling: check cache first, then current value, then default
      // When disabling: preserve the customUrl value
      customUrl: enabled
        ? (cached?.customUrl || tempSettings.customUrl || (backendOption?.defaultUrl || ''))
        : tempSettings.customUrl
    });
  };

  const handleCustomUrlChange = (newUrl: string) => {
    // Update the cache with the new custom URL
    const updatedCache = {
      ...backendCache,
      [tempSettings.backend]: {
        customUrl: newUrl,
        model: tempSettings.model,
        useCustomModel: tempSettings.useCustomModel || false
      }
    };
    setBackendCache(updatedCache);

    setTempSettings({
      ...tempSettings,
      customUrl: newUrl
    });
  };

  const handleModelSelect = (selectedModel: string) => {
    // Update the cache with the selected model
    const updatedCache = {
      ...backendCache,
      [tempSettings.backend]: {
        customUrl: tempSettings.customUrl || '',
        model: selectedModel,
        useCustomModel: tempSettings.useCustomModel || false
      }
    };
    setBackendCache(updatedCache);

    setTempSettings({
      ...tempSettings,
      model: selectedModel
    });
  };

  const handleCustomModelToggle = (enabled: boolean) => {
    const backendOption = BACKEND_OPTIONS.find(opt => opt.value === tempSettings.backend);
    const cached = backendCache[tempSettings.backend];

    let modelToUse: string;

    if (enabled) {
      // When enabling custom model, keep current model
      modelToUse = tempSettings.model;
    } else {
      // When disabling custom model, find a valid preset model
      // First check if cached model is in the preset list
      const cachedModelInList = cached?.model && backendOption?.models?.includes(cached.model);
      // Then check if current temp model is in the preset list
      const tempModelInList = backendOption?.models?.includes(tempSettings.model);

      if (cachedModelInList) {
        // Use cached model if it's a valid preset
        modelToUse = cached.model;
      } else if (tempModelInList) {
        // Use current temp model if it's a valid preset
        modelToUse = tempSettings.model;
      } else {
        // Fall back to first model in the preset list
        modelToUse = backendOption?.models?.[0] || tempSettings.model;
      }
    }

    const updatedCache = {
      ...backendCache,
      [tempSettings.backend]: {
        customUrl: tempSettings.customUrl || '',
        model: modelToUse,
        useCustomModel: enabled
      }
    };
    setBackendCache(updatedCache);

    setTempSettings({
      ...tempSettings,
      useCustomModel: enabled,
      model: modelToUse
    });
  };

  const handleCustomModelChange = (newModel: string) => {
    // Update the cache with the custom model
    const updatedCache = {
      ...backendCache,
      [tempSettings.backend]: {
        customUrl: tempSettings.customUrl || '',
        model: newModel,
        useCustomModel: tempSettings.useCustomModel || false
      }
    };
    setBackendCache(updatedCache);

    setTempSettings({
      ...tempSettings,
      model: newModel
    });
  };

  // Tool Server handlers
  const handleAddServer = async () => {
    if (!newServerUrl.trim()) return;

    const url = newServerUrl.trim();

    // Set to checking state
    setConnectivityStatus(prev => ({
      ...prev,
      [url]: { status: 'checking' }
    }));

    try {
      // Call backend to validate and register
      const response = await fetch(HTTP_SERVER + '/validate-mcp-server', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          url: url,
          name: `Server ${Date.now()}`
        })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const result = await response.json();

      if (result.status === 'connected') {
        // Server validated and registered successfully
        const newServer: ToolServer = {
          id: Date.now().toString(),
          url: result.url || url
        };

        setTempSettings({
          ...tempSettings,
          toolServers: [...(tempSettings.toolServers || []), newServer]
        });

        // Update connectivity status with tools
        setConnectivityStatus(prev => ({
          ...prev,
          [newServer.url]: {
            status: 'connected',
            tools: result.tools?.length > 0 ? result.tools : undefined
          }
        }));

        setNewServerUrl('');
        setAddingServer(false);

        // Notify parent that server was added
        if (onServerAdded) {
           onServerAdded();
        }
      } else {
        // Validation failed
        setConnectivityStatus(prev => ({
          ...prev,
          [url]: {
            status: 'disconnected',
            error: result.error || 'Validation failed'
          }
        }));

        // Show error to user but keep the input visible
        alert(`Failed to add server: ${result.error || 'Validation failed'}`);
      }
    } catch (error: any) {
      setConnectivityStatus(prev => ({
        ...prev,
        [url]: {
          status: 'disconnected',
          error: error.message || 'Connection failed'
        }
      }));

      alert(`Failed to add server: ${error.message || 'Connection failed'}`);
    }
  };

  const handleEditServer = (serverId: string) => {
    const server = tempSettings.toolServers?.find(s => s.id === serverId);
    if (server) {
      setEditingServer(serverId);
      setEditServerUrl(server.url);
    }
  };

  const handleSaveServerEdit = async (serverId: string) => {
    if (!editServerUrl.trim()) return;

    const oldServer = tempSettings.toolServers?.find(s => s.id === serverId);
    if (!oldServer) return;

    const newUrl = editServerUrl.trim();

    // Set to checking state
    setConnectivityStatus(prev => ({
      ...prev,
      [newUrl]: { status: 'checking' }
    }));

    try {
      // Call backend to validate and register
      const response = await fetch(HTTP_SERVER + '/validate-mcp-server', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          url: newUrl,
          name: oldServer.name || `Server ${serverId}`
        })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const result = await response.json();

      if (result.status === 'connected') {
        // Clean up old URL status
        if (oldServer.url !== newUrl) {
          setConnectivityStatus(prev => {
            const newStatus = { ...prev };
            delete newStatus[oldServer.url];
            return newStatus;
          });
        }

        // Update server URL
        setTempSettings({
          ...tempSettings,
          toolServers: tempSettings.toolServers?.map(s =>
            s.id === serverId ? { ...s, url: result.url || newUrl } : s
          ) || []
        });

        // Update connectivity status
        setConnectivityStatus(prev => ({
          ...prev,
          [result.url || newUrl]: {
            status: 'connected',
            tools: result.tools?.length > 0 ? result.tools : undefined
          }
        }));

        setEditingServer(null);
        setEditServerUrl('');
      } else {
        // Validation failed
        setConnectivityStatus(prev => ({
          ...prev,
          [newUrl]: {
            status: 'disconnected',
            error: result.error || 'Validation failed'
          }
        }));

        alert(`Failed to update server: ${result.error || 'Validation failed'}`);
      }
    } catch (error: any) {
      setConnectivityStatus(prev => ({
        ...prev,
        [newUrl]: {
          status: 'disconnected',
          error: error.message || 'Connection failed'
        }
      }));

      alert(`Failed to update server: ${error.message || 'Connection failed'}`);
    }
  };

  const handleDeleteServer = async (serverId: string) => {
    const server = tempSettings.toolServers?.find(s => s.id === serverId);
    if (!server) return;

    console.log('🗑️ Deleting server:', server.url);

    try {
      const response = await fetch(HTTP_SERVER + '/delete-mcp-server', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ url: server.url })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }

      const result = await response.json();
      console.log('🗑️ Backend response:', result);

      if (result.status === 'deleted' || result.status === 'not_found') {
        // Cancel any active connection
        const controller = activeConnectionsRef.current.get(server.url);
        if (controller) {
          controller.abort();
          activeConnectionsRef.current.delete(server.url);
        }

        // Clean up status
        setConnectivityStatus(prev => {
          const newStatus = { ...prev };
          delete newStatus[server.url];
          return newStatus;
        });

        if (pinnedServer === serverId) {
          setPinnedServer(null);
        }

        // Remove from tempSettings
        setTempSettings({
          ...tempSettings,
          toolServers: tempSettings.toolServers?.filter(s => s.id !== serverId) || []
        });

        // Notify parent
        if (onServerRemoved) {
          onServerRemoved();
        }

        console.log('✅ Server deleted successfully');
      } else {
        alert(`Failed to delete server: ${result.message}`);
      }
    } catch (error: any) {
      alert(`Failed to delete server: ${error.message || 'Connection failed'}`);
    }
  };

  const handleClearAllServers = async () => {
    if (!tempSettings.toolServers || tempSettings.toolServers.length === 0) {
      return;
    }

    if (!confirm(`Remove all ${tempSettings.toolServers.length} tool servers?`)) {
      return;
    }

    // Delete each server from backend
    const deletePromises = tempSettings.toolServers.map(async (server) => {
      try {
        await fetch(HTTP_SERVER + '/delete-mcp-server', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ url: server.url })
        });
      } catch (error) {
        console.warn(`Error deleting ${server.url}:`, error);
      }
    });

    await Promise.all(deletePromises);

    // Cancel all active connections
    activeConnectionsRef.current.forEach(controller => controller.abort());
    activeConnectionsRef.current.clear();

    // Clear frontend state
    setTempSettings({
      ...tempSettings,
      toolServers: []
    });
    setConnectivityStatus({});
    setPinnedServer(null);

    if (onServerRemoved) {
      onServerRemoved();
    }
  };

  const currentBackendOption = BACKEND_OPTIONS.find(opt => opt.value === tempSettings.backend);

  return (
    <>
      <button
        onClick={handleOpenModal}
        className={`btn btn-secondary btn-sm ${className}`}
      >
        <User className="w-4 h-4" />
        <span>Profile</span>
        <svg className="w-3 h-3 ml-1" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Profile Settings Modal */}
      {isModalOpen && (
        <div className="modal-overlay">
          <div className="modal-content modal-content-lg">
            <div className="modal-header">
              <div>
                <h2 className="modal-title">{username}'s Profile Settings</h2>
                <p className="modal-subtitle">Configure your connection and tools</p>
              </div>
              <button
                onClick={handleCancel}
                className="btn-icon"
              >
                <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
                </svg>
              </button>
            </div>

            {/* Tab Navigation */}
            <div className="border-b border-border">
              <div className="flex gap-1 px-6">
                <button
                  onClick={() => setActiveTab('orchestrator')}
                  className={`flex items-center gap-2 px-4 py-3 border-b-2 transition-colors ${
                    activeTab === 'orchestrator'
                      ? 'border-primary text-primary font-medium'
                      : 'border-transparent text-muted hover:text-secondary'
                  }`}
                >
                  <Settings className="w-4 h-4" />
                  <span>Orchestrator</span>
                </button>
                <button
                  onClick={() => setActiveTab('tools')}
                  className={`flex items-center gap-2 px-4 py-3 border-b-2 transition-colors ${
                    activeTab === 'tools'
                      ? 'border-primary text-primary font-medium'
                      : 'border-transparent text-muted hover:text-secondary'
                  }`}
                >
                  <Wrench className="w-4 h-4" />
                  <span>Tool Servers</span>
                  {tempSettings.toolServers && tempSettings.toolServers.length > 0 && (
                    <span className="bg-primary text-white text-xs px-2 py-0.5 rounded-full">
                      {tempSettings.toolServers.length}
                    </span>
                  )}
                </button>
              </div>
            </div>

            <div className="modal-body space-y-4">
              {/* Orchestrator Tab */}
              {activeTab === 'orchestrator' && (
                <div className="space-y-4">
                {/* Backend Selector */}
                <div className="form-group">
                  <label className="form-label">
                    Backend
                  </label>
                  <select
                    value={tempSettings.backend}
                    onChange={(e) => handleBackendChange(e.target.value)}
                    className="form-select"
                  >
                    {BACKEND_OPTIONS.map(option => (
                      <option key={option.value} value={option.value}>
                        {option.label}
                      </option>
                    ))}
                  </select>
                </div>

                {/* Use Custom URL Checkbox */}
                <div>
                  <label className="form-label flex items-center gap-2 cursor-pointer">
                    <input
                      type="checkbox"
                      checked={tempSettings.useCustomUrl}
                      onChange={(e) => handleCustomUrlToggle(e.target.checked)}
                      className="form-checkbox"
                    />
                    <span>Use custom URL for this backend</span>
                  </label>
                  <p className="text-xs text-muted mt-1 ml-6">
                    Override the default endpoint with a custom server URL
                  </p>
                </div>

                {/* Custom URL Field (conditional) */}
                {tempSettings.useCustomUrl && (
                  <div className="form-group animate-fadeIn">
                    <label className="form-label">
                      Custom URL
                      <span className="text-muted text-xs ml-2">
                        {tempSettings.backend === 'vllm' && '(vLLM endpoint)'}
                        {tempSettings.backend === 'ollama' && '(Ollama endpoint)'}
                        {(tempSettings.backend === 'livai') && '(LivAI base URL)'}
                        {(tempSettings.backend === 'llamame') && '(LLamaMe base URL)'}
                        {(tempSettings.backend === 'alcf') && '(ALCF Sophia base URL)'}
                        {tempSettings.backend === 'openai' && '(OpenAI-compatible endpoint)'}
                        {tempSettings.backend === 'gemini' && '(Gemini API endpoint)'}
                      </span>
                    </label>
                    <input
                      type="text"
                      value={tempSettings.customUrl || ''}
                      onChange={(e) => handleCustomUrlChange(e.target.value)}
                      placeholder={currentBackendOption?.defaultUrl || 'http://localhost:8000'}
                      className="form-input"
                    />
                    <p className="text-xs text-muted mt-1">
                      Default: {currentBackendOption?.defaultUrl || 'Not set'}
                    </p>
                  </div>
                )}

                {/* Model Selection */}
                <div className="form-group">
                  <label className="form-label">
                    Model
                    <span className="text-muted text-xs ml-2">
                      {tempSettings.backend === 'openai' && '(GPT models)'}
                      {tempSettings.backend === 'livai' && '(LLNL Enterprise models)'}
                      {tempSettings.backend === 'llamame' && '(LLNL Internal models)'}
                      {tempSettings.backend === 'alcf' && '(ACLF Internal models)'}
                      {tempSettings.backend === 'gemini' && '(Gemini models)'}
                      {tempSettings.backend === 'ollama' && '(Local models)'}
                      {tempSettings.backend === 'vllm' && '(vLLM models)'}
                    </span>
                  </label>

                  {!tempSettings.useCustomModel ? (
                    <select
                      value={tempSettings.model}
                      onChange={(e) => handleModelSelect(e.target.value)}
                      className="form-select text-mono"
                    >
                      {currentBackendOption?.models?.map(model => (
                        <option key={model} value={model}>
                          {model}
                        </option>
                      ))}
                    </select>
                  ) : (
                    <input
                      type="text"
                      value={tempSettings.model}
                      onChange={(e) => handleCustomModelChange(e.target.value)}
                      placeholder="Enter custom model name"
                      className="form-input text-mono"
                    />
                  )}
                </div>

                {/* Use Custom Model Checkbox */}
                <div>
                  <label className="form-label flex items-center gap-2 cursor-pointer">
                    <input
                      type="checkbox"
                      checked={tempSettings.useCustomModel || false}
                      onChange={(e) => handleCustomModelToggle(e.target.checked)}
                      className="form-checkbox"
                    />
                    <span>Use custom model name</span>
                  </label>
                  <p className="text-xs text-muted mt-1 ml-6">
                    Enter a custom model identifier not in the preset list
                  </p>
                </div>

                {/* API Key Field */}
                <div className="form-group">
                  <label className="form-label">
                    API Key
                    <span className="text-muted text-xs ml-2">
                      {(tempSettings.backend === 'ollama' || tempSettings.backend === 'huggingface' || tempSettings.backend === 'vllm') && '(Optional for local backends)'}
                      {tempSettings.backend === 'openai' && '(OPENAI_API_KEY)'}
                      {tempSettings.backend === 'livai' && '(LIVAI_API_KEY)'}
                      {tempSettings.backend === 'llamame' && '(LLAMAME_API_KEY)'}
                      {tempSettings.backend === 'alcf' && '(ALCF_API_KEY)'}
                      {tempSettings.backend === 'gemini' && '(GOOGLE_API_KEY)'}
                    </span>
                  </label>
                  <input
                    type="password"
                    value={tempSettings.apiKey}
                    onChange={(e) => setTempSettings({...tempSettings, apiKey: e.target.value})}
                    placeholder="Enter your API key"
                    className="form-input text-mono"
                  />
                </div>

                {/* Divider */}
                <div className="border-t border-border my-4"></div>

                {/* Molecule Name Preference */}
                <div className="form-group">
                  <label className="form-label">
                    Preferred Molecule Name Format
                  </label>
                  <select
                    value={tempSettings.moleculeName || 'brand'}
                    onChange={(e) => setTempSettings({...tempSettings, moleculeName: e.target.value as any})}
                    className="form-select"
                  >
                    {MOLECULE_NAME_OPTIONS.map(option => (
                      <option key={option.value} value={option.value}>
                        {option.label}
                      </option>
                    ))}
                  </select>
                  <p className="text-xs text-muted mt-1">
                    Choose how molecule names are displayed throughout the application
                  </p>
                </div>
                </div>
              )}

              {/* Tool Servers Tab */}
              {activeTab === 'tools' && (
                <div className="space-y-4">
                  <div className="flex-between">
                    <div>
                      <h3 className="text-lg font-medium text-primary">Custom Tool Servers (MCP)</h3>
                      <p className="text-xs text-muted mt-1">
                        Configure external MCP tool servers for extended functionality
                      </p>
                    </div>
                  {tempSettings.toolServers && tempSettings.toolServers.length > 0 && (
                    <button
                      onClick={handleClearAllServers}
                      className="btn btn-tertiary btn-sm text-xs"
                    >
                      <Trash2 className="w-3 h-3" />
                      Clear All
                    </button>
                  )}
                </div>

                {/* Server List */}
                <div className="space-y-2">
                  {tempSettings.toolServers?.map(server => (
                    <div key={server.id} className="border border-border rounded p-3 bg-surface-secondary hover:bg-surface-tertiary transition-colors">
                      {editingServer === server.id ? (
                        <div className="space-y-2">
                          <input
                            type="text"
                            value={editServerUrl}
                            onChange={(e) => setEditServerUrl(e.target.value)}
                            onKeyDown={(e) => {
                              if (e.key === 'Enter') handleSaveServerEdit(server.id);
                              if (e.key === 'Escape') {
                                setEditingServer(null);
                                setEditServerUrl('');
                              }
                            }}
                            placeholder="https://example.com/sse"
                            className="form-input text-sm text-mono"
                            autoFocus
                          />
                          <div className="flex gap-2">
                            <button
                              onClick={() => handleSaveServerEdit(server.id)}
                              className="btn btn-secondary btn-sm text-xs flex-1"
                            >
                              Save
                            </button>
                            <button
                              onClick={() => {
                                setEditingServer(null);
                                setEditServerUrl('');
                              }}
                              className="btn btn-tertiary btn-sm text-xs flex-1"
                            >
                              Cancel
                            </button>
                          </div>
                        </div>
                      ) : (
                        <div className="flex items-center gap-3">
                          {/* Connectivity Indicator */}
                          <div
                            className="relative flex-shrink-0 connectivity-indicator"
                            onMouseEnter={() => setHoveredServer(server.id)}
                            onMouseLeave={() => setHoveredServer(null)}
                            onClick={() => setPinnedServer(pinnedServer === server.id ? null : server.id)}
                            style={{ cursor: 'pointer' }}
                          >
                            {connectivityStatus[server.url]?.status === 'checking' ? (
                              <Loader2 className="w-4 h-4 text-muted animate-spin" />
                            ) : (
                              <div
                                className={`w-3 h-3 rounded-full ${
                                  connectivityStatus[server.url]?.status === 'connected'
                                    ? 'bg-green-500'
                                    : 'bg-red-500'
                                }`}
                                title={
                                  connectivityStatus[server.url]?.status === 'connected'
                                    ? 'Connected (click for tools)'
                                    : connectivityStatus[server.url]?.error || 'Disconnected'
                                }
                              />
                            )}

                            {/* Tools Tooltip */}
                            {(hoveredServer === server.id || pinnedServer === server.id) &&
                             connectivityStatus[server.url]?.status === 'connected' &&
                             connectivityStatus[server.url]?.tools &&
                             connectivityStatus[server.url].tools!.length > 0 && (
                              <div ref={tooltipRef} className="ws-tooltip" style={{ left: 'auto', right: 0, top: '2.5rem' }}>
                                <p className="text-xs font-semibold mb-2 text-primary">Available Tools:</p>
                                <ul className="text-xs space-y-1.5 max-h-[200px] overflow-y-auto custom-scrollbar">
                                  {connectivityStatus[server.url].tools!.map((tool, idx) => (
                                    <li key={idx} className="text-secondary">
                                      <span className="font-mono text-primary font-semibold">• {tool.name}</span>
                                      {tool.description && (
                                        <p className="text-secondary text-[10px] ml-3 mt-0.5">
                                          {tool.description.length > 80
                                            ? `${tool.description.substring(0, 80)}...`
                                            : tool.description}
                                        </p>
                                      )}
                                    </li>
                                  ))}
                                </ul>
                              </div>
                            )}
                          </div>

                          {/* URL */}
                          <div className="flex-1 min-w-0">
                            <p className="text-sm text-mono truncate text-primary">{server.url}</p>
                          </div>

                          {/* Action Buttons */}
                          <div className="flex items-center gap-2 flex-shrink-0">
                            <button
                              onClick={() => handleEditServer(server.id)}
                              className="btn-icon p-1.5"
                              title="Edit server"
                            >
                              <Edit2 className="w-4 h-4" />
                            </button>
                            <button
                              onClick={() => handleDeleteServer(server.id)}
                              className="btn-icon p-1.5 hover:text-red-500"
                              title="Delete server"
                            >
                              <Trash2 className="w-4 h-4" />
                            </button>
                          </div>
                        </div>
                      )}
                    </div>
                  ))}

                  {/* Empty State */}
                  {(!tempSettings.toolServers || tempSettings.toolServers.length === 0) && !addingServer && (
                    <div className="border border-border border-dashed rounded p-8 text-center">
                      <Wrench className="w-12 h-12 text-muted mx-auto mb-3" />
                      <p className="text-sm text-muted font-medium">No tool servers configured</p>
                      <p className="text-xs text-tertiary mt-1">Add MCP servers to extend functionality</p>
                    </div>
                  )}
                </div>

                {/* Add New Server */}
                {addingServer ? (
                  <div className="border border-primary rounded p-3 bg-surface-secondary space-y-2">
                    <label className="text-xs text-muted font-medium">MCP Server URL</label>
                    <input
                      type="text"
                      value={newServerUrl}
                      onChange={(e) => setNewServerUrl(e.target.value)}
                      onKeyDown={(e) => {
                        if (e.key === 'Enter') handleAddServer();
                        if (e.key === 'Escape') {
                          setAddingServer(false);
                          setNewServerUrl('');
                        }
                      }}
                      placeholder="https://example.com/sse"
                      className="form-input text-sm text-mono"
                      autoFocus
                    />
                    <div className="flex gap-2">
                      <button
                        onClick={handleAddServer}
                        className="btn btn-secondary btn-sm flex-1 text-xs"
                      >
                        <Plus className="w-3 h-3" />
                        Add Server
                      </button>
                      <button
                        onClick={() => {
                          setAddingServer(false);
                          setNewServerUrl('');
                        }}
                        className="btn btn-tertiary btn-sm flex-1 text-xs"
                      >
                        Cancel
                      </button>
                    </div>
                  </div>
                ) : (
                  <button
                    onClick={() => setAddingServer(true)}
                    className="btn btn-secondary btn-sm w-full group"
                  >
                    <Plus className="w-4 h-4 group-hover:scale-110 transition-transform" />
                    <span>Add Tool Server</span>
                  </button>
                )}
              </div>
              )}
            </div>

            <div className="modal-footer">
              <button
                onClick={handleSave}
                className="btn btn-primary flex-1"
              >
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                </svg>
                Save Settings
              </button>
              <button
                onClick={handleCancel}
                className="btn btn-tertiary"
              >
                Cancel
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
};
