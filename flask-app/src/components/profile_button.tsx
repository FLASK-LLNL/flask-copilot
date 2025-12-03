import React from 'react';
import { User } from 'lucide-react';
import { ProfileSettings } from '../types';

interface ProfileButtonProps {
  onClick?: () => void;
  onSettingsChange?: (settings: ProfileSettings) => void;
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

export const ProfileButton: React.FC<ProfileButtonProps> = ({
  onClick,
  onSettingsChange,
  initialSettings,
  username,
  className = ''
}) => {
  const [isModalOpen, setIsModalOpen] = React.useState(false);
  // Cache for storing backend-specific settings
  const [backendCache, setBackendCache] = React.useState<Record<string, {
    customUrl: string;
    model: string;
    useCustomModel: boolean;
  }>>({});

  // Merge provided initial settings with defaults
  const defaultSettings: ProfileSettings = {
    backend: 'openai',
    useCustomUrl: false,
    customUrl: '',
    model: 'gpt-5-nano',
    useCustomModel: false,
    apiKey: '',
    ...initialSettings
  };

  const [settings, setSettings] = React.useState<ProfileSettings>(defaultSettings);
  const [tempSettings, setTempSettings] = React.useState<ProfileSettings>(settings);

  // Update settings when initialSettings prop changes
  React.useEffect(() => {
    if (initialSettings) {
      const backendOption = BACKEND_OPTIONS.find(opt => opt.value === initialSettings.backend);

      // Check if the model is in the predefined list for this backend
      const modelInList = backendOption?.models?.includes(initialSettings.model || '');

      const updatedSettings = {
        ...settings,
        ...initialSettings,
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

  const handleOpenModal = () => {
    setTempSettings(settings);
    setIsModalOpen(true);
    onClick?.();
  };

  const handleSave = () => {
    setSettings(tempSettings);
    setIsModalOpen(false);
    console.log('Settings saved:', tempSettings);

    // Call the callback with the saved settings
    if (onSettingsChange) {
      onSettingsChange(tempSettings);
    }
  };

  const handleCancel = () => {
    setTempSettings(settings);
    setIsModalOpen(false);
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
      useCustomModel: useCustomModelToUse
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
                <p className="modal-subtitle">Configure your connection and model settings</p>
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

            <div className="modal-body space-y-4">
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
