
import React from 'react';
import { User } from 'lucide-react';

interface ProfileButtonProps {
  onClick?: () => void;
  onSettingsChange?: (settings: ProfileSettings) => void;
  initialSettings?: Partial<ProfileSettings>;
  className?: string;
}

interface ProfileSettings {
  backend: string;
  customUrl?: string;
  model: string;
  apiKey: string;
}

const BACKEND_OPTIONS = [
  { value: 'openai', label: 'OpenAI' },
  { value: 'gemini', label: 'Google Gemini' },
  { value: 'livai', label: 'LivAI (LivChat)' },
  { value: 'ollama', label: 'Ollama' },
  { value: 'huggingface', label: 'HuggingFace Local' },
  { value: 'vllm', label: 'vLLM' },
  { value: 'custom', label: 'Custom URL' },
];

const BACKENDS_REQUIRING_URL = ['livai', 'vllm', 'ollama', 'custom'];

export const ProfileButton: React.FC<ProfileButtonProps> = ({
  onClick,
  onSettingsChange,
  initialSettings,
  className = ''
}) => {
  const [isModalOpen, setIsModalOpen] = React.useState(false);

  // Merge provided initial settings with defaults
  const defaultSettings: ProfileSettings = {
    backend: 'openai',
    customUrl: '',
    model: 'gpt-5-nano',
    apiKey: '',
    ...initialSettings
  };

  const [settings, setSettings] = React.useState<ProfileSettings>(defaultSettings);
  const [tempSettings, setTempSettings] = React.useState<ProfileSettings>(settings);

  // Update settings when initialSettings prop changes
  React.useEffect(() => {
    if (initialSettings) {
      const updatedSettings = {
        ...settings,
        ...initialSettings
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
    setTempSettings({
      ...tempSettings,
      backend: newBackend,
      // Clear custom URL if switching away from a backend that needs it
      customUrl: BACKENDS_REQUIRING_URL.includes(newBackend) ? tempSettings.customUrl : ''
    });
  };

  const showUrlField = BACKENDS_REQUIRING_URL.includes(tempSettings.backend);

  return (
    <>
      <button
        onClick={handleOpenModal}
        className={`px-4 py-2 bg-purple-500/30 text-white rounded-lg text-sm font-semibold hover:bg-purple-500/50 transition-all flex items-center gap-2 ${className}`}
      >
        <User className="w-4 h-4" />
        <span>Profile</span>
        <svg className="w-3 h-3 ml-1" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>

      {/* Profile Settings Modal */}
      {isModalOpen && (
        <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-50 flex items-center justify-center p-4">
          <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-2xl w-full p-6">
            <div className="flex items-center justify-between mb-6">
              <div>
                <h2 className="text-xl font-bold text-white">Profile Settings</h2>
                <p className="text-sm text-purple-300">Configure your connection and model settings</p>
              </div>
              <button
                onClick={handleCancel}
                className="text-purple-300 hover:text-white transition-colors"
              >
              </button>
            </div>

            <div className="space-y-4">
              {/* Backend Selector */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Backend
                </label>
                <select
                  value={tempSettings.backend}
                  onChange={(e) => handleBackendChange(e.target.value)}
                  className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white cursor-pointer"
                >
                  {BACKEND_OPTIONS.map(option => (
                    <option key={option.value} value={option.value} className="bg-slate-800">
                      {option.label}
                    </option>
                  ))}
                </select>
              </div>

              {/* Custom URL Field (conditional) */}
              {showUrlField && (
                <div>
                  <label className="block text-sm font-medium text-purple-200 mb-2">
                    {tempSettings.backend === 'custom' ? 'Custom URL' : 'Server URL'}
                    <span className="text-purple-400 text-xs ml-2">
                      {tempSettings.backend === 'vllm' && '(vLLM endpoint)'}
                      {tempSettings.backend === 'ollama' && '(Ollama endpoint)'}
                      {(tempSettings.backend === 'livai') && '(LivAI base URL)'}
                    </span>
                  </label>
                  <input
                    type="text"
                    value={tempSettings.customUrl || ''}
                    onChange={(e) => setTempSettings({...tempSettings, customUrl: e.target.value})}
                    placeholder={
                      tempSettings.backend === 'vllm' ? 'http://localhost:8000/v1' :
                      tempSettings.backend === 'ollama' ? 'http://localhost:11434' :
                      'http://localhost:8000'
                    }
                    className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50"
                  />
                </div>
              )}

              {/* Model Field */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Model
                  <span className="text-purple-400 text-xs ml-2">
                    {tempSettings.backend === 'openai' && '(e.g., gpt-4, gpt-5)'}
                    {tempSettings.backend === 'gemini' && '(e.g., gemini-flash-latest)'}
                    {tempSettings.backend === 'ollama' && '(e.g., gpt-oss:latest)'}
                    {tempSettings.backend === 'vllm' && '(e.g., gpt-oss)'}
                  </span>
                </label>
                <input
                  type="text"
                  value={tempSettings.model}
                  onChange={(e) => setTempSettings({...tempSettings, model: e.target.value})}
                  placeholder="claude-sonnet-4-20250514"
                  className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 font-mono text-sm"
                />
              </div>

              {/* API Key Field */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  API Key
                  <span className="text-purple-400 text-xs ml-2">
                    {(tempSettings.backend === 'ollama' || tempSettings.backend === 'huggingface' || tempSettings.backend === 'vllm') && '(Optional for local backends)'}
                    {(tempSettings.backend === 'openai' || tempSettings.backend === 'livai') && '(OPENAI_API_KEY)'}
                    {tempSettings.backend === 'gemini' && '(GOOGLE_API_KEY)'}
                  </span>
                </label>
                <input
                  type="password"
                  value={tempSettings.apiKey}
                  onChange={(e) => setTempSettings({...tempSettings, apiKey: e.target.value})}
                  placeholder="Enter your API key"
                  className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50 font-mono text-sm"
                />
              </div>
            </div>

            <div className="flex gap-3 mt-6">
              <button
                onClick={handleSave}
                className="flex-1 px-6 py-3 bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded-lg font-semibold shadow-lg hover:shadow-xl transform hover:-translate-y-0.5 transition-all flex items-center justify-center gap-2"
              >
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                </svg>
                Save Settings
              </button>
              <button
                onClick={handleCancel}
                className="px-6 py-3 bg-white/20 text-white rounded-lg font-semibold hover:bg-white/30 transition-all"
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
