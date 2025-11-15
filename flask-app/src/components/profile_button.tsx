import React from 'react';
import { User } from 'lucide-react';

interface ProfileButtonProps {
  onClick?: () => void;
  className?: string;
}

export const ProfileButton: React.FC<ProfileButtonProps> = ({
  onClick,
  className = ''
}) => {
  const [isModalOpen, setIsModalOpen] = React.useState(false);
  const [settings, setSettings] = React.useState<ProfileSettings>({
    url: 'http://localhost:8000',
    orchestrator: 'default',
    model: 'claude-sonnet-4-20250514',
    apiKey: ''
  });

  const [tempSettings, setTempSettings] = React.useState<ProfileSettings>(settings);

  const handleOpenModal = () => {
    setTempSettings(settings);
    setIsModalOpen(true);
    onClick?.();
  };

  const handleSave = () => {
    setSettings(tempSettings);
    setIsModalOpen(false);
    console.log('Settings saved:', tempSettings);
  };

  const handleCancel = () => {
    setTempSettings(settings);
    setIsModalOpen(false);
  };

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
              {/* URL Field */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Server URL
                </label>
                <input
                  type="text"
                  value={tempSettings.url}
                  onChange={(e) => setTempSettings({...tempSettings, url: e.target.value})}
                  placeholder="http://localhost:8000"
                  className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50"
                />
              </div>

              {/* Orchestrator Field */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Orchestrator
                </label>
                <input
                  type="text"
                  value={tempSettings.orchestrator}
                  onChange={(e) => setTempSettings({...tempSettings, orchestrator: e.target.value})}
                  placeholder="default"
                  className="w-full px-4 py-3 bg-white/10 border-2 border-purple-400/50 rounded-lg focus:border-purple-400 focus:outline-none text-white placeholder-purple-300/50"
                />
              </div>

              {/* Model Field */}
              <div>
                <label className="block text-sm font-medium text-purple-200 mb-2">
                  Model
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
