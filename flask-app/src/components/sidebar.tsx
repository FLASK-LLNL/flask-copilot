// Reasoning sidebar component

import { useRef, useEffect, useState } from "react";
import { SidebarMessage, SidebarProps, SidebarState, VisibleSources } from "../types";
import { X } from "lucide-react";
import { MarkdownText } from "./markdown";
import { MoleculeSVG } from "./molecule";

export const useSidebarState = (): SidebarState => {
    const [messages, setMessages] = useState<SidebarMessage[]>([]);
    const [sourceFilterOpen, setSourceFilterOpen] = useState<boolean>(false);
    const [visibleSources, setVisibleSources] = useState<VisibleSources>({
        'System': true,
        'Reasoning': true,
        'Logger (Error)': true,
        'Logger (Warning)': true,
        'Logger (Info)': false,
        'Logger (Debug)': false
    });


    return {messages, setMessages, sourceFilterOpen, setSourceFilterOpen, visibleSources, setVisibleSources };
};

export const ReasoningSidebar: React.FC<SidebarProps> = ({messages, rdkitModule, setSidebarOpen, sourceFilterOpen, setSourceFilterOpen, visibleSources, setVisibleSources}) => {
    const sidebarRef = useRef<HTMLDivElement>(null);

    // Auto-scroll sidebar to bottom when new messages arrive
    useEffect(() => {
        if (sidebarRef.current && messages.length > 0) {
            // Delay scroll to allow molecules to render
            setTimeout(() => {
            if (sidebarRef.current) {
                sidebarRef.current.scrollTo({
                top: sidebarRef.current.scrollHeight,
                behavior: 'smooth'
                });
            }
            }, 100);
        }
    }, [messages]);


    return (
      <div className="w-96 bg-slate-900 border-l-2 border-purple-400 shadow-2xl flex-shrink-0 sticky top-0 h-screen overflow-hidden animate-slideInSidebar">
        <div className="flex flex-col h-full">
          <div className="flex items-center justify-between p-4 border-b border-purple-400/30">
            <h3 className="text-lg font-semibold text-white">Reasoning</h3>
            <div className="flex items-center gap-2">
              <div className="relative">
                <button 
                  onClick={(e) => { e.stopPropagation(); setSourceFilterOpen(!sourceFilterOpen); }}
                  onMouseDown={(e) => e.stopPropagation()}
                  className="px-3 py-1 bg-purple-500/30 text-white rounded-lg text-xs font-semibold hover:bg-purple-500/50 transition-all flex items-center gap-2"
                >
                  <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z" />
                  </svg>
                  Filter
                  <svg className={`w-3 h-3 transition-transform ${sourceFilterOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                  </svg>
                </button>
                {sourceFilterOpen && (
                  <div className="absolute top-full mt-2 right-0 bg-slate-800 border-2 border-purple-400 rounded-lg shadow-2xl py-2 min-w-48 z-[100]" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                    <div className="px-3 py-2 border-b border-purple-400/30 text-xs font-semibold text-purple-300">
                      Message Sources
                    </div>
                    {Object.keys(visibleSources).map(source => (
                      <label key={source} className="flex items-center gap-2 px-3 py-2 hover:bg-purple-600/30 cursor-pointer transition-colors">
                        <input
                          type="checkbox"
                          checked={visibleSources[source]}
                          onChange={() => setVisibleSources(prev => ({ ...prev, [source]: !prev[source] }))}
                          className="w-4 h-4 rounded border-purple-400/50 bg-white/20 text-purple-600 focus:ring-purple-500 focus:ring-offset-0"
                        />
                        <span className="text-sm text-white">{source}</span>
                        <span className="ml-auto text-xs text-purple-400">
                          ({messages.filter(m => m.source === source).length})
                        </span>
                      </label>
                    ))}
                    <div className="px-3 py-2 border-t border-purple-400/30 flex gap-2">
                      <button
                        onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: true }), {}))}
                        className="flex-1 px-2 py-1 bg-purple-500/30 text-white rounded text-xs font-semibold hover:bg-purple-500/50 transition-all"
                      >
                        All
                      </button>
                      <button
                        onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: false }), {}))}
                        className="flex-1 px-2 py-1 bg-purple-500/30 text-white rounded text-xs font-semibold hover:bg-purple-500/50 transition-all"
                      >
                        None
                      </button>
                    </div>
                  </div>
                )}
              </div>
              <button onClick={() => setSidebarOpen(false)} className="text-purple-300 hover:text-white transition-colors">
                <X className="w-5 h-5" />
              </button>
            </div>
          </div>

          <div className="flex-1 overflow-y-auto p-4 space-y-3" ref={sidebarRef} style={{ overflowX: 'hidden' }}>
            {messages.filter(msg => visibleSources[msg.source]).length === 0 ? (
              <div className="text-center text-purple-400 mt-8">
                <svg className="w-12 h-12 mx-auto mb-2 opacity-50" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z" />
                </svg>
                <p className="text-sm">
                  {messages.length === 0 ? 'No messages yet' : 'No messages match the selected filters'}
                </p>
              </div>
            ) : (
              messages.filter(msg => visibleSources[msg.source]).map((msg, idx) => (
                <div key={msg.id} className="bg-white/5 rounded-lg p-4 border border-purple-400/30 animate-slideIn opacity-0" style={{ animationDelay: `${idx * 50}ms`, animationFillMode: 'forwards' }}>
                  <div className="flex items-center justify-between mb-2">
                    <div className="text-xs text-purple-400">
                      {new Date(msg.timestamp).toLocaleTimeString()}
                    </div>
                    <div className="px-2 py-0.5 bg-purple-500/30 rounded text-xs font-semibold text-purple-200">
                      {msg.source}
                    </div>
                  </div>
                  <div className="text-sm text-purple-100">
                    <MarkdownText text={msg.message} />
                  </div>
                  {msg.smiles && (
                    <div className="mt-3 bg-white/50 rounded-lg p-2 flex justify-center">
                      <MoleculeSVG smiles={msg.smiles} height={120} rdkitModule={rdkitModule} />
                    </div>
                  )}
                </div>
              ))
            )}
          </div>
        </div>
      </div>
    );
};