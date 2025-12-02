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
      <div className="reasoning-sidebar flex-col">
        <div className="card-header">
          <h3 className="heading-3">Reasoning</h3>
          <div className="flex items-center gap-2">
            <div className="filter-control">
              <button 
                onClick={(e) => { e.stopPropagation(); setSourceFilterOpen(!sourceFilterOpen); }}
                onMouseDown={(e) => e.stopPropagation()}
                className="btn btn-secondary btn-sm"
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
                <div className="filter-menu" onClick={(e) => e.stopPropagation()} onMouseDown={(e) => e.stopPropagation()}>
                  <div className="filter-menu-header">
                    Message Sources
                  </div>
                  {Object.keys(visibleSources).map(source => (
                    <label key={source} className="filter-menu-item">
                      <input
                        type="checkbox"
                        checked={visibleSources[source]}
                        onChange={() => setVisibleSources(prev => ({ ...prev, [source]: !prev[source] }))}
                        className="form-checkbox"
                      />
                      <span className="text-sm text-primary">{source}</span>
                      <span className="ml-auto text-xs text-muted">
                        ({messages.filter(m => m.source === source).length})
                      </span>
                    </label>
                  ))}
                  <div className="filter-menu-footer">
                    <button
                      onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: true }), {}))}
                      className="btn btn-secondary btn-sm btn-filter flex-1"
                    >
                      All
                    </button>
                    <button
                      onClick={() => setVisibleSources(Object.keys(visibleSources).reduce((acc, key) => ({ ...acc, [key]: false }), {}))}
                      className="btn btn-secondary btn-sm btn-filter flex-1"
                    >
                      None
                    </button>
                  </div>
                </div>
              )}
            </div>
            <button onClick={() => setSidebarOpen(false)} className="btn-icon">
              <X className="w-5 h-5" />
            </button>
          </div>
        </div>

        <div className="reasoning-messages space-y-3 custom-scrollbar" ref={sidebarRef}>
          {messages.filter(msg => visibleSources[msg.source]).length === 0 ? (
            <div className="empty-state">
              <svg className="empty-state-icon" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z" />
              </svg>
              <p className="text-sm">
                {messages.length === 0 ? 'No messages yet' : 'No messages match the selected filters'}
              </p>
            </div>
          ) : (
            messages.filter(msg => visibleSources[msg.source]).map((msg, idx) => (
              <div key={`${msg.id}-${idx}`} className="message-card" style={{ animationDelay: `${idx * 50}ms` }}>
                <div className="flex items-center justify-between mb-2">
                  <div className="text-xs text-muted">
                    {new Date(msg.timestamp).toLocaleTimeString()}
                  </div>
                  <div className="badge badge-primary">
                    {msg.source}
                  </div>
                </div>
                <div className="text-sm text-secondary">
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
    );
};
