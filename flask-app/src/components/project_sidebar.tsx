import React, { useState } from 'react';
import { ChevronRight, ChevronDown, ChevronLeft, Plus, FolderOpen, FileText, Loader2, Edit2, Trash2 } from 'lucide-react';
import { Project, Experiment, ProjectSelection } from '../types';
import { useProjectData } from '../hooks/useProjectData';

interface ProjectSidebarProps {
  isOpen: boolean;
  onToggle: () => void;
  selection: ProjectSelection;
  onSelectionChange: (selection: ProjectSelection) => void;
  onLoadContext: (projectId: string, experimentId: string | null) => void;
  onSaveContext: () => void;
  onReset: () => void;
  isComputing?: boolean;  // Track if main app is currently computing
  onCreateProjectAndExperiment?: (projectName: string, experimentName: string) => Promise<{ projectId: string; experimentId: string }>;
}

export const ProjectSidebar: React.FC<ProjectSidebarProps> = ({
  isOpen,
  onToggle,
  selection,
  onSelectionChange,
  onLoadContext,
  onSaveContext,
  onReset,
  isComputing = false,
  onCreateProjectAndExperiment
}) => {
  const {
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    setExperimentRunning
  } = useProjectData();

  const SIDEBAR_WIDTH_STORAGE_KEY = 'flask_copilot_sidebar_width';
  const MIN_WIDTH = 200;
  const MAX_WIDTH = 1600;
  const DEFAULT_WIDTH = 320;
  const COLLAPSE_THRESHOLD = 5;

  const [expandedProjects, setExpandedProjects] = useState<Set<string>>(new Set());
  const [creatingProject, setCreatingProject] = useState(false);
  const [newProjectName, setNewProjectName] = useState('');
  const [creatingExperimentFor, setCreatingExperimentFor] = useState<string | null>(null);
  const [newExperimentName, setNewExperimentName] = useState('');
  const [editingProject, setEditingProject] = useState<string | null>(null);
  const [editProjectName, setEditProjectName] = useState('');
  const [editingExperiment, setEditingExperiment] = useState<{ projectId: string; experimentId: string } | null>(null);
  const [editExperimentName, setEditExperimentName] = useState('');
  const [deletingItem, setDeletingItem] = useState<{ type: 'project' | 'experiment'; projectId: string; experimentId?: string } | null>(null);
  const [hasRestoredSelection, setHasRestoredSelection] = useState(false);
  const [sidebarWidth, setSidebarWidth] = useState<number>(() => {
    const saved = localStorage.getItem(SIDEBAR_WIDTH_STORAGE_KEY);
    return saved ? parseInt(saved, 10) : DEFAULT_WIDTH;
  });
  const [isResizing, setIsResizing] = useState(false);

  // Save width to localStorage when it changes
  React.useEffect(() => {
    localStorage.setItem(SIDEBAR_WIDTH_STORAGE_KEY, sidebarWidth.toString());
  }, [sidebarWidth]);

  // Handle resize
  React.useEffect(() => {
    if (!isResizing) return;

    // Add cursor style to body
    document.body.style.cursor = 'col-resize';
    document.body.style.userSelect = 'none';

    const handleMouseMove = (e: MouseEvent) => {
      const newWidth = e.clientX;
      
      // Check if width falls below collapse threshold
      if (newWidth < COLLAPSE_THRESHOLD) {
        onToggle(); // Collapse the sidebar
        setIsResizing(false);
        return;
      }
      
      // Constrain width between min and max
      if (newWidth >= MIN_WIDTH && newWidth <= MAX_WIDTH) {
        setSidebarWidth(newWidth);
      }
    };

    const handleMouseUp = () => {
      setIsResizing(false);
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
    };

    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
    };
  }, [isResizing, onToggle]);

  // Update experiment running status when isComputing changes
  React.useEffect(() => {
    if (selection.projectId && selection.experimentId) {
      setExperimentRunning(selection.projectId, selection.experimentId, isComputing);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isComputing, selection.projectId, selection.experimentId]);

  const toggleProjectExpanded = (projectId: string) => {
    const newExpanded = new Set(expandedProjects);
    if (newExpanded.has(projectId)) {
      newExpanded.delete(projectId);
    } else {
      newExpanded.add(projectId);
    }
    setExpandedProjects(newExpanded);
  };

  const handleProjectClick = (project: Project) => {
    // No change
    if (project.id == selection.projectId) {
      return;
    }

    // Save before selecting away
    onSaveContext();
    onReset();

    // Auto-select the last experiment if the project has any
    const lastExperiment = project.experiments.length > 0 
      ? project.experiments[project.experiments.length - 1] 
      : null;
    
    const newSelection: ProjectSelection = {
      projectId: project.id,
      experimentId: lastExperiment?.id || null
    };
    
    onSelectionChange(newSelection);
    onLoadContext(project.id, lastExperiment?.id || null);
    
    // Auto-expand when clicking a project
    if (!expandedProjects.has(project.id)) {
      toggleProjectExpanded(project.id);
    }
  };

  const handleExperimentClick = (project: Project, experiment: Experiment) => {
    // No change
    if (project.id == selection.projectId && experiment.id == selection.experimentId) {
      return;
    }

    // Save before selecting away
    onSaveContext();
    onReset();
    
    const newSelection: ProjectSelection = {
      projectId: project.id,
      experimentId: experiment.id
    };
    onSelectionChange(newSelection);
    onLoadContext(project.id, experiment.id);
  };

  const handleCreateProject = async () => {
    if (!newProjectName.trim()) return;
    
    // Save before selecting away
    onSaveContext();

    // Reset UI state
    onReset();

    try {
      const project = await createProject(newProjectName);
      setNewProjectName('');
      setCreatingProject(false);
      // Auto-select and expand the new project
      onSelectionChange({ projectId: project.id, experimentId: null });
      setExpandedProjects(prev => new Set(prev).add(project.id));
    } catch (error) {
      console.error('Error creating project:', error);
      alert('Failed to create project');
    }
  };

  const handleCreateExperiment = async (projectId: string) => {
    if (!newExperimentName.trim()) return;

    // Save before selecting away
    onSaveContext();

    // Reset UI state
    onReset();
    
    try {
      const experiment = await createExperiment(projectId, newExperimentName);
      setNewExperimentName('');
      setCreatingExperimentFor(null);
      // Auto-select the new experiment
      onSelectionChange({ projectId, experimentId: experiment.id });
      onLoadContext(projectId, experiment.id);
    } catch (error) {
      console.error('Error creating experiment:', error);
      alert('Failed to create experiment');
    }
  };

  const handleEditProject = (project: Project) => {
    setEditingProject(project.id);
    setEditProjectName(project.name);
  };

  const handleSaveProjectEdit = async (project: Project) => {
    if (!editProjectName.trim()) return;
    
    try {
      await updateProject({ ...project, name: editProjectName });
      setEditingProject(null);
      setEditProjectName('');
    } catch (error) {
      console.error('Error updating project:', error);
      alert('Failed to update project');
    }
  };

  const handleDeleteProject = async (projectId: string) => {
    try {
      await deleteProject(projectId);
      setDeletingItem(null);
      // Clear selection if deleted project was selected
      if (selection.projectId === projectId) {
        onSelectionChange({ projectId: null, experimentId: null });
      }
    } catch (error) {
      console.error('Error deleting project:', error);
      alert('Failed to delete project');
    }
  };

  const handleEditExperiment = (projectId: string, experiment: Experiment) => {
    setEditingExperiment({ projectId, experimentId: experiment.id });
    setEditExperimentName(experiment.name);
  };

  const handleSaveExperimentEdit = async (projectId: string, experiment: Experiment) => {
    if (!editExperimentName.trim()) return;
    
    try {
      await updateExperiment(projectId, { ...experiment, name: editExperimentName });
      setEditingExperiment(null);
      setEditExperimentName('');
    } catch (error) {
      console.error('Error updating experiment:', error);
      alert('Failed to update experiment');
    }
  };

  const handleDeleteExperiment = async (projectId: string, experimentId: string) => {
    try {
      await deleteExperiment(projectId, experimentId);
      setDeletingItem(null);
      // Clear selection if deleted experiment was selected
      if (selection.experimentId === experimentId) {
        onSelectionChange({ projectId: selection.projectId, experimentId: null });
      }
    } catch (error) {
      console.error('Error deleting experiment:', error);
      alert('Failed to delete experiment');
    }
  };

  if (!isOpen) {
    return (
      <>
        <div className="w-12 bg-slate-800/50 backdrop-blur-sm border-r border-purple-400/30 flex flex-col items-center py-4">
          <button
            onClick={onToggle}
            className="p-2 text-purple-300 hover:text-white hover:bg-purple-600/30 rounded-lg transition-all"
            title="Open Projects"
          >
            <FolderOpen className="w-5 h-5" />
          </button>
        </div>

        {/* Delete Confirmation Modal (still needs to show when sidebar is collapsed) */}
        {deletingItem && (
          <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-[60] flex items-center justify-center p-4">
            <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-md w-full p-6">
              <div className="mb-4">
                <h3 className="text-lg font-bold text-white mb-2">Confirm Delete</h3>
                <p className="text-purple-200 text-sm">
                  {deletingItem.type === 'project' 
                    ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.experiments.length || 0} experiments inside will also be deleted.`
                    : 'Are you sure you want to delete this experiment?'}
                </p>
                <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
              </div>
              
              <div className="flex gap-3">
                <button
                  onClick={() => {
                    if (deletingItem.type === 'project') {
                      handleDeleteProject(deletingItem.projectId);
                    } else {
                      handleDeleteExperiment(deletingItem.projectId, deletingItem.experimentId!);
                    }
                  }}
                  className="flex-1 px-4 py-2 bg-red-600 hover:bg-red-500 text-white rounded-lg font-semibold transition-all"
                >
                  Delete
                </button>
                <button
                  onClick={() => setDeletingItem(null)}
                  className="flex-1 px-4 py-2 bg-white/20 hover:bg-white/30 text-white rounded-lg font-semibold transition-all"
                >
                  Cancel
                </button>
              </div>
            </div>
          </div>
        )}
      </>
    );
  }

  return (
    <>
    <div 
      className="bg-slate-800/50 backdrop-blur-sm border-r border-purple-400/30 flex flex-col relative"
      style={{ width: `${sidebarWidth}px` }}
    >
      {/* Header */}
      <div className="p-4 border-b border-purple-400/30">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-2">
            <FolderOpen className="w-5 h-5 text-purple-400" />
            <h2 className="text-lg font-semibold text-white">Projects</h2>
          </div>
          <button
            onClick={onToggle}
            className="p-1.5 text-purple-300 hover:text-white hover:bg-purple-600/30 rounded-lg transition-all"
            title="Collapse"
          >
            <ChevronLeft className="w-5 h-5" />
          </button>
        </div>
        
        {/* No Project/Experiment Selected Warning */}
        {(!selection.projectId || !selection.experimentId) && (
          <div className="mt-2 px-3 py-2 bg-amber-500/20 border border-amber-400/50 rounded-lg">
            <p className="text-amber-200 text-xs">
              {!selection.projectId 
                ? '⚠️ No project selected. A new project will be created when you run.'
                : '⚠️ No experiment selected. A new experiment will be created when you run.'}
            </p>
          </div>
        )}
      </div>

      {/* Content */}
      <div className="flex-1 overflow-y-auto p-4">
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <Loader2 className="w-6 h-6 text-purple-400 animate-spin" />
          </div>
        ) : (
          <div className="space-y-1">
            {/* Projects List */}
            {projects.map((project) => (
              <div key={project.id} className="space-y-0.5">
                {/* Project Item */}
                <div className="flex items-center group relative">
                  <button
                    onClick={() => toggleProjectExpanded(project.id)}
                    className="p-1 text-purple-300 hover:text-white transition-colors flex-shrink-0"
                  >
                    {expandedProjects.has(project.id) ? (
                      <ChevronDown className="w-4 h-4" />
                    ) : (
                      <ChevronRight className="w-4 h-4" />
                    )}
                  </button>
                  
                  {editingProject === project.id ? (
                    <div className="flex-1 flex items-center gap-1 px-2">
                      <input
                        type="text"
                        value={editProjectName}
                        onChange={(e) => setEditProjectName(e.target.value)}
                        onKeyDown={(e) => {
                          if (e.key === 'Enter') handleSaveProjectEdit(project);
                          if (e.key === 'Escape') {
                            setEditingProject(null);
                            setEditProjectName('');
                          }
                        }}
                        className="flex-1 px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-sm focus:outline-none focus:border-purple-400"
                        autoFocus
                      />
                      <button
                        onClick={() => handleSaveProjectEdit(project)}
                        className="px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                      >
                        Save
                      </button>
                      <button
                        onClick={() => {
                          setEditingProject(null);
                          setEditProjectName('');
                        }}
                        className="px-2 py-1 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                      >
                        Cancel
                      </button>
                    </div>
                  ) : (
                    <>
                      <button
                        disabled={isComputing}
                        onClick={() => handleProjectClick(project)}
                        className={`flex-1 px-3 py-2 text-left text-sm rounded-lg transition-all flex items-center gap-2 min-w-0 disabled:cursor-not-allowed disabled:hover:bg-white/30 ${
                          selection.projectId === project.id
                            ? 'bg-purple-600/50 text-white font-medium'
                            : 'text-purple-200 hover:text-white hover:bg-purple-600/30'
                        }`}
                      >
                        <FolderOpen className="w-4 h-4 flex-shrink-0" />
                        <span className="truncate flex-1">{project.name}</span>
                        {project.experiments.length > 0 && (
                          <span className="text-xs text-purple-400 flex-shrink-0">
                            {project.experiments.length}
                          </span>
                        )}
                      </button>
                      {/* Edit/Delete buttons with background overlay */}
                      <div className="absolute right-1 flex items-center gap-0.5 opacity-0 group-hover:opacity-100 transition-opacity">
                        <div className="bg-gradient-to-l from-slate-800 via-slate-800 to-transparent pl-8 pr-1 py-1 flex items-center gap-0.5">
                          <button
                            onClick={(e) => {
                              e.stopPropagation();
                              handleEditProject(project);
                            }}
                            className="p-1 text-purple-300 hover:text-white hover:bg-purple-600/50 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                            title="Edit project"
                          >
                            <Edit2 className="w-3.5 h-3.5" />
                          </button>
                          <button
                            disabled={isComputing}
                            onClick={(e) => {
                              e.stopPropagation();
                              setDeletingItem({ type: 'project', projectId: project.id });
                            }}
                            className="p-1 text-purple-300 hover:text-red-400 hover:bg-red-500/30 rounded transition-all backdrop-blur-sm bg-slate-800/80 disabled:opacity-50 disabled:cursor-not-allowed"
                            title="Delete project"
                          >
                            <Trash2 className="w-3.5 h-3.5" />
                          </button>
                        </div>
                      </div>
                    </>
                  )}
                </div>

                {/* Experiments List */}
                {expandedProjects.has(project.id) && (
                  <div className="ml-6 space-y-0.5">
                    {/* Experiment Items */}
                    {project.experiments.map((experiment) => (
                      <div key={experiment.id} className="flex items-center group relative">
                        {editingExperiment?.experimentId === experiment.id ? (
                          <div className="flex-1 flex items-center gap-1 px-2 py-1">
                            <input
                              type="text"
                              value={editExperimentName}
                              onChange={(e) => setEditExperimentName(e.target.value)}
                              onKeyDown={(e) => {
                                if (e.key === 'Enter') handleSaveExperimentEdit(project.id, experiment);
                                if (e.key === 'Escape') {
                                  setEditingExperiment(null);
                                  setEditExperimentName('');
                                }
                              }}
                              className="flex-1 px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-xs focus:outline-none focus:border-purple-400"
                              autoFocus
                            />
                            <button
                              onClick={() => handleSaveExperimentEdit(project.id, experiment)}
                              className="px-2 py-0.5 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                            >
                              Save
                            </button>
                            <button
                              onClick={() => {
                                setEditingExperiment(null);
                                setEditExperimentName('');
                              }}
                              className="px-2 py-0.5 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                            >
                              Cancel
                            </button>
                          </div>
                        ) : (
                          <>
                            <button
                              disabled={isComputing}
                              onClick={() => handleExperimentClick(project, experiment)}
                              className={`flex-1 px-3 py-1.5 text-left text-xs rounded transition-all flex items-center gap-2 min-w-0 ${
                                selection.experimentId === experiment.id
                                  ? 'bg-purple-600/50 text-white font-medium'
                                  : 'text-purple-200 hover:text-white hover:bg-purple-600/20 disabled:cursor-not-allowed disabled:hover:bg-white/0'
                              }`}
                            >
                              {experiment.isRunning ? (
                                <Loader2 className="w-3 h-3 flex-shrink-0 animate-spin text-purple-400" />
                              ) : (
                                <FileText className="w-3 h-3 flex-shrink-0" />
                              )}
                              <span className="truncate flex-1">{experiment.name}</span>
                            </button>
                            {/* Edit/Delete buttons with background overlay */}
                            <div className="absolute right-1 flex items-center gap-0.5 opacity-0 group-hover:opacity-100 transition-opacity">
                              <div className="bg-gradient-to-l from-slate-800 via-slate-800 to-transparent pl-6 pr-1 py-0.5 flex items-center gap-0.5">
                                <button
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    handleEditExperiment(project.id, experiment);
                                  }}
                                  className="p-0.5 text-purple-300 hover:text-white hover:bg-purple-600/50 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                                  title="Edit experiment"
                                >
                                  <Edit2 className="w-3 h-3" />
                                </button>
                                <button
                                  disabled={isComputing}
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    setDeletingItem({ type: 'experiment', projectId: project.id, experimentId: experiment.id });
                                  }}
                                  className="p-0.5 text-purple-300 hover:text-red-400 hover:bg-red-500/30 rounded transition-all backdrop-blur-sm bg-slate-800/80 disabled:opacity-50 disabled:cursor-not-allowed"
                                  title="Delete experiment"
                                >
                                  <Trash2 className="w-3 h-3" />
                                </button>
                              </div>
                            </div>
                          </>
                        )}
                      </div>
                    ))}

                    {/* Create New Experiment (at bottom) */}
                    {creatingExperimentFor !== project.id ? (
                      <button
                        onClick={() => setCreatingExperimentFor(project.id)}
                        className="w-full px-3 py-1.5 text-left text-xs text-purple-300 hover:text-white hover:bg-purple-600/20 rounded transition-all flex items-center gap-2 group"
                      >
                        <Plus className="w-3 h-3 group-hover:scale-110 transition-transform" />
                        <span>New Experiment</span>
                      </button>
                    ) : (
                      <div className="px-3 py-1.5 space-y-2">
                        <input
                          type="text"
                          value={newExperimentName}
                          onChange={(e) => setNewExperimentName(e.target.value)}
                          onKeyDown={(e) => {
                            if (e.key === 'Enter' && !isComputing) handleCreateExperiment(project.id);
                            if (e.key === 'Escape') {
                              setCreatingExperimentFor(null);
                              setNewExperimentName('');
                            }
                          }}
                          placeholder="Experiment name..."
                          className="w-full px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-xs focus:outline-none focus:border-purple-400"
                          autoFocus
                        />
                        <div className="flex gap-2">
                          <button
                            disabled={isComputing}
                            onClick={() => handleCreateExperiment(project.id)}
                            className="flex-1 px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
                          >
                            Create
                          </button>
                          <button
                            onClick={() => {
                              setCreatingExperimentFor(null);
                              setNewExperimentName('');
                            }}
                            className="flex-1 px-2 py-1 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                          >
                            Cancel
                          </button>
                        </div>
                      </div>
                    )}
                  </div>
                )}
              </div>
            ))}

            {/* Create New Project Button (at bottom) */}
            {!creatingProject ? (
              <button
                onClick={() => setCreatingProject(true)}
                className="w-full px-3 py-2 text-left text-sm text-purple-300 hover:text-white hover:bg-purple-600/30 rounded-lg transition-all flex items-center gap-2 group"
              >
                <Plus className="w-4 h-4 group-hover:scale-110 transition-transform" />
                <span>New Project</span>
              </button>
            ) : (
              <div className="px-3 py-2 space-y-2">
                <input
                  type="text"
                  value={newProjectName}
                  onChange={(e) => setNewProjectName(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === 'Enter' && !isComputing) handleCreateProject();
                    if (e.key === 'Escape') {
                      setCreatingProject(false);
                      setNewProjectName('');
                    }
                  }}
                  placeholder="Project name..."
                  className="w-full px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-sm focus:outline-none focus:border-purple-400"
                  autoFocus
                />
                <div className="flex gap-2">
                  <button
                    disabled={isComputing}
                    onClick={handleCreateProject}
                    className="flex-1 px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
                  >
                    Create
                  </button>
                  <button
                    onClick={() => {
                      setCreatingProject(false);
                      setNewProjectName('');
                    }}
                    className="flex-1 px-2 py-1 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                  >
                    Cancel
                  </button>
                </div>
              </div>
            )}

            {projects.length === 0 && !creatingProject && (
              <div className="text-center py-8 text-purple-300 text-sm">
                <p>No projects yet</p>
                <p className="text-xs mt-1">Click "New Project" to get started</p>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Resize Handle */}
      <div
        className={`absolute top-0 right-0 bottom-0 w-2 cursor-col-resize group ${
          isResizing ? 'bg-purple-400/30' : ''
        }`}
        onMouseDown={(e) => {
          e.preventDefault();
          setIsResizing(true);
        }}
        title="Drag to resize (drag left to collapse)"
      >
        {/* Visual indicator line */}
        <div className="absolute top-0 bottom-0 right-0 w-0.5 bg-purple-400/20 group-hover:bg-purple-400/50 transition-colors" />
        {/* Center grip indicator */}
        <div className="absolute top-1/2 right-0.5 -translate-y-1/2 flex flex-col gap-1">
          <div className="w-0.5 h-1 bg-purple-400/40 group-hover:bg-purple-400/70 transition-colors" />
          <div className="w-0.5 h-1 bg-purple-400/40 group-hover:bg-purple-400/70 transition-colors" />
          <div className="w-0.5 h-1 bg-purple-400/40 group-hover:bg-purple-400/70 transition-colors" />
        </div>
      </div>
    </div>

    {/* Delete Confirmation Modal */}
    {deletingItem && (
      <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-[60] flex items-center justify-center p-4">
        <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-md w-full p-6">
          <div className="mb-4">
            <h3 className="text-lg font-bold text-white mb-2">Confirm Delete</h3>
            <p className="text-purple-200 text-sm">
              {deletingItem.type === 'project' 
                ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.experiments.length || 0} experiments inside will also be deleted.`
                : 'Are you sure you want to delete this experiment?'}
            </p>
            <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
          </div>
          
          <div className="flex gap-3">
            <button
              onClick={() => {
                if (deletingItem.type === 'project') {
                  handleDeleteProject(deletingItem.projectId);
                } else {
                  handleDeleteExperiment(deletingItem.projectId, deletingItem.experimentId!);
                }
              }}
              className="flex-1 px-4 py-2 bg-red-600 hover:bg-red-500 text-white rounded-lg font-semibold transition-all"
            >
              Delete
            </button>
            <button
              onClick={() => setDeletingItem(null)}
              className="flex-1 px-4 py-2 bg-white/20 hover:bg-white/30 text-white rounded-lg font-semibold transition-all"
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

// Custom hook to manage sidebar state
export const useProjectSidebar = () => {
  const SELECTION_STORAGE_KEY = 'flask_copilot_last_selection';
  
  const [isOpen, setIsOpen] = useState(true);
  const [selection, setSelectionState] = useState<ProjectSelection>({
    projectId: null,
    experimentId: null
  });

  // Save selection to localStorage whenever it changes
  const setSelection = React.useCallback((newSelection: ProjectSelection) => {
    setSelectionState(newSelection);
    localStorage.setItem(SELECTION_STORAGE_KEY, JSON.stringify(newSelection));
  }, []);

  return {
    isOpen,
    setIsOpen,
    toggleSidebar: () => setIsOpen(prev => !prev),
    selection,
    setSelection,
  };
};

// Separate hook for project management functions to be used in the main app
export const useProjectManagement = () => {
  const { projects, createProject, createExperiment, updateExperiment } = useProjectData();

  const createProjectAndExperiment = React.useCallback(async (
    projectName: string, 
    experimentName: string
  ): Promise<{ projectId: string; experimentId: string }> => {
    const project = await createProject(projectName);
    const experiment = await createExperiment(project.id, experimentName);
    
    return {
      projectId: project.id,
      experimentId: experiment.id
    };
  }, [createProject, createExperiment]);

  return {
    projects,  // Expose projects list
    createProjectAndExperiment,
    createExperiment,
    updateExperiment
  };
};
