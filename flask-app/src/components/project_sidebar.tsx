import React, { useState } from 'react';
import { ChevronRight, ChevronDown, ChevronLeft, Plus, FolderOpen, FileText, Loader2, Edit2, Trash2 } from 'lucide-react';
import { Project, Experiment, ProjectSelection } from '../types';
import { ProjectData } from '../hooks/useProjectData';

interface ProjectSidebarProps {
  projectData: ProjectData;
  isOpen: boolean;
  onToggle: () => void;
  selection: ProjectSelection;
  onSelectionChange: (selection: ProjectSelection) => void;
  onLoadContext: (projectId: string, experimentId: string | null) => void;
  onSaveContext: () => boolean;
  onReset: () => void;
  isComputing?: boolean;  // Track if main app is currently computing
  onCreateProjectAndExperiment?: (projectName: string, experimentName: string) => Promise<{ projectId: string; experimentId: string }>;
}

export const ProjectSidebar: React.FC<ProjectSidebarProps> = ({
  projectData,
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
    projectsRef,
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    setExperimentRunning
  } = projectData;

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

  const editingProjectRef = React.useRef<HTMLDivElement>(null);
  const editingExperimentRef = React.useRef<HTMLDivElement>(null);

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

  React.useEffect(() => {
    const handleClickOutside = (event: MouseEvent): void => {
      // Check if editing a project
      if (editingProject && editingProjectRef.current) {
        if (!editingProjectRef.current.contains(event.target as Node)) {
          setEditingProject(null);
          setEditProjectName('');
        }
      }

      // Check if editing an experiment
      if (editingExperiment && editingExperimentRef.current) {
        if (!editingExperimentRef.current.contains(event.target as Node)) {
          setEditingExperiment(null);
          setEditExperimentName('');
        }
      }
    };

    if (editingProject || editingExperiment) {
      window.addEventListener('mousedown', handleClickOutside);
      return () => window.removeEventListener('mousedown', handleClickOutside);
    }
  }, [editingProject, editingExperiment]);


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
    if (onSaveContext())
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
    if (onSaveContext())
      onReset();

    const newSelection: ProjectSelection = {
      projectId: project.id,
      experimentId: experiment.id
    };
    onSelectionChange(newSelection);
    onLoadContext(project.id, experiment.id);
  };

  const handleCreateProject = async () => {
    const now = new Date();
    const month = String(now.getMonth() + 1).padStart(2, '0');
    const day = String(now.getDate()).padStart(2, '0');
    const year = String(now.getFullYear()).slice(-2);
    const hours = String(now.getHours()).padStart(2, '0');
    const minutes = String(now.getMinutes()).padStart(2, '0');
    const timestamp = `${month}/${day}/${year} ${hours}:${minutes}`;

    const newProjectName = `Project ${timestamp}`;

    // Save before selecting away
    if (onSaveContext()) {
      // Reset UI state
      onReset();
    }

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
    const curProject = projectsRef.current!.find(p => p.id === projectId)!;

    let i = curProject.experiments.length + 1;
    let newExperimentName = `Experiment ${i}`
    while (curProject.experiments.findIndex(e => e.name == newExperimentName) !== -1) {
      i += 1;
      newExperimentName = `Experiment ${i}`;
    }

    // Save before selecting away
    if (onSaveContext()) {
      // Reset UI state
      onReset();
    }

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
    console.log('handleEditProject called with:', project);
    if (!project || !project.id) {
      console.error('Invalid project passed to handleEditProject');
      return;
    }
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
        <div className="sidebar sidebar-collapsed">
          <button
            onClick={onToggle}
            className="btn-icon"
            title="Open Projects"
          >
            <FolderOpen className="w-5 h-5" />
          </button>
        </div>

        {/* Delete Confirmation Modal (still needs to show when sidebar is collapsed) */}
        {deletingItem && (
          <div className="modal-overlay z-[60]">
            <div className="modal-content modal-content-sm">
              <div className="mb-4">
                <h3 className="modal-title mb-2">Confirm Delete</h3>
                <p className="text-secondary text-sm">
                  {deletingItem.type === 'project'
                    ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.experiments.length || 0} experiments inside will also be deleted.`
                    : 'Are you sure you want to delete this experiment?'}
                </p>
                <p className="text-tertiary text-xs mt-2">This action cannot be undone.</p>
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
                  className="btn btn-danger flex-1"
                >
                  Delete
                </button>
                <button
                  onClick={() => setDeletingItem(null)}
                  className="btn btn-tertiary flex-1"
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
      className={`sidebar ${isResizing ? 'resizing' : ''}`}
      style={{ width: `${sidebarWidth}px` }}
    >
      {/* Header */}
      <div className="sidebar-header">
        <div className="flex-between mb-2">
          <div className="flex items-center gap-2">
            <FolderOpen className="w-5 h-5 text-muted" />
            <h2 className="heading-3">Projects</h2>
          </div>
          <button
            onClick={onToggle}
            className="btn-icon"
            title="Collapse"
          >
            <ChevronLeft className="w-5 h-5" />
          </button>
        </div>

        {/* No Project/Experiment Selected Warning */}
        {(!selection.projectId || !selection.experimentId) && (
          <div className="alert alert-warning mt-2">
            <p className="text-warning text-xs">
              {!selection.projectId
                ? 'Warning: No project selected. A new project will be created when you run.'
                : 'Warning: No experiment selected. A new experiment will be created when you run.'}
            </p>
          </div>
        )}
      </div>

      {/* Content */}
      <div className="sidebar-content custom-scrollbar">
        {loading ? (
          <div className="flex-center py-8">
            <Loader2 className="w-6 h-6 text-muted animate-spin" />
          </div>
        ) : (
          <div className="space-y-1">
            {/* Projects List */}
            {projects.map((project) => (
              <div key={project.id} className="space-y-0.5">
                {/* Project Item */}
                <div className="project-item">
                  <button
                    onClick={() => toggleProjectExpanded(project.id)}
                    className="btn-icon p-1"
                  >
                    {expandedProjects.has(project.id) ? (
                      <ChevronDown className="w-4 h-4" />
                    ) : (
                      <ChevronRight className="w-4 h-4" />
                    )}
                  </button>

                  {editingProject === project.id ? (
                    <div ref={editingProjectRef} className="flex-1 flex items-center gap-1 px-2">
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
                        className="form-input flex-1 py-1 text-sm"
                        autoFocus
                      />
                      <button
                        onClick={() => handleSaveProjectEdit(project)}
                        className="btn btn-secondary btn-sm"
                      >
                        Save
                      </button>
                      <button
                        onClick={() => {
                          setEditingProject(null);
                          setEditProjectName('');
                        }}
                        className="btn btn-tertiary btn-sm"
                      >
                        Cancel
                      </button>
                    </div>
                  ) : (
                    <>
                      <button
                        disabled={isComputing}
                        onClick={() => handleProjectClick(project)}
                        className={`project-button ${
                          selection.projectId === project.id ? 'project-button-active' : ''
                        }`}
                      >
                        <FolderOpen className="w-4 h-4 flex-shrink-0" />
                        <span className="truncate flex-1">{project.name}</span>
                        {project.experiments.length > 0 && (
                          <span className="text-xs text-muted flex-shrink-0">
                            {project.experiments.length}
                          </span>
                        )}
                      </button>
                      {/* Edit/Delete buttons with background overlay */}
                      <div className="project-actions">
                        <div className="project-actions-bg">
                          <button
                            onClick={(e) => {
                              e.stopPropagation();
                              handleEditProject(project);
                            }}
                            className="action-button"
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
                            className="action-button action-button-danger"
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
                  <div className="experiment-item space-y-0.5">
                    {/* Experiment Items */}
                    {project.experiments.map((experiment) => (
                      <div key={experiment.id} className="project-item">
                        {editingExperiment?.experimentId === experiment.id ? (
                          <div ref={editingExperimentRef} className="flex-1 flex items-center gap-1 px-2 py-1">
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
                              className="form-input flex-1 py-1 text-xs"
                              autoFocus
                            />
                            <button
                              onClick={() => handleSaveExperimentEdit(project.id, experiment)}
                              className="btn btn-secondary btn-sm text-xs"
                            >
                              Save
                            </button>
                            <button
                              onClick={() => {
                                setEditingExperiment(null);
                                setEditExperimentName('');
                              }}
                              className="btn btn-tertiary btn-sm text-xs"
                            >
                              Cancel
                            </button>
                          </div>
                        ) : (
                          <>
                            <button
                              disabled={isComputing}
                              onClick={() => handleExperimentClick(project, experiment)}
                              className={`experiment-button ${
                                selection.experimentId === experiment.id ? 'experiment-button-active' : ''
                              }`}
                            >
                              {experiment.isRunning ? (
                                <Loader2 className="w-3 h-3 flex-shrink-0 animate-spin text-muted" />
                              ) : (
                                <FileText className="w-3 h-3 flex-shrink-0" />
                              )}
                              <span className="truncate flex-1">{experiment.name}</span>
                            </button>
                            {/* Edit/Delete buttons with background overlay */}
                            <div className="project-actions">
                              <div className="project-actions-bg">
                                <button
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    handleEditExperiment(project.id, experiment);
                                  }}
                                  className="action-button"
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
                                  className="action-button action-button-danger"
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
                        disabled={isComputing}
                        onClick={() => handleCreateExperiment(project.id)}
                        className="experiment-button group w-full"
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
                          className="form-input py-1 text-xs"
                          autoFocus
                        />
                        <div className="flex gap-2">
                          <button
                            disabled={isComputing}
                            onClick={() => handleCreateExperiment(project.id)}
                            className="btn btn-secondary btn-sm flex-1 text-xs"
                          >
                            Create
                          </button>
                          <button
                            onClick={() => {
                              setCreatingExperimentFor(null);
                              setNewExperimentName('');
                            }}
                            className="btn btn-tertiary btn-sm flex-1 text-xs"
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
                disabled={isComputing}
                onClick={handleCreateProject}
                className="project-button w-full group"
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
                  className="form-input py-1 text-sm"
                  autoFocus
                />
                <div className="flex gap-2">
                  <button
                    disabled={isComputing}
                    onClick={handleCreateProject}
                    className="btn btn-secondary btn-sm flex-1 text-xs"
                  >
                    Create
                  </button>
                  <button
                    onClick={() => {
                      setCreatingProject(false);
                      setNewProjectName('');
                    }}
                    className="btn btn-tertiary btn-sm flex-1 text-xs"
                  >
                    Cancel
                  </button>
                </div>
              </div>
            )}

            {projects.length === 0 && !creatingProject && (
              <div className="empty-state py-8">
                <p>No projects yet</p>
                <p className="text-xs mt-1">Click "New Project" to get started</p>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Resize Handle */}
      <div
        className={`sidebar-resize-handle ${isResizing ? 'bg-secondary' : ''}`}
        onMouseDown={(e) => {
          e.preventDefault();
          setIsResizing(true);
        }}
        title="Drag to resize (drag left to collapse)"
      >
        {/* Visual indicator line */}
        <div className="absolute top-0 bottom-0 right-0 w-0.5 bg-secondary group-hover:bg-primary transition-colors" />
        {/* Center grip indicator */}
        <div className="absolute top-1/2 right-0.5 -translate-y-1/2 flex flex-col gap-1">
          <div className="w-0.5 h-1 bg-secondary group-hover:bg-primary transition-colors" />
          <div className="w-0.5 h-1 bg-secondary group-hover:bg-primary transition-colors" />
          <div className="w-0.5 h-1 bg-secondary group-hover:bg-primary transition-colors" />
        </div>
      </div>
    </div>

    {/* Delete Confirmation Modal */}
    {deletingItem && (
      <div className="modal-overlay z-[60]">
        <div className="modal-content modal-content-sm">
          <div className="mb-4">
            <h3 className="modal-title mb-2">Confirm Delete</h3>
            <p className="text-secondary text-sm">
              {deletingItem.type === 'project'
                ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.experiments.length || 0} experiments inside will also be deleted.`
                : 'Are you sure you want to delete this experiment?'}
            </p>
            <p className="text-tertiary text-xs mt-2">This action cannot be undone.</p>
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
              className="btn btn-danger flex-1"
            >
              Delete
            </button>
            <button
              onClick={() => setDeletingItem(null)}
              className="btn btn-tertiary flex-1"
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

  // Initialize selection from localStorage so it survives browser restarts
  const [isOpen, setIsOpen] = useState(true);
  const [selection, setSelectionState] = useState<ProjectSelection>(() => {
    try {
      const stored = localStorage.getItem(SELECTION_STORAGE_KEY);
      if (stored) {
        const parsed = JSON.parse(stored);
        if (parsed && (parsed.projectId || parsed.experimentId)) {
          console.log('Restored sidebar selection from localStorage:', parsed);
          return parsed as ProjectSelection;
        }
      }
    } catch (e) {
      console.error('Error reading sidebar selection from localStorage:', e);
    }
    return { projectId: null, experimentId: null };
  });

  // Save selection to localStorage whenever it changes
  const setSelection = React.useCallback((newSelection: ProjectSelection) => {
    setSelectionState(newSelection);
    localStorage.setItem(SELECTION_STORAGE_KEY, JSON.stringify(newSelection));
  }, []);

  const selectionRef = React.useRef(selection);
  React.useEffect(() => {
    selectionRef.current = selection;
  }, [selection]);

  return {
    isOpen,
    setIsOpen,
    toggleSidebar: () => setIsOpen(prev => !prev),
    selection,
    setSelection,
    selectionRef
  };
};

// Separate hook for project management functions to be used in the main app
export const useProjectManagement = (projectData: ProjectData) => {
  const { projectsRef, projects, createProject, createExperiment, updateExperiment } = projectData;

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
    projectsRef,  // Expose projects list
    projects,
    createProjectAndExperiment,
    createExperiment,
    updateExperiment
  };
};
