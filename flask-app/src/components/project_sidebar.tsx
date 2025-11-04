import React, { useState } from 'react';
import { ChevronRight, ChevronDown, ChevronLeft, Plus, FolderOpen, FileText, Loader2, Edit2, Trash2 } from 'lucide-react';
import { Project, Task, ProjectSelection } from '../types';
import { useProjectData } from '../hooks/useProjectData';

interface ProjectSidebarProps {
  isOpen: boolean;
  onToggle: () => void;
  selection: ProjectSelection;
  onSelectionChange: (selection: ProjectSelection) => void;
  onLoadContext: (projectId: string, taskId: string | null) => void;
  isComputing?: boolean;  // Track if main app is currently computing
  hasLoadedInitialSelection?: boolean;  // Track if initial selection from localStorage has been loaded
  onCreateProjectAndTask?: (projectName: string, taskName: string) => Promise<{ projectId: string; taskId: string }>;
}

export const ProjectSidebar: React.FC<ProjectSidebarProps> = ({
  isOpen,
  onToggle,
  selection,
  onSelectionChange,
  onLoadContext,
  isComputing = false,
  hasLoadedInitialSelection = false,
  onCreateProjectAndTask
}) => {
  const {
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createTask,
    updateTask,
    deleteTask,
    setTaskRunning
  } = useProjectData();

  const [expandedProjects, setExpandedProjects] = useState<Set<string>>(new Set());
  const [creatingProject, setCreatingProject] = useState(false);
  const [newProjectName, setNewProjectName] = useState('');
  const [creatingTaskFor, setCreatingTaskFor] = useState<string | null>(null);
  const [newTaskName, setNewTaskName] = useState('');
  const [editingProject, setEditingProject] = useState<string | null>(null);
  const [editProjectName, setEditProjectName] = useState('');
  const [editingTask, setEditingTask] = useState<{ projectId: string; taskId: string } | null>(null);
  const [editTaskName, setEditTaskName] = useState('');
  const [deletingItem, setDeletingItem] = useState<{ type: 'project' | 'task'; projectId: string; taskId?: string } | null>(null);
  const [hasRestoredSelection, setHasRestoredSelection] = useState(false);

  // Auto-expand and load context when initial selection is loaded from localStorage
  React.useEffect(() => {
    if (hasLoadedInitialSelection && !loading && !hasRestoredSelection && projects.length > 0 && selection.projectId) {
      // Verify the selection still exists in the projects
      const project = projects.find(p => p.id === selection.projectId);
      if (project) {
        // Auto-expand the project
        setExpandedProjects(prev => new Set(prev).add(selection.projectId!));
        
        // Load the context
        onLoadContext(selection.projectId, selection.taskId);
        setHasRestoredSelection(true);
      }
    }
  }, [hasLoadedInitialSelection, loading, projects, selection, onLoadContext, hasRestoredSelection]);

  // Update task running status when isComputing changes
  React.useEffect(() => {
    if (selection.projectId && selection.taskId) {
      setTaskRunning(selection.projectId, selection.taskId, isComputing);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isComputing, selection.projectId, selection.taskId]);

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
    // Auto-select the last task if the project has any
    const lastTask = project.tasks.length > 0 
      ? project.tasks[project.tasks.length - 1] 
      : null;
    
    const newSelection: ProjectSelection = {
      projectId: project.id,
      taskId: lastTask?.id || null
    };
    
    onSelectionChange(newSelection);
    onLoadContext(project.id, lastTask?.id || null);
    
    // Auto-expand when clicking a project
    if (!expandedProjects.has(project.id)) {
      toggleProjectExpanded(project.id);
    }
  };

  const handleTaskClick = (project: Project, task: Task) => {
    const newSelection: ProjectSelection = {
      projectId: project.id,
      taskId: task.id
    };
    onSelectionChange(newSelection);
    onLoadContext(project.id, task.id);
  };

  const handleCreateProject = async () => {
    if (!newProjectName.trim()) return;
    
    try {
      const project = await createProject(newProjectName);
      setNewProjectName('');
      setCreatingProject(false);
      // Auto-select and expand the new project
      onSelectionChange({ projectId: project.id, taskId: null });
      setExpandedProjects(prev => new Set(prev).add(project.id));
    } catch (error) {
      console.error('Error creating project:', error);
      alert('Failed to create project');
    }
  };

  const handleCreateTask = async (projectId: string) => {
    if (!newTaskName.trim()) return;
    
    try {
      const task = await createTask(projectId, newTaskName);
      setNewTaskName('');
      setCreatingTaskFor(null);
      // Auto-select the new task
      onSelectionChange({ projectId, taskId: task.id });
      onLoadContext(projectId, task.id);
    } catch (error) {
      console.error('Error creating task:', error);
      alert('Failed to create task');
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
        onSelectionChange({ projectId: null, taskId: null });
      }
    } catch (error) {
      console.error('Error deleting project:', error);
      alert('Failed to delete project');
    }
  };

  const handleEditTask = (projectId: string, task: Task) => {
    setEditingTask({ projectId, taskId: task.id });
    setEditTaskName(task.name);
  };

  const handleSaveTaskEdit = async (projectId: string, task: Task) => {
    if (!editTaskName.trim()) return;
    
    try {
      await updateTask(projectId, { ...task, name: editTaskName });
      setEditingTask(null);
      setEditTaskName('');
    } catch (error) {
      console.error('Error updating task:', error);
      alert('Failed to update task');
    }
  };

  const handleDeleteTask = async (projectId: string, taskId: string) => {
    try {
      await deleteTask(projectId, taskId);
      setDeletingItem(null);
      // Clear selection if deleted task was selected
      if (selection.taskId === taskId) {
        onSelectionChange({ projectId: selection.projectId, taskId: null });
      }
    } catch (error) {
      console.error('Error deleting task:', error);
      alert('Failed to delete task');
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
                    ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.tasks.length || 0} tasks inside will also be deleted.`
                    : 'Are you sure you want to delete this task?'}
                </p>
                <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
              </div>
              
              <div className="flex gap-3">
                <button
                  onClick={() => {
                    if (deletingItem.type === 'project') {
                      handleDeleteProject(deletingItem.projectId);
                    } else {
                      handleDeleteTask(deletingItem.projectId, deletingItem.taskId!);
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
    <div className="w-80 bg-slate-800/50 backdrop-blur-sm border-r border-purple-400/30 flex flex-col">
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
        
        {/* No Project/Task Selected Warning */}
        {(!selection.projectId || !selection.taskId) && (
          <div className="mt-2 px-3 py-2 bg-amber-500/20 border border-amber-400/50 rounded-lg">
            <p className="text-amber-200 text-xs">
              {!selection.projectId 
                ? '⚠️ No project selected. A new project will be created when you run.'
                : '⚠️ No task selected. A new task will be created when you run.'}
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
                        onClick={() => handleProjectClick(project)}
                        className={`flex-1 px-3 py-2 text-left text-sm rounded-lg transition-all flex items-center gap-2 min-w-0 ${
                          selection.projectId === project.id
                            ? 'bg-purple-600/50 text-white font-medium'
                            : 'text-purple-200 hover:text-white hover:bg-purple-600/30'
                        }`}
                      >
                        <FolderOpen className="w-4 h-4 flex-shrink-0" />
                        <span className="truncate flex-1">{project.name}</span>
                        {project.tasks.length > 0 && (
                          <span className="text-xs text-purple-400 flex-shrink-0">
                            {project.tasks.length}
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
                            onClick={(e) => {
                              e.stopPropagation();
                              setDeletingItem({ type: 'project', projectId: project.id });
                            }}
                            className="p-1 text-purple-300 hover:text-red-400 hover:bg-red-500/30 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                            title="Delete project"
                          >
                            <Trash2 className="w-3.5 h-3.5" />
                          </button>
                        </div>
                      </div>
                    </>
                  )}
                </div>

                {/* Tasks List */}
                {expandedProjects.has(project.id) && (
                  <div className="ml-6 space-y-0.5">
                    {/* Task Items */}
                    {project.tasks.map((task) => (
                      <div key={task.id} className="flex items-center group relative">
                        {editingTask?.taskId === task.id ? (
                          <div className="flex-1 flex items-center gap-1 px-2 py-1">
                            <input
                              type="text"
                              value={editTaskName}
                              onChange={(e) => setEditTaskName(e.target.value)}
                              onKeyDown={(e) => {
                                if (e.key === 'Enter') handleSaveTaskEdit(project.id, task);
                                if (e.key === 'Escape') {
                                  setEditingTask(null);
                                  setEditTaskName('');
                                }
                              }}
                              className="flex-1 px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-xs focus:outline-none focus:border-purple-400"
                              autoFocus
                            />
                            <button
                              onClick={() => handleSaveTaskEdit(project.id, task)}
                              className="px-2 py-0.5 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                            >
                              Save
                            </button>
                            <button
                              onClick={() => {
                                setEditingTask(null);
                                setEditTaskName('');
                              }}
                              className="px-2 py-0.5 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                            >
                              Cancel
                            </button>
                          </div>
                        ) : (
                          <>
                            <button
                              onClick={() => handleTaskClick(project, task)}
                              className={`flex-1 px-3 py-1.5 text-left text-xs rounded transition-all flex items-center gap-2 min-w-0 ${
                                selection.taskId === task.id
                                  ? 'bg-purple-600/50 text-white font-medium'
                                  : 'text-purple-200 hover:text-white hover:bg-purple-600/20'
                              }`}
                            >
                              {task.isRunning ? (
                                <Loader2 className="w-3 h-3 flex-shrink-0 animate-spin text-purple-400" />
                              ) : (
                                <FileText className="w-3 h-3 flex-shrink-0" />
                              )}
                              <span className="truncate flex-1">{task.name}</span>
                            </button>
                            {/* Edit/Delete buttons with background overlay */}
                            <div className="absolute right-1 flex items-center gap-0.5 opacity-0 group-hover:opacity-100 transition-opacity">
                              <div className="bg-gradient-to-l from-slate-800 via-slate-800 to-transparent pl-6 pr-1 py-0.5 flex items-center gap-0.5">
                                <button
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    handleEditTask(project.id, task);
                                  }}
                                  className="p-0.5 text-purple-300 hover:text-white hover:bg-purple-600/50 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                                  title="Edit task"
                                >
                                  <Edit2 className="w-3 h-3" />
                                </button>
                                <button
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    setDeletingItem({ type: 'task', projectId: project.id, taskId: task.id });
                                  }}
                                  className="p-0.5 text-purple-300 hover:text-red-400 hover:bg-red-500/30 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                                  title="Delete task"
                                >
                                  <Trash2 className="w-3 h-3" />
                                </button>
                              </div>
                            </div>
                          </>
                        )}
                      </div>
                    ))}

                    {/* Create New Task (at bottom) */}
                    {creatingTaskFor !== project.id ? (
                      <button
                        onClick={() => setCreatingTaskFor(project.id)}
                        className="w-full px-3 py-1.5 text-left text-xs text-purple-300 hover:text-white hover:bg-purple-600/20 rounded transition-all flex items-center gap-2 group"
                      >
                        <Plus className="w-3 h-3 group-hover:scale-110 transition-transform" />
                        <span>New Task</span>
                      </button>
                    ) : (
                      <div className="px-3 py-1.5 space-y-2">
                        <input
                          type="text"
                          value={newTaskName}
                          onChange={(e) => setNewTaskName(e.target.value)}
                          onKeyDown={(e) => {
                            if (e.key === 'Enter') handleCreateTask(project.id);
                            if (e.key === 'Escape') {
                              setCreatingTaskFor(null);
                              setNewTaskName('');
                            }
                          }}
                          placeholder="Task name..."
                          className="w-full px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-xs focus:outline-none focus:border-purple-400"
                          autoFocus
                        />
                        <div className="flex gap-2">
                          <button
                            onClick={() => handleCreateTask(project.id)}
                            className="flex-1 px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                          >
                            Create
                          </button>
                          <button
                            onClick={() => {
                              setCreatingTaskFor(null);
                              setNewTaskName('');
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
                    if (e.key === 'Enter') handleCreateProject();
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
                    onClick={handleCreateProject}
                    className="flex-1 px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
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
    </div>

    {/* Delete Confirmation Modal */}
    {deletingItem && (
      <div className="fixed inset-0 bg-black/70 backdrop-blur-sm z-[60] flex items-center justify-center p-4">
        <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-2xl shadow-2xl max-w-md w-full p-6">
          <div className="mb-4">
            <h3 className="text-lg font-bold text-white mb-2">Confirm Delete</h3>
            <p className="text-purple-200 text-sm">
              {deletingItem.type === 'project' 
                ? `Are you sure you want to delete this project? All ${projects.find(p => p.id === deletingItem.projectId)?.tasks.length || 0} tasks inside will also be deleted.`
                : 'Are you sure you want to delete this task?'}
            </p>
            <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
          </div>
          
          <div className="flex gap-3">
            <button
              onClick={() => {
                if (deletingItem.type === 'project') {
                  handleDeleteProject(deletingItem.projectId);
                } else {
                  handleDeleteTask(deletingItem.projectId, deletingItem.taskId!);
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
    taskId: null
  });
  const [hasLoadedInitialSelection, setHasLoadedInitialSelection] = useState(false);

  // Load last selection from localStorage on mount
  React.useEffect(() => {
    const savedSelection = localStorage.getItem(SELECTION_STORAGE_KEY);
    if (savedSelection) {
      try {
        const parsed = JSON.parse(savedSelection);
        setSelectionState(parsed);
      } catch (e) {
        console.error('Error loading last selection:', e);
      }
    }
    setHasLoadedInitialSelection(true);
  }, []);

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
    hasLoadedInitialSelection
  };
};

// Separate hook for project management functions to be used in the main app
export const useProjectManagement = () => {
  const { projects, createProject, createTask } = useProjectData();

  const createProjectAndTask = React.useCallback(async (
    projectName: string, 
    taskName: string
  ): Promise<{ projectId: string; taskId: string }> => {
    const project = await createProject(projectName);
    const task = await createTask(project.id, taskName);
    
    return {
      projectId: project.id,
      taskId: task.id
    };
  }, [createProject, createTask]);

  return {
    projects,  // Expose projects list
    createProjectAndTask,
    createTask  // Export this so it can be used separately
  };
};
