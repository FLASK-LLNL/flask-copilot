import React, { useState } from 'react';
import { ChevronRight, ChevronDown, ChevronLeft, Plus, FolderOpen, FileText, Loader2, Edit2, Trash2 } from 'lucide-react';
import { Experiment, Task, ExperimentSelection } from '../types';
import { useExperimentData } from '../hooks/useExperimentData';

interface ExperimentSidebarProps {
  isOpen: boolean;
  onToggle: () => void;
  selection: ExperimentSelection;
  onSelectionChange: (selection: ExperimentSelection) => void;
  onLoadContext: (experimentId: string, taskId: string | null) => void;
  isComputing?: boolean;  // Track if main app is currently computing
  hasLoadedInitialSelection?: boolean;  // Track if initial selection from localStorage has been loaded
  onCreateExperimentAndTask?: (experimentName: string, taskName: string) => Promise<{ experimentId: string; taskId: string }>;
}

export const ExperimentSidebar: React.FC<ExperimentSidebarProps> = ({
  isOpen,
  onToggle,
  selection,
  onSelectionChange,
  onLoadContext,
  isComputing = false,
  hasLoadedInitialSelection = false,
  onCreateExperimentAndTask
}) => {
  const {
    experiments,
    loading,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    createTask,
    updateTask,
    deleteTask,
    setTaskRunning
  } = useExperimentData();

  const SIDEBAR_WIDTH_STORAGE_KEY = 'flask_copilot_sidebar_width';
  const MIN_WIDTH = 200;
  const MAX_WIDTH = 1600;
  const DEFAULT_WIDTH = 320;
  const COLLAPSE_THRESHOLD = 5;

  const [expandedExperiments, setExpandedExperiments] = useState<Set<string>>(new Set());
  const [creatingExperiment, setCreatingExperiment] = useState(false);
  const [newExperimentName, setNewExperimentName] = useState('');
  const [creatingTaskFor, setCreatingTaskFor] = useState<string | null>(null);
  const [newTaskName, setNewTaskName] = useState('');
  const [editingExperiment, setEditingExperiment] = useState<string | null>(null);
  const [editExperimentName, setEditExperimentName] = useState('');
  const [editingTask, setEditingTask] = useState<{ experimentId: string; taskId: string } | null>(null);
  const [editTaskName, setEditTaskName] = useState('');
  const [deletingItem, setDeletingItem] = useState<{ type: 'experiment' | 'task'; experimentId: string; taskId?: string } | null>(null);
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

  // Auto-expand and load context when initial selection is loaded from localStorage
  React.useEffect(() => {
    if (hasLoadedInitialSelection && !loading && !hasRestoredSelection && experiments.length > 0 && selection.experimentId) {
      // Verify the selection still exists in the experiments
      const experiment = experiments.find(p => p.id === selection.experimentId);
      if (experiment) {
        // Auto-expand the experiment
        setExpandedExperiments(prev => new Set(prev).add(selection.experimentId!));
        
        // Load the context
        onLoadContext(selection.experimentId, selection.taskId);
        setHasRestoredSelection(true);
      }
    }
  }, [hasLoadedInitialSelection, loading, experiments, selection, onLoadContext, hasRestoredSelection]);

  // Update task running status when isComputing changes
  React.useEffect(() => {
    if (selection.experimentId && selection.taskId) {
      setTaskRunning(selection.experimentId, selection.taskId, isComputing);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isComputing, selection.experimentId, selection.taskId]);

  const toggleExperimentExpanded = (experimentId: string) => {
    const newExpanded = new Set(expandedExperiments);
    if (newExpanded.has(experimentId)) {
      newExpanded.delete(experimentId);
    } else {
      newExpanded.add(experimentId);
    }
    setExpandedExperiments(newExpanded);
  };

  const handleExperimentClick = (experiment: Experiment) => {
    // Auto-select the last task if the experiment has any
    const lastTask = experiment.tasks.length > 0 
      ? experiment.tasks[experiment.tasks.length - 1] 
      : null;
    
    const newSelection: ExperimentSelection = {
      experimentId: experiment.id,
      taskId: lastTask?.id || null
    };
    
    onSelectionChange(newSelection);
    onLoadContext(experiment.id, lastTask?.id || null);
    
    // Auto-expand when clicking a experiment
    if (!expandedExperiments.has(experiment.id)) {
      toggleExperimentExpanded(experiment.id);
    }
  };

  const handleTaskClick = (experiment: Experiment, task: Task) => {
    const newSelection: ExperimentSelection = {
      experimentId: experiment.id,
      taskId: task.id
    };
    onSelectionChange(newSelection);
    onLoadContext(experiment.id, task.id);
  };

  const handleCreateExperiment = async () => {
    if (!newExperimentName.trim()) return;
    
    try {
      const experiment = await createExperiment(newExperimentName);
      setNewExperimentName('');
      setCreatingExperiment(false);
      // Auto-select and expand the new experiment
      onSelectionChange({ experimentId: experiment.id, taskId: null });
      setExpandedExperiments(prev => new Set(prev).add(experiment.id));
    } catch (error) {
      console.error('Error creating experiment:', error);
      alert('Failed to create experiment');
    }
  };

  const handleCreateTask = async (experimentId: string) => {
    if (!newTaskName.trim()) return;
    
    try {
      const task = await createTask(experimentId, newTaskName);
      setNewTaskName('');
      setCreatingTaskFor(null);
      // Auto-select the new task
      onSelectionChange({ experimentId, taskId: task.id });
      onLoadContext(experimentId, task.id);
    } catch (error) {
      console.error('Error creating task:', error);
      alert('Failed to create task');
    }
  };

  const handleEditExperiment = (experiment: Experiment) => {
    setEditingExperiment(experiment.id);
    setEditExperimentName(experiment.name);
  };

  const handleSaveExperimentEdit = async (experiment: Experiment) => {
    if (!editExperimentName.trim()) return;
    
    try {
      await updateExperiment({ ...experiment, name: editExperimentName });
      setEditingExperiment(null);
      setEditExperimentName('');
    } catch (error) {
      console.error('Error updating experiment:', error);
      alert('Failed to update experiment');
    }
  };

  const handleDeleteExperiment = async (experimentId: string) => {
    try {
      await deleteExperiment(experimentId);
      setDeletingItem(null);
      // Clear selection if deleted experiment was selected
      if (selection.experimentId === experimentId) {
        onSelectionChange({ experimentId: null, taskId: null });
      }
    } catch (error) {
      console.error('Error deleting experiment:', error);
      alert('Failed to delete experiment');
    }
  };

  const handleEditTask = (experimentId: string, task: Task) => {
    setEditingTask({ experimentId, taskId: task.id });
    setEditTaskName(task.name);
  };

  const handleSaveTaskEdit = async (experimentId: string, task: Task) => {
    if (!editTaskName.trim()) return;
    
    try {
      await updateTask(experimentId, { ...task, name: editTaskName });
      setEditingTask(null);
      setEditTaskName('');
    } catch (error) {
      console.error('Error updating task:', error);
      alert('Failed to update task');
    }
  };

  const handleDeleteTask = async (experimentId: string, taskId: string) => {
    try {
      await deleteTask(experimentId, taskId);
      setDeletingItem(null);
      // Clear selection if deleted task was selected
      if (selection.taskId === taskId) {
        onSelectionChange({ experimentId: selection.experimentId, taskId: null });
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
            title="Open Experiments"
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
                  {deletingItem.type === 'experiment' 
                    ? `Are you sure you want to delete this experiment? All ${experiments.find(p => p.id === deletingItem.experimentId)?.tasks.length || 0} tasks inside will also be deleted.`
                    : 'Are you sure you want to delete this task?'}
                </p>
                <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
              </div>
              
              <div className="flex gap-3">
                <button
                  onClick={() => {
                    if (deletingItem.type === 'experiment') {
                      handleDeleteExperiment(deletingItem.experimentId);
                    } else {
                      handleDeleteTask(deletingItem.experimentId, deletingItem.taskId!);
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
            <h2 className="text-lg font-semibold text-white">Experiments</h2>
          </div>
          <button
            onClick={onToggle}
            className="p-1.5 text-purple-300 hover:text-white hover:bg-purple-600/30 rounded-lg transition-all"
            title="Collapse"
          >
            <ChevronLeft className="w-5 h-5" />
          </button>
        </div>
        
        {/* No Experiment/Task Selected Warning */}
        {(!selection.experimentId || !selection.taskId) && (
          <div className="mt-2 px-3 py-2 bg-amber-500/20 border border-amber-400/50 rounded-lg">
            <p className="text-amber-200 text-xs">
              {!selection.experimentId 
                ? '⚠️ No experiment selected. A new experiment will be created when you run.'
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
            {/* Experiments List */}
            {experiments.map((experiment) => (
              <div key={experiment.id} className="space-y-0.5">
                {/* Experiment Item */}
                <div className="flex items-center group relative">
                  <button
                    onClick={() => toggleExperimentExpanded(experiment.id)}
                    className="p-1 text-purple-300 hover:text-white transition-colors flex-shrink-0"
                  >
                    {expandedExperiments.has(experiment.id) ? (
                      <ChevronDown className="w-4 h-4" />
                    ) : (
                      <ChevronRight className="w-4 h-4" />
                    )}
                  </button>
                  
                  {editingExperiment === experiment.id ? (
                    <div className="flex-1 flex items-center gap-1 px-2">
                      <input
                        type="text"
                        value={editExperimentName}
                        onChange={(e) => setEditExperimentName(e.target.value)}
                        onKeyDown={(e) => {
                          if (e.key === 'Enter') handleSaveExperimentEdit(experiment);
                          if (e.key === 'Escape') {
                            setEditingExperiment(null);
                            setEditExperimentName('');
                          }
                        }}
                        className="flex-1 px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-sm focus:outline-none focus:border-purple-400"
                        autoFocus
                      />
                      <button
                        onClick={() => handleSaveExperimentEdit(experiment)}
                        className="px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                      >
                        Save
                      </button>
                      <button
                        onClick={() => {
                          setEditingExperiment(null);
                          setEditExperimentName('');
                        }}
                        className="px-2 py-1 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                      >
                        Cancel
                      </button>
                    </div>
                  ) : (
                    <>
                      <button
                        onClick={() => handleExperimentClick(experiment)}
                        className={`flex-1 px-3 py-2 text-left text-sm rounded-lg transition-all flex items-center gap-2 min-w-0 ${
                          selection.experimentId === experiment.id
                            ? 'bg-purple-600/50 text-white font-medium'
                            : 'text-purple-200 hover:text-white hover:bg-purple-600/30'
                        }`}
                      >
                        <FolderOpen className="w-4 h-4 flex-shrink-0" />
                        <span className="truncate flex-1">{experiment.name}</span>
                        {experiment.tasks.length > 0 && (
                          <span className="text-xs text-purple-400 flex-shrink-0">
                            {experiment.tasks.length}
                          </span>
                        )}
                      </button>
                      {/* Edit/Delete buttons with background overlay */}
                      <div className="absolute right-1 flex items-center gap-0.5 opacity-0 group-hover:opacity-100 transition-opacity">
                        <div className="bg-gradient-to-l from-slate-800 via-slate-800 to-transparent pl-8 pr-1 py-1 flex items-center gap-0.5">
                          <button
                            onClick={(e) => {
                              e.stopPropagation();
                              handleEditExperiment(experiment);
                            }}
                            className="p-1 text-purple-300 hover:text-white hover:bg-purple-600/50 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                            title="Edit experiment"
                          >
                            <Edit2 className="w-3.5 h-3.5" />
                          </button>
                          <button
                            onClick={(e) => {
                              e.stopPropagation();
                              setDeletingItem({ type: 'experiment', experimentId: experiment.id });
                            }}
                            className="p-1 text-purple-300 hover:text-red-400 hover:bg-red-500/30 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                            title="Delete experiment"
                          >
                            <Trash2 className="w-3.5 h-3.5" />
                          </button>
                        </div>
                      </div>
                    </>
                  )}
                </div>

                {/* Tasks List */}
                {expandedExperiments.has(experiment.id) && (
                  <div className="ml-6 space-y-0.5">
                    {/* Task Items */}
                    {experiment.tasks.map((task) => (
                      <div key={task.id} className="flex items-center group relative">
                        {editingTask?.taskId === task.id ? (
                          <div className="flex-1 flex items-center gap-1 px-2 py-1">
                            <input
                              type="text"
                              value={editTaskName}
                              onChange={(e) => setEditTaskName(e.target.value)}
                              onKeyDown={(e) => {
                                if (e.key === 'Enter') handleSaveTaskEdit(experiment.id, task);
                                if (e.key === 'Escape') {
                                  setEditingTask(null);
                                  setEditTaskName('');
                                }
                              }}
                              className="flex-1 px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-xs focus:outline-none focus:border-purple-400"
                              autoFocus
                            />
                            <button
                              onClick={() => handleSaveTaskEdit(experiment.id, task)}
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
                              onClick={() => handleTaskClick(experiment, task)}
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
                                    handleEditTask(experiment.id, task);
                                  }}
                                  className="p-0.5 text-purple-300 hover:text-white hover:bg-purple-600/50 rounded transition-all backdrop-blur-sm bg-slate-800/80"
                                  title="Edit task"
                                >
                                  <Edit2 className="w-3 h-3" />
                                </button>
                                <button
                                  onClick={(e) => {
                                    e.stopPropagation();
                                    setDeletingItem({ type: 'task', experimentId: experiment.id, taskId: task.id });
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
                    {creatingTaskFor !== experiment.id ? (
                      <button
                        onClick={() => setCreatingTaskFor(experiment.id)}
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
                            if (e.key === 'Enter') handleCreateTask(experiment.id);
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
                            onClick={() => handleCreateTask(experiment.id)}
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

            {/* Create New Experiment Button (at bottom) */}
            {!creatingExperiment ? (
              <button
                onClick={() => setCreatingExperiment(true)}
                className="w-full px-3 py-2 text-left text-sm text-purple-300 hover:text-white hover:bg-purple-600/30 rounded-lg transition-all flex items-center gap-2 group"
              >
                <Plus className="w-4 h-4 group-hover:scale-110 transition-transform" />
                <span>New Experiment</span>
              </button>
            ) : (
              <div className="px-3 py-2 space-y-2">
                <input
                  type="text"
                  value={newExperimentName}
                  onChange={(e) => setNewExperimentName(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === 'Enter') handleCreateExperiment();
                    if (e.key === 'Escape') {
                      setCreatingExperiment(false);
                      setNewExperimentName('');
                    }
                  }}
                  placeholder="Experiment name..."
                  className="w-full px-2 py-1 bg-white/10 border border-purple-400/50 rounded text-white text-sm focus:outline-none focus:border-purple-400"
                  autoFocus
                />
                <div className="flex gap-2">
                  <button
                    onClick={handleCreateExperiment}
                    className="flex-1 px-2 py-1 bg-purple-600 hover:bg-purple-500 text-white text-xs rounded transition-colors"
                  >
                    Create
                  </button>
                  <button
                    onClick={() => {
                      setCreatingExperiment(false);
                      setNewExperimentName('');
                    }}
                    className="flex-1 px-2 py-1 bg-white/10 hover:bg-white/20 text-white text-xs rounded transition-colors"
                  >
                    Cancel
                  </button>
                </div>
              </div>
            )}

            {experiments.length === 0 && !creatingExperiment && (
              <div className="text-center py-8 text-purple-300 text-sm">
                <p>No experiments yet</p>
                <p className="text-xs mt-1">Click "New Experiment" to get started</p>
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
              {deletingItem.type === 'experiment' 
                ? `Are you sure you want to delete this experiment? All ${experiments.find(p => p.id === deletingItem.experimentId)?.tasks.length || 0} tasks inside will also be deleted.`
                : 'Are you sure you want to delete this task?'}
            </p>
            <p className="text-purple-300 text-xs mt-2">This action cannot be undone.</p>
          </div>
          
          <div className="flex gap-3">
            <button
              onClick={() => {
                if (deletingItem.type === 'experiment') {
                  handleDeleteExperiment(deletingItem.experimentId);
                } else {
                  handleDeleteTask(deletingItem.experimentId, deletingItem.taskId!);
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
export const useExperimentSidebar = () => {
  const SELECTION_STORAGE_KEY = 'flask_copilot_last_selection';
  
  const [isOpen, setIsOpen] = useState(true);
  const [selection, setSelectionState] = useState<ExperimentSelection>({
    experimentId: null,
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
  const setSelection = React.useCallback((newSelection: ExperimentSelection) => {
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

// Separate hook for experiment management functions to be used in the main app
export const useExperimentManagement = () => {
  const { experiments, createExperiment, createTask } = useExperimentData();

  const createExperimentAndTask = React.useCallback(async (
    experimentName: string, 
    taskName: string
  ): Promise<{ experimentId: string; taskId: string }> => {
    const experiment = await createExperiment(experimentName);
    const task = await createTask(experiment.id, taskName);
    
    return {
      experimentId: experiment.id,
      taskId: task.id
    };
  }, [createExperiment, createTask]);

  return {
    experiments,  // Expose experiments list
    createExperimentAndTask,
    createTask  // Export this so it can be used separately
  };
};
