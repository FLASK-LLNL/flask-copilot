import { useState, useEffect, useRef } from 'react';
import { Project, Experiment } from '../types';
import { HTTP_SERVER } from '../config';

const STORAGE_KEY = 'flask_copilot_projects';

// This interface defines the contract for data sources (local storage, database, etc.)
interface ProjectDataSource {
  loadProjects: () => Promise<Project[]>;
  saveProjects: (projects: Project[]) => Promise<void>;
  createProject: (name: string) => Promise<Project>;
  updateProject: (project: Project) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  createExperiment: (projectId: string, name: string) => Promise<Experiment>;
  updateExperiment: (projectId: string, experiment: Experiment) => Promise<void>;
  deleteExperiment: (projectId: string, experimentId: string) => Promise<void>;
  setExperimentRunning: (
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ) => Promise<void>;
}

export interface ProjectData {
  projectsRef: React.RefObject<Project[]>;
  projects: Project[];
  loading: boolean;
  createProject: (name: string) => Promise<Project>;
  updateProject: (project: Project) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  createExperiment: (projectId: string, name: string) => Promise<Experiment>;
  updateExperiment: (projectId: string, experiment: Experiment) => Promise<void>;
  deleteExperiment: (projectId: string, experimentId: string) => Promise<void>;
  setExperimentRunning: (
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ) => Promise<void>;
  refreshProjects: () => Promise<void>;
  clearAllProjects: () => Promise<void>;
}

// LocalStorage implementation
class LocalStorageDataSource implements ProjectDataSource {
  async loadProjects(): Promise<Project[]> {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (!stored) return [];
    try {
      return JSON.parse(stored);
    } catch (e) {
      console.error('Error loading projects from localStorage:', e);
      return [];
    }
  }

  async saveProjects(projects: Project[]): Promise<void> {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(projects));
  }

  async createProject(name: string): Promise<Project> {
    const newProject: Project = {
      id: `project_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString(),
      experiments: [],
    };

    const projects = await this.loadProjects();
    projects.push(newProject);
    await this.saveProjects(projects);

    return newProject;
  }

  async updateProject(project: Project): Promise<void> {
    const projects = await this.loadProjects();
    const index = projects.findIndex((p) => p.id === project.id);
    if (index !== -1) {
      projects[index] = { ...project, lastModified: new Date().toISOString() };
      await this.saveProjects(projects);
    }
  }

  async deleteProject(projectId: string): Promise<void> {
    const projects = await this.loadProjects();
    const filtered = projects.filter((p) => p.id !== projectId);
    await this.saveProjects(filtered);
  }

  async createExperiment(projectId: string, name: string): Promise<Experiment> {
    const newExperiment: Experiment = {
      id: `exp_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString(),
    };

    const projects = await this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      project.experiments.push(newExperiment);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }

    return newExperiment;
  }

  async updateExperiment(projectId: string, experiment: Experiment): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      const expIndex = project.experiments.findIndex((e) => e.id === experiment.id);
      if (expIndex !== -1) {
        project.experiments[expIndex] = { ...experiment, lastModified: new Date().toISOString() };
        project.lastModified = new Date().toISOString();
        await this.saveProjects(projects);
      }
    }
  }

  async deleteExperiment(projectId: string, experimentId: string): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      project.experiments = project.experiments.filter((e) => e.id !== experimentId);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }
  }

  async setExperimentRunning(
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      const experiment = project.experiments.find((e) => e.id === experimentId);
      if (experiment) {
        experiment.isRunning = isRunning;
        experiment.lastModified = new Date().toISOString();
        project.lastModified = new Date().toISOString();
        await this.saveProjects(projects);
      }
    }
  }
}

// Database implementation with localStorage fallback
class DatabaseDataSource implements ProjectDataSource {
  private useLocalStorageFallback = false;
  private localStorageSource = new LocalStorageDataSource();

  /** Load all projects from the database. */
  async loadProjects(): Promise<Project[]> {
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects`);
        if (response.ok) {
          const projects = await response.json();
          console.log(`[DatabaseDataSource] Loaded ${projects.length} projects from database`);
          localStorage.setItem(STORAGE_KEY, JSON.stringify(projects));
          return projects;
        } else {
          console.warn(`[DatabaseDataSource] Database request failed: ${response.status}`);
          this.useLocalStorageFallback = true;
        }
      } catch (error) {
        console.warn('[DatabaseDataSource] Database unavailable, falling back to localStorage:', error);
        this.useLocalStorageFallback = true;
      }
    }

    // Fallback to localStorage
    console.log('[DatabaseDataSource] Using localStorage fallback');
    return this.localStorageSource.loadProjects();
  }

  /** Delete all projects and experiments from the database and clear localStorage. */
  async clearAllProjects(): Promise<void> {
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/all`, {
          method: 'DELETE',
        });
        if (!response.ok) {
          console.warn('[DatabaseDataSource] Failed to delete all projects from database');
        }
      } catch (error) {
        console.warn('[DatabaseDataSource] Database unavailable during clear:', error);
      }
    }
    // Clear localStorage caches
    localStorage.removeItem(STORAGE_KEY);
    localStorage.removeItem('flask_copilot_last_selection');
  }

  async saveProjects(projects: Project[]): Promise<void> {
    // Always save to localStorage as cache
    await this.localStorageSource.saveProjects(projects);
  }

  async createProject(name: string): Promise<Project> {
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ name }),
        });
        if (response.ok) {
          const project = await response.json();
          // Update localStorage cache
          const projects = await this.localStorageSource.loadProjects();
          projects.push(project);
          await this.localStorageSource.saveProjects(projects);
          return project;
        }
      } catch (error) {
        console.warn('Database unavailable, using localStorage:', error);
        this.useLocalStorageFallback = true;
      }
    }
    
    return this.localStorageSource.createProject(name);
  }

  async updateProject(project: Project): Promise<void> {
    // Update localStorage cache first
    await this.localStorageSource.updateProject(project);
    
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${project.id}`, {
          method: 'PUT',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ name: project.name }),
        });
        if (!response.ok) {
          console.warn('Failed to update project in database');
        }
      } catch (error) {
        console.warn('Database unavailable:', error);
      }
    }
  }

  async deleteProject(projectId: string): Promise<void> {
    // Delete from localStorage cache
    await this.localStorageSource.deleteProject(projectId);
    
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${projectId}`, {
          method: 'DELETE',
        });
        if (!response.ok) {
          console.warn('Failed to delete project from database');
        }
      } catch (error) {
        console.warn('Database unavailable:', error);
      }
    }
  }

  async createExperiment(projectId: string, name: string): Promise<Experiment> {
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${projectId}/experiments`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ name }),
        });
        if (response.ok) {
          const experiment = await response.json();
          // Update localStorage cache
          const projects = await this.localStorageSource.loadProjects();
          const project = projects.find(p => p.id === projectId);
          if (project) {
            project.experiments.push(experiment);
            project.lastModified = new Date().toISOString();
            await this.localStorageSource.saveProjects(projects);
          }
          return experiment;
        } else {
          // HTTP error from the server (e.g. 500 from concurrent DB write).
          // Do NOT fall back to localStorage: the experiment would appear locally
          // but be invisible to the DB, causing data loss on next page load.
          const text = await response.text().catch(() => response.status.toString());
          throw new Error(`Failed to create experiment in database (HTTP ${response.status}): ${text}`);
        }
      } catch (error) {
        if (error instanceof Error && error.message.startsWith('Failed to create experiment')) {
          // Re-throw DB HTTP errors so the caller can surface them to the user.
          throw error;
        }
        // Network error (DB unreachable) – fall back to localStorage.
        console.warn('Database unavailable, using localStorage for experiment creation:', error);
        this.useLocalStorageFallback = true;
      }
    }
    
    return this.localStorageSource.createExperiment(projectId, name);
  }

  async updateExperiment(projectId: string, experiment: Experiment): Promise<void> {
    // Update localStorage cache first
    await this.localStorageSource.updateExperiment(projectId, experiment);
    
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${projectId}/experiments/${experiment.id}`, {
          method: 'PUT',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ experiment }),
        });
        if (!response.ok) {
          console.warn('Failed to update experiment in database');
        }
      } catch (error) {
        console.warn('Database unavailable:', error);
      }
    }
  }

  async deleteExperiment(projectId: string, experimentId: string): Promise<void> {
    // Delete from localStorage cache
    await this.localStorageSource.deleteExperiment(projectId, experimentId);
    
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${projectId}/experiments/${experimentId}`, {
          method: 'DELETE',
        });
        if (!response.ok) {
          console.warn('Failed to delete experiment from database');
        }
      } catch (error) {
        console.warn('Database unavailable:', error);
      }
    }
  }

  async setExperimentRunning(projectId: string, experimentId: string, isRunning: boolean): Promise<void> {
    // Update localStorage cache
    await this.localStorageSource.setExperimentRunning(projectId, experimentId, isRunning);
    
    if (!this.useLocalStorageFallback) {
      try {
        const response = await fetch(`${HTTP_SERVER}/api/projects/${projectId}/experiments/${experimentId}/running?is_running=${isRunning}`, {
          method: 'PUT',
        });
        if (!response.ok) {
          console.warn('Failed to update experiment running status in database');
        }
      } catch (error) {
        console.warn('Database unavailable:', error);
      }
    }
  }
}

// TODO: Future server implementation
// class ServerDataSource implements ProjectDataSource {
//   async loadProjects(): Promise<Project[]> {
//     const response = await fetch('/api/projects');
//     return response.json();
//   }
//   // ...
// }

// Synchronous read of projects from localStorage so the data is available
// before any effects run (avoids race with WebSocket onopen).
const loadProjectsSync = (): Project[] => {
  try {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (stored) {
      return JSON.parse(stored);
    }
  } catch (e) {
    console.error('Error reading projects from localStorage:', e);
  }
  return [];
};

// Hook for managing project data
export const useProjectData = () => {
  const [projects, setProjects] = useState<Project[]>(loadProjectsSync);
  const [loading, setLoading] = useState(false);

  const projectsRef = useRef(projects);
  useEffect(() => {
    projectsRef.current = projects;
  }, [projects]);

  // Use database-backed data source with localStorage fallback
  const dataSource = useRef(new DatabaseDataSource()).current;

  // Track whether the initial load from the database has completed.
  // Subsequent 2-second polls should NOT flip the loading flag to avoid
  // the sidebar flickering a spinner every cycle.
  const initialLoadDoneRef = useRef(false);

  const loadProjects = async () => {
    if (!initialLoadDoneRef.current) {
      setLoading(true);
    }
    try {
      const data = await dataSource.loadProjects();

      // Only update state if data actually changed (avoids re-render cascades)
      const newJson = JSON.stringify(data);
      const oldJson = JSON.stringify(projectsRef.current);
      if (newJson !== oldJson) {
        setProjects(data);
        projectsRef.current = data;
      }
    } catch (error) {
      console.error('Error loading projects:', error);
    } finally {
      if (!initialLoadDoneRef.current) {
        setLoading(false);
        initialLoadDoneRef.current = true;
      }
    }
  };

  const createProject = async (name: string): Promise<Project> => {
    const newProject = await dataSource.createProject(name);
    await loadProjects();
    return newProject;
  };

  const updateProject = async (project: Project) => {
    await dataSource.updateProject(project);
  };

  const deleteProject = async (projectId: string) => {
    await dataSource.deleteProject(projectId);
    await loadProjects();
  };

  const createExperiment = async (projectId: string, name: string): Promise<Experiment> => {
    const newExperiment = await dataSource.createExperiment(projectId, name);
    await loadProjects();
    return newExperiment;
  };

  const updateExperiment = async (projectId: string, experiment: Experiment) => {
    // Immediately update projectsRef so any synchronous reads in the same
    // event loop tick (e.g., loadContextFromExperiment after a sidebar
    // click) see the freshly-saved data instead of waiting for the next
    // 2-second poll.
    const updatedProjects = projectsRef.current.map(p => {
      if (p.id !== projectId) return p;
      return {
        ...p,
        experiments: p.experiments.map(e =>
          e.id === experiment.id
            ? { ...experiment, lastModified: new Date().toISOString() }
            : e
        ),
        lastModified: new Date().toISOString(),
      };
    });
    projectsRef.current = updatedProjects;
    setProjects(updatedProjects);

    await dataSource.updateExperiment(projectId, experiment);
  };

  const deleteExperiment = async (projectId: string, experimentId: string) => {
    await dataSource.deleteExperiment(projectId, experimentId);
    await loadProjects();
  };

  const setExperimentRunning = async (
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ) => {
    await dataSource.setExperimentRunning(projectId, experimentId, isRunning);
  };

  const clearAllProjects = async () => {
    await dataSource.clearAllProjects();
    setProjects([]);
    projectsRef.current = [];
  };

  // Auto-refresh projects from database every 2 seconds to pick up changes
  // from other browsers.
  useEffect(() => {
    console.log('[useProjectData] Initial load of projects');
    loadProjects();

    const refreshInterval = setInterval(() => {
      loadProjects();
    }, 2000);

    return () => clearInterval(refreshInterval);
  }, []); // Empty deps - only run once on mount

  return {
    projectsRef,
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    setExperimentRunning,
    refreshProjects: loadProjects,
    clearAllProjects
  };
};
