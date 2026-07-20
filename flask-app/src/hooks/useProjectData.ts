import { Mutex } from 'async-mutex';
import { useState, useEffect, useRef } from 'react';
import { Project, Experiment } from '../types';

import { HTTP_SERVER } from '../config';
const STORAGE_KEY = 'flask_copilot_projects';

// Types for interacting with the database
interface ExperimentResponse {
  id: string;
  projectId: string;
  name: string;
  createdAt: string;
  lastModified: string;
  data: any;
}

interface ProjectResponse {
  id: string;
  name: string;
  userId: string;
  createdAt: string;
  lastModified: string;
  experiments: ExperimentResponse[];
}

export interface ExperimentUpdate {
  id?: string;
  name?: string;
  isRunning?: boolean; // Track if experiment is currently computing
  [key: string]: any;
}

// This interface defines the contract for data sources (local storage, database, etc.)
interface ProjectDataSource {
  loadProjects: () => Promise<Project[]>;
  createProject: (name: string) => Promise<Project | null>;
  updateProject: (project: Project) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  createExperiment: (projectId: string, name: string) => Promise<Experiment | null>;
  updateExperiment: (projectId: string, experiment: ExperimentUpdate) => Promise<void>;
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
  createProject: (name: string) => Promise<Project | null>;
  updateProject: (project: Project) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  createExperiment: (projectId: string, name: string) => Promise<Experiment | null>;
  updateExperiment: (projectId: string, experiment: Experiment) => Promise<void>;
  deleteExperiment: (projectId: string, experimentId: string) => Promise<void>;
  setExperimentRunning: (
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ) => Promise<void>;
  refreshProjects: () => Promise<void>;
}

// Just need to flatten fields from 'data'.
function flattenExperiments(serverExperiments: ExperimentResponse[]): Experiment[] {
  return serverExperiments.map(({ data, projectId, ...sqlFields }) => {
    // problemId does NOT come along for this ride.
    const experiment = {
      ...data,
      ...sqlFields,
    };
    return experiment;
  });
}

// Reformat from the ProjectResponseWithExperiments to the frontend Project.
function flattenProjects(serverProjects: ProjectResponse[]): Project[] {
  return serverProjects.map(({ userId, experiments, ...rest }: ProjectResponse) => {
    // userId does NOT come along for this ride.
    return {
      ...rest,
      experiments: flattenExperiments(experiments),
    };
  });
}

const httpServerUrl = HTTP_SERVER;
class ServerDataSource implements ProjectDataSource {
  migrateDone: boolean = false;
  private readonly mutex: Mutex = new Mutex();

  async migrateProjectsFromLocalStorage(): Promise<void> {
    const stored = localStorage.getItem(STORAGE_KEY);
    // Nothing to migrate
    if (!stored || this.migrateDone) {
      this.migrateDone = true;
      return;
    }

    try {
      const ls_projects = JSON.parse(stored);
      const send_projects = ls_projects.map(({ name, experiments }: Project) => {
        return {
          name: name,
          experiments: experiments.map((exp: Experiment) => {
            return {
              data: exp,
            };
          }),
        };
      });
      const response = await fetch(`${httpServerUrl}/projects/migrate`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(send_projects),
      });
      if (!response.ok) {
        throw new Error(`migrateProjects response status: ${response.status}`);
      }
      // Everything is ok, blow up the localStorage.
      // FIXME (trb): UNCOMMENT THIS LINE WHEN YOU'RE READY TO GO LIVE
      // localStorage.removeItem(STORAGE_KEY);
      this.migrateDone = true;
    } catch (e) {
      console.error('Error migrating projects from localStorage:', e);
      return;
    }
  }

  async loadProjects(): Promise<Project[]> {
    try {
      // FIXME (trb): Perhaps there's a better spot to hook this is
      // in?? It will short-circuit, but ugh.

      await this.mutex.runExclusive(async () => {
        if (!this.migrateDone) await this.migrateProjectsFromLocalStorage();
      });

      const response = await fetch(`${httpServerUrl}/projects`);
      if (!response.ok) {
        throw new Error(`loadProjects response status: ${response.status}`);
      }
      const result = flattenProjects(await response.json());
      console.log('TRB LOADPROJECTS:', result);
      return result;
    } catch (e) {
      console.error('Error loading projects from server:', e);
      return [];
    }
  }

  async createProject(name: string): Promise<Project | null> {
    try {
      const response = await fetch(`${httpServerUrl}/projects`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          name: name,
        }),
      });
      if (!response.ok) {
        throw new Error(`createProject response status: ${await response.text()}`);
      }
      const project = await response.json();
      return project;
    } catch (e) {
      console.error('Error creating project:', e);
      return null;
    }
  }

  // FIXME (trb): This could be more efficiently implemented in the backend.
  async updateProject(project: Project): Promise<void> {
    try {
      // First we update the project metadata
      const response = await fetch(`${httpServerUrl}/projects/${project.id}`, {
        method: 'PUT',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(project),
      });
      if (!response.ok) {
        throw new Error(`updateProject response status: ${response.status}`);
      }

      // Then we update each experiment
      for (const exp of project.experiments) {
        await this.updateExperiment(project.id, exp);
      }
    } catch (e) {
      console.error('Error updating project:', e);
    }
  }

  async deleteProject(projectId: string): Promise<void> {
    try {
      const response = await fetch(`${httpServerUrl}/projects/${projectId}`, {
        method: 'DELETE',
      });
      if (!response.ok) {
        throw new Error(`deleteProject response status: ${response.status}`);
      }
    } catch (e) {
      console.error('Error deleting project:', e);
    }
  }

  async createExperiment(projectId: string, name: string): Promise<Experiment | null> {
    try {
      const response = await fetch(`${httpServerUrl}/projects/${projectId}/experiments`, {
        method: 'POST',
        headers: {
          Accept: 'application/json',
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ name: name }),
      });
      if (!response.ok) {
        throw new Error(`createExperiment response status: ${response.status}`);
      }
      const result = await response.json();
      return result;
    } catch (e) {
      console.error('Error creating experiment:', e);
    }
    return null;
  }

  async updateExperiment(projectId: string, experiment: ExperimentUpdate): Promise<void> {
    try {
      const response = await fetch(
        `${httpServerUrl}/projects/${projectId}/experiments/${experiment.id}`,
        {
          method: 'PUT',
          headers: {
            Accept: 'application/json',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ data: experiment }),
        }
      );
      if (!response.ok) {
        throw new Error(`updateExperiment response status: ${response.status}`);
      }
    } catch (e) {
      console.error('Error updating experiment:', e);
    }
  }

  async deleteExperiment(projectId: string, experimentId: string): Promise<void> {
    try {
      const response = await fetch(
        `${httpServerUrl}/projects/${projectId}/experiments/${experimentId}`,
        {
          method: 'DELETE',
        }
      );
      if (!response.ok) {
        throw new Error(`deleteExperiment response status: ${response.status}`);
      }
    } catch (e) {
      console.error('Error deleting experiment:', e);
    }
  }

  async setExperimentRunning(
    projectId: string,
    experimentId: string,
    isRunning: boolean
  ): Promise<void> {
    await this.updateExperiment(projectId, { id: experimentId, isRunning: isRunning });
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

// Hook for managing project data
export const useProjectData = () => {
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(true);

  const projectsRef = useRef(projects);
  useEffect(() => {
    projectsRef.current = projects;
  }, [projects]);

  // TODO(later): Swap this to use ServerDataSource
  const dataSource: ProjectDataSource = new ServerDataSource();

  useEffect(() => {
    loadProjects();
  }, []);

  const loadProjects = async () => {
    setLoading(true);
    try {
      const data = await dataSource.loadProjects();
      setProjects(data);
      projectsRef.current = data;
    } catch (error) {
      console.error('Error loading projects:', error);
    } finally {
      setLoading(false);
    }
  };

  const createProject = async (name: string): Promise<Project | null> => {
    const newProject = await dataSource.createProject(name);
    await loadProjects();
    return newProject;
  };

  const updateProject = async (project: Project): Promise<void> => {
    await dataSource.updateProject(project);
    await loadProjects();
  };

  const deleteProject = async (projectId: string): Promise<void> => {
    await dataSource.deleteProject(projectId);
    await loadProjects();
  };

  const createExperiment = async (projectId: string, name: string): Promise<Experiment | null> => {
    const newExperiment = await dataSource.createExperiment(projectId, name);
    await loadProjects();
    return newExperiment;
  };

  const updateExperiment = async (projectId: string, experiment: Experiment) => {
    await dataSource.updateExperiment(projectId, experiment);
    await loadProjects();
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
    await loadProjects();
  };

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
  };
};
