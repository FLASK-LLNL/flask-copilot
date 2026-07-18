import { useState, useEffect, useRef } from 'react';
import { Project, Experiment } from '../types';

import { HTTP_SERVER } from '../config';
const STORAGE_KEY = 'flask_copilot_projects';

// This interface defines the contract for data sources (local storage, database, etc.)
interface ProjectDataSource {
  loadProjects: () => Project[];
  saveProjects: (projects: Project[]) => void;
  createProject: (name: string) => Project;
  updateProject: (project: Project) => void;
  deleteProject: (projectId: string) => void;
  createExperiment: (projectId: string, name: string) => Experiment;
  updateExperiment: (projectId: string, experiment: Experiment) => void;
  deleteExperiment: (projectId: string, experimentId: string) => void;
  setExperimentRunning: (projectId: string, experimentId: string, isRunning: boolean) => void;
}

export interface ProjectData {
  projectsRef: React.RefObject<Project[]>;
  projects: Project[];
  loading: boolean;
  createProject: (name: string) => Project;
  updateProject: (project: Project) => void;
  deleteProject: (projectId: string) => void;
  createExperiment: (projectId: string, name: string) => Experiment;
  updateExperiment: (projectId: string, experiment: Experiment) => void;
  deleteExperiment: (projectId: string, experimentId: string) => void;
  setExperimentRunning: (projectId: string, experimentId: string, isRunning: boolean) => void;
  refreshProjects: () => void;
}

// LocalStorage implementation
class LocalStorageDataSource implements ProjectDataSource {
  loadProjects(): Project[] {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (!stored) return [];
    try {
      return JSON.parse(stored);
    } catch (e) {
      console.error('Error loading projects from localStorage:', e);
      return [];
    }
  }

  saveProjects(projects: Project[]): void {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(projects));
    } catch (e) {
      console.log('saveProjects error:', e);
    }
  }

  createProject(name: string): Project {
    const newProject: Project = {
      id: `project_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString(),
      experiments: [],
    };

    const projects = this.loadProjects();
    projects.push(newProject);
    this.saveProjects(projects);

    return newProject;
  }

  updateProject(project: Project): void {
    const projects = this.loadProjects();
    const index = projects.findIndex((p) => p.id === project.id);
    if (index !== -1) {
      projects[index] = { ...project, lastModified: new Date().toISOString() };
      this.saveProjects(projects);
    }
  }

  deleteProject(projectId: string): void {
    const projects = this.loadProjects();
    const filtered = projects.filter((p) => p.id !== projectId);
    this.saveProjects(filtered);
  }

  createExperiment(projectId: string, name: string): Experiment {
    const newExperiment: Experiment = {
      id: `exp_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString(),
    };

    const projects = this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      project.experiments.push(newExperiment);
      project.lastModified = new Date().toISOString();
      this.saveProjects(projects);
    }

    return newExperiment;
  }

  updateExperiment(projectId: string, experiment: Experiment): void {
    const projects = this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      const expIndex = project.experiments.findIndex((e) => e.id === experiment.id);
      if (expIndex !== -1) {
        project.experiments[expIndex] = { ...experiment, lastModified: new Date().toISOString() };
        project.lastModified = new Date().toISOString();
        this.saveProjects(projects);
      }
    }
  }

  deleteExperiment(projectId: string, experimentId: string): void {
    const projects = this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      project.experiments = project.experiments.filter((e) => e.id !== experimentId);
      project.lastModified = new Date().toISOString();
      this.saveProjects(projects);
    }
  }

  setExperimentRunning(projectId: string, experimentId: string, isRunning: boolean): void {
    const projects = this.loadProjects();
    const project = projects.find((p) => p.id === projectId);
    if (project) {
      const experiment = project.experiments.find((e) => e.id === experimentId);
      if (experiment) {
        experiment.isRunning = isRunning;
        experiment.lastModified = new Date().toISOString();
        project.lastModified = new Date().toISOString();
        this.saveProjects(projects);
      }
    }
  }
}

// Need to map python names to camel case, need to flatten fields from 'data'.
function flattenExperiment(serverExperiment): Experiment {
  const { data, ...sqlFields } = serverExperiment;
  // problem_id does NOT come along for this ride.
  const experiment = {
    ...data,
    createdAt: sqlFields.created_at,
    lastModified: sqlFields.last_modified,
    id: sqlFields.id,
  };
  return experiment;
}

function flattenExperiments(serverExperiments): Experiment[] {
  return serverExperiments.map((db_experiment) => {
    const { data, ...sqlFields } = db_experiment;
    // problem_id does NOT come along for this ride.
    const experiment = {
      ...data,
      createdAt: sqlFields.created_at,
      lastModified: sqlFields.last_modified,
      id: sqlFields.id,
    };
    return experiment;
  });
}

// Reformat from the ProjectResponseWithExperiments to the frontend Project.
function flattenProjects(serverProjects): Project[] {
  return serverProjects.map((db_project) => {
    return {
      id: db_project.id,
      name: db_project.name,
      createdAt: db_project.created_at,
      lastModified: db_project.last_modified,
      experiments: flattenExperiments(db_project.experiments),
    };
  });
}

const httpServerUrl = HTTP_SERVER;
class ServerDataSource implements ProjectDataSource {
  async loadProjects(): Project[] {
    try {
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

  saveProjects(projects: Project[]): void {
    throw new Error('We should not call saveProjects anymore');
  }

  async createProject(name: string): Project {
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
  async updateProject(project: Project): void {
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

  async deleteProject(projectId: string): void {
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

  async createExperiment(projectId: string, name: string): Experiment {
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
        return null;
      }
      const result = await response.json();
      return result;
    } catch (e) {
      console.error('Error creating experiment:', e);
    }
    return null;
  }

  async updateExperiment(projectId: string, experiment: Experiment): void {
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

  async deleteExperiment(projectId: string, experimentId: string): void {
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

  async setExperimentRunning(projectId: string, experimentId: string, isRunning: boolean): void {
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

  const createProject = async (name: string): Project => {
    const newProject = await dataSource.createProject(name);
    await loadProjects();
    return newProject;
  };

  const updateProject = async (project: Project) => {
    await dataSource.updateProject(project);
    await loadProjects();
  };

  const deleteProject = async (projectId: string) => {
    await dataSource.deleteProject(projectId);
    await loadProjects();
  };

  const createExperiment = async (projectId: string, name: string): Experiment => {
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
