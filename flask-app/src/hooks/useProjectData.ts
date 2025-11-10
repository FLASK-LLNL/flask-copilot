import { useState, useEffect } from 'react';
import { Project, Experiment } from '../types';

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
  setExperimentRunning: (projectId: string, experimentId: string, isRunning: boolean) => Promise<void>;
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
      experiments: []
    };
    
    const projects = await this.loadProjects();
    projects.push(newProject);
    await this.saveProjects(projects);
    
    return newProject;
  }

  async updateProject(project: Project): Promise<void> {
    const projects = await this.loadProjects();
    const index = projects.findIndex(p => p.id === project.id);
    if (index !== -1) {
      projects[index] = { ...project, lastModified: new Date().toISOString() };
      await this.saveProjects(projects);
    }
  }

  async deleteProject(projectId: string): Promise<void> {
    const projects = await this.loadProjects();
    const filtered = projects.filter(p => p.id !== projectId);
    await this.saveProjects(filtered);
  }

  async createExperiment(projectId: string, name: string): Promise<Experiment> {
    const newExperiment: Experiment = {
      id: `exp_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString()
    };

    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      project.experiments.push(newExperiment);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }

    return newExperiment;
  }

  async updateExperiment(projectId: string, experiment: Experiment): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      const expIndex = project.experiments.findIndex(e => e.id === experiment.id);
      if (expIndex !== -1) {
        project.experiments[expIndex] = { ...experiment, lastModified: new Date().toISOString() };
        project.lastModified = new Date().toISOString();
        await this.saveProjects(projects);
      }
    }
  }

  async deleteExperiment(projectId: string, experimentId: string): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      project.experiments = project.experiments.filter(e => e.id !== experimentId);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }
  }

  async setExperimentRunning(projectId: string, experimentId: string, isRunning: boolean): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      const experiment = project.experiments.find(e => e.id === experimentId);
      if (experiment) {
        experiment.isRunning = isRunning;
        experiment.lastModified = new Date().toISOString();
        project.lastModified = new Date().toISOString();
        await this.saveProjects(projects);
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

// Hook for managing project data
export const useProjectData = () => {
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(true);
  
  // TODO(later): Swap this to use ServerDataSource
  const dataSource: ProjectDataSource = new LocalStorageDataSource();

  useEffect(() => {
    loadProjects();
  }, []);

  const loadProjects = async () => {
    setLoading(true);
    try {
      const data = await dataSource.loadProjects();
      setProjects(data);
    } catch (error) {
      console.error('Error loading projects:', error);
    } finally {
      setLoading(false);
    }
  };

  const createProject = async (name: string): Promise<Project> => {
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

  const createExperiment = async (projectId: string, name: string): Promise<Experiment> => {
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

  const setExperimentRunning = async (projectId: string, experimentId: string, isRunning: boolean) => {
    await dataSource.setExperimentRunning(projectId, experimentId, isRunning);
    await loadProjects();
  };

  return {
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    setExperimentRunning,
    refreshProjects: loadProjects
  };
};
