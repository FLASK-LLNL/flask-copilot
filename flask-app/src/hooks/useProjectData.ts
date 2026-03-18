import { useState, useEffect, useRef } from 'react';
import { Project, Experiment } from '../types';

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
  const dataSource: ProjectDataSource = new LocalStorageDataSource();

  useEffect(() => {
    loadProjects();
  }, []);

  const loadProjects = () => {
    setLoading(true);
    try {
      const data = dataSource.loadProjects();
      setProjects(data);
      projectsRef.current = data;
    } catch (error) {
      console.error('Error loading projects:', error);
    } finally {
      setLoading(false);
    }
  };

  const createProject = (name: string): Project => {
    const newProject = dataSource.createProject(name);
    loadProjects();
    return newProject;
  };

  const updateProject = (project: Project) => {
    dataSource.updateProject(project);
    loadProjects();
  };

  const deleteProject = (projectId: string) => {
    dataSource.deleteProject(projectId);
    loadProjects();
  };

  const createExperiment = (projectId: string, name: string): Experiment => {
    const newExperiment = dataSource.createExperiment(projectId, name);
    loadProjects();
    return newExperiment;
  };

  const updateExperiment = (projectId: string, experiment: Experiment) => {
    dataSource.updateExperiment(projectId, experiment);
    loadProjects();
  };

  const deleteExperiment = (projectId: string, experimentId: string) => {
    dataSource.deleteExperiment(projectId, experimentId);
    loadProjects();
  };

  const setExperimentRunning = (projectId: string, experimentId: string, isRunning: boolean) => {
    dataSource.setExperimentRunning(projectId, experimentId, isRunning);
    loadProjects();
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
