import { useState, useEffect } from 'react';
import { Project, Task } from '../types';

const STORAGE_KEY = 'flask_copilot_projects';

// This interface defines the contract for data sources (local storage, database, etc.)
interface ProjectDataSource {
  loadProjects: () => Promise<Project[]>;
  saveProjects: (projects: Project[]) => Promise<void>;
  createProject: (name: string) => Promise<Project>;
  updateProject: (project: Project) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  createTask: (projectId: string, name: string) => Promise<Task>;
  updateTask: (projectId: string, task: Task) => Promise<void>;
  deleteTask: (projectId: string, taskId: string) => Promise<void>;
  setTaskRunning: (projectId: string, taskId: string, isRunning: boolean) => Promise<void>;
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
      tasks: []
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

  async createTask(projectId: string, name: string): Promise<Task> {
    const newTask: Task = {
      id: `exp_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString()
    };

    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      project.tasks.push(newTask);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }

    return newTask;
  }

  async updateTask(projectId: string, task: Task): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      const expIndex = project.tasks.findIndex(e => e.id === task.id);
      if (expIndex !== -1) {
        project.tasks[expIndex] = { ...task, lastModified: new Date().toISOString() };
        project.lastModified = new Date().toISOString();
        await this.saveProjects(projects);
      }
    }
  }

  async deleteTask(projectId: string, taskId: string): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      project.tasks = project.tasks.filter(e => e.id !== taskId);
      project.lastModified = new Date().toISOString();
      await this.saveProjects(projects);
    }
  }

  async setTaskRunning(projectId: string, taskId: string, isRunning: boolean): Promise<void> {
    const projects = await this.loadProjects();
    const project = projects.find(p => p.id === projectId);
    if (project) {
      const task = project.tasks.find(e => e.id === taskId);
      if (task) {
        task.isRunning = isRunning;
        task.lastModified = new Date().toISOString();
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

  const createTask = async (projectId: string, name: string): Promise<Task> => {
    const newTask = await dataSource.createTask(projectId, name);
    await loadProjects();
    return newTask;
  };

  const updateTask = async (projectId: string, task: Task) => {
    await dataSource.updateTask(projectId, task);
    await loadProjects();
  };

  const deleteTask = async (projectId: string, taskId: string) => {
    await dataSource.deleteTask(projectId, taskId);
    await loadProjects();
  };

  const setTaskRunning = async (projectId: string, taskId: string, isRunning: boolean) => {
    await dataSource.setTaskRunning(projectId, taskId, isRunning);
    await loadProjects();
  };

  return {
    projects,
    loading,
    createProject,
    updateProject,
    deleteProject,
    createTask,
    updateTask,
    deleteTask,
    setTaskRunning,
    refreshProjects: loadProjects
  };
};
