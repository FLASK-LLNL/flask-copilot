import { useState, useEffect } from 'react';
import { Experiment, Task } from '../types';

const STORAGE_KEY = 'flask_copilot_experiments';

// This interface defines the contract for data sources (local storage, database, etc.)
interface ExperimentDataSource {
  loadExperiments: () => Promise<Experiment[]>;
  saveExperiments: (experiments: Experiment[]) => Promise<void>;
  createExperiment: (name: string) => Promise<Experiment>;
  updateExperiment: (experiment: Experiment) => Promise<void>;
  deleteExperiment: (experimentId: string) => Promise<void>;
  createTask: (experimentId: string, name: string) => Promise<Task>;
  updateTask: (experimentId: string, task: Task) => Promise<void>;
  deleteTask: (experimentId: string, taskId: string) => Promise<void>;
  setTaskRunning: (experimentId: string, taskId: string, isRunning: boolean) => Promise<void>;
}

// LocalStorage implementation
class LocalStorageDataSource implements ExperimentDataSource {
  async loadExperiments(): Promise<Experiment[]> {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (!stored) return [];
    try {
      return JSON.parse(stored);
    } catch (e) {
      console.error('Error loading experiments from localStorage:', e);
      return [];
    }
  }

  async saveExperiments(experiments: Experiment[]): Promise<void> {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(experiments));
  }

  async createExperiment(name: string): Promise<Experiment> {
    const newExperiment: Experiment = {
      id: `experiment_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString(),
      tasks: []
    };
    
    const experiments = await this.loadExperiments();
    experiments.push(newExperiment);
    await this.saveExperiments(experiments);
    
    return newExperiment;
  }

  async updateExperiment(experiment: Experiment): Promise<void> {
    const experiments = await this.loadExperiments();
    const index = experiments.findIndex(p => p.id === experiment.id);
    if (index !== -1) {
      experiments[index] = { ...experiment, lastModified: new Date().toISOString() };
      await this.saveExperiments(experiments);
    }
  }

  async deleteExperiment(experimentId: string): Promise<void> {
    const experiments = await this.loadExperiments();
    const filtered = experiments.filter(p => p.id !== experimentId);
    await this.saveExperiments(filtered);
  }

  async createTask(experimentId: string, name: string): Promise<Task> {
    const newTask: Task = {
      id: `exp_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      name,
      createdAt: new Date().toISOString(),
      lastModified: new Date().toISOString()
    };

    const experiments = await this.loadExperiments();
    const experiment = experiments.find(p => p.id === experimentId);
    if (experiment) {
      experiment.tasks.push(newTask);
      experiment.lastModified = new Date().toISOString();
      await this.saveExperiments(experiments);
    }

    return newTask;
  }

  async updateTask(experimentId: string, task: Task): Promise<void> {
    const experiments = await this.loadExperiments();
    const experiment = experiments.find(p => p.id === experimentId);
    if (experiment) {
      const expIndex = experiment.tasks.findIndex(e => e.id === task.id);
      if (expIndex !== -1) {
        experiment.tasks[expIndex] = { ...task, lastModified: new Date().toISOString() };
        experiment.lastModified = new Date().toISOString();
        await this.saveExperiments(experiments);
      }
    }
  }

  async deleteTask(experimentId: string, taskId: string): Promise<void> {
    const experiments = await this.loadExperiments();
    const experiment = experiments.find(p => p.id === experimentId);
    if (experiment) {
      experiment.tasks = experiment.tasks.filter(e => e.id !== taskId);
      experiment.lastModified = new Date().toISOString();
      await this.saveExperiments(experiments);
    }
  }

  async setTaskRunning(experimentId: string, taskId: string, isRunning: boolean): Promise<void> {
    const experiments = await this.loadExperiments();
    const experiment = experiments.find(p => p.id === experimentId);
    if (experiment) {
      const task = experiment.tasks.find(e => e.id === taskId);
      if (task) {
        task.isRunning = isRunning;
        task.lastModified = new Date().toISOString();
        experiment.lastModified = new Date().toISOString();
        await this.saveExperiments(experiments);
      }
    }
  }
}

// TODO: Future server implementation
// class ServerDataSource implements ExperimentDataSource {
//   async loadExperiments(): Promise<Experiment[]> {
//     const response = await fetch('/api/experiments');
//     return response.json();
//   }
//   // ...
// }

// Hook for managing experiment data
export const useExperimentData = () => {
  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [loading, setLoading] = useState(true);
  
  // TODO(later): Swap this to use ServerDataSource
  const dataSource: ExperimentDataSource = new LocalStorageDataSource();

  useEffect(() => {
    loadExperiments();
  }, []);

  const loadExperiments = async () => {
    setLoading(true);
    try {
      const data = await dataSource.loadExperiments();
      setExperiments(data);
    } catch (error) {
      console.error('Error loading experiments:', error);
    } finally {
      setLoading(false);
    }
  };

  const createExperiment = async (name: string): Promise<Experiment> => {
    const newExperiment = await dataSource.createExperiment(name);
    await loadExperiments();
    return newExperiment;
  };

  const updateExperiment = async (experiment: Experiment) => {
    await dataSource.updateExperiment(experiment);
    await loadExperiments();
  };

  const deleteExperiment = async (experimentId: string) => {
    await dataSource.deleteExperiment(experimentId);
    await loadExperiments();
  };

  const createTask = async (experimentId: string, name: string): Promise<Task> => {
    const newTask = await dataSource.createTask(experimentId, name);
    await loadExperiments();
    return newTask;
  };

  const updateTask = async (experimentId: string, task: Task) => {
    await dataSource.updateTask(experimentId, task);
    await loadExperiments();
  };

  const deleteTask = async (experimentId: string, taskId: string) => {
    await dataSource.deleteTask(experimentId, taskId);
    await loadExperiments();
  };

  const setTaskRunning = async (experimentId: string, taskId: string, isRunning: boolean) => {
    await dataSource.setTaskRunning(experimentId, taskId, isRunning);
    await loadExperiments();
  };

  return {
    experiments,
    loading,
    createExperiment,
    updateExperiment,
    deleteExperiment,
    createTask,
    updateTask,
    deleteTask,
    setTaskRunning,
    refreshExperiments: loadExperiments
  };
};
