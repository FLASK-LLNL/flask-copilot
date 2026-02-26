################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from pydantic import BaseModel, Field
from datetime import datetime
from typing import List, Optional, Any


class ExperimentBase(BaseModel):
    name: str


class ExperimentCreate(ExperimentBase):
    pass


class ExperimentUpdate(BaseModel):
    name: Optional[str] = None
    is_running: Optional[bool] = None

    # System state fields
    smiles: Optional[str] = None
    problem_type: Optional[str] = Field(None, alias="problemType")
    problem_name: Optional[str] = Field(None, alias="problemName")
    system_prompt: Optional[str] = Field(None, alias="systemPrompt")
    problem_prompt: Optional[str] = Field(None, alias="problemPrompt")
    property_type: Optional[str] = Field(None, alias="propertyType")
    custom_property_name: Optional[str] = Field(None, alias="customPropertyName")
    custom_property_desc: Optional[str] = Field(None, alias="customPropertyDesc")
    custom_property_ascending: Optional[bool] = Field(None, alias="customPropertyAscending")

    # Complex nested data (stored as JSON)
    tree_nodes: Optional[Any] = Field(None, alias="treeNodes")
    edges: Optional[Any] = None
    metrics_history: Optional[Any] = Field(None, alias="metricsHistory")
    visible_metrics: Optional[Any] = Field(None, alias="visibleMetrics")
    graph_state: Optional[Any] = Field(None, alias="graphState")
    auto_zoom: Optional[bool] = Field(None, alias="autoZoom")
    sidebar_state: Optional[Any] = Field(None, alias="sidebarState")

    # Experiment context
    experiment_context: Optional[str] = Field(None, alias="experimentContext")

    class Config:
        populate_by_name = True


class Experiment(ExperimentBase):
    id: str
    project_id: str = Field(..., alias="projectId")
    user: str
    created_at: datetime = Field(..., alias="createdAt")
    last_modified: datetime = Field(..., alias="lastModified")
    is_running: Optional[bool] = Field(None, alias="isRunning")

    # System state fields
    smiles: Optional[str] = None
    problem_type: Optional[str] = Field(None, alias="problemType")
    problem_name: Optional[str] = Field(None, alias="problemName")
    system_prompt: Optional[str] = Field(None, alias="systemPrompt")
    problem_prompt: Optional[str] = Field(None, alias="problemPrompt")
    property_type: Optional[str] = Field(None, alias="propertyType")
    custom_property_name: Optional[str] = Field(None, alias="customPropertyName")
    custom_property_desc: Optional[str] = Field(None, alias="customPropertyDesc")
    custom_property_ascending: Optional[bool] = Field(None, alias="customPropertyAscending")

    # Complex nested data
    tree_nodes: Optional[Any] = Field(None, alias="treeNodes")
    edges: Optional[Any] = None
    metrics_history: Optional[Any] = Field(None, alias="metricsHistory")
    visible_metrics: Optional[Any] = Field(None, alias="visibleMetrics")
    graph_state: Optional[Any] = Field(None, alias="graphState")
    auto_zoom: Optional[bool] = Field(None, alias="autoZoom")
    sidebar_state: Optional[Any] = Field(None, alias="sidebarState")

    # Experiment context
    experiment_context: Optional[str] = Field(None, alias="experimentContext")

    class Config:
        from_attributes = True
        populate_by_name = True


class ProjectBase(BaseModel):
    name: str


class ProjectCreate(ProjectBase):
    pass


class ProjectUpdate(ProjectBase):
    pass


class Project(ProjectBase):
    id: str
    user: str
    created_at: datetime = Field(..., alias="createdAt")
    last_modified: datetime = Field(..., alias="lastModified")
    experiments: List[Experiment] = []

    class Config:
        from_attributes = True
        populate_by_name = True
