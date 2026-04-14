################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from __future__ import annotations

from datetime import datetime
from typing import List, Optional, TypeAlias

from pydantic import BaseModel, Field

JsonObject: TypeAlias = dict[str, object]
JsonArray: TypeAlias = list[JsonObject]
JsonValue: TypeAlias = object


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
    custom_property_ascending: Optional[bool] = Field(
        None, alias="customPropertyAscending"
    )

    # Complex nested data (stored as JSON)
    tree_nodes: Optional[JsonArray] = Field(None, alias="treeNodes")
    edges: Optional[JsonArray] = None
    metrics_history: Optional[JsonArray] = Field(None, alias="metricsHistory")
    visible_metrics: Optional[dict[str, bool]] = Field(None, alias="visibleMetrics")
    graph_state: Optional[JsonObject] = Field(None, alias="graphState")
    auto_zoom: Optional[bool] = Field(None, alias="autoZoom")
    sidebar_state: Optional[JsonObject] = Field(None, alias="sidebarState")

    # Experiment context
    experiment_context: Optional[JsonValue] = Field(None, alias="experimentContext")

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
    custom_property_ascending: Optional[bool] = Field(
        None, alias="customPropertyAscending"
    )

    # Complex nested data
    tree_nodes: Optional[JsonArray] = Field(None, alias="treeNodes")
    edges: Optional[JsonArray] = None
    metrics_history: Optional[JsonArray] = Field(None, alias="metricsHistory")
    visible_metrics: Optional[dict[str, bool]] = Field(None, alias="visibleMetrics")
    graph_state: Optional[JsonObject] = Field(None, alias="graphState")
    auto_zoom: Optional[bool] = Field(None, alias="autoZoom")
    sidebar_state: Optional[JsonObject] = Field(None, alias="sidebarState")

    # Experiment context
    experiment_context: Optional[JsonValue] = Field(None, alias="experimentContext")

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


class ExperimentData(BaseModel):
    id: str
    name: str
    createdAt: str
    lastModified: str
    isRunning: Optional[bool] = None
    smiles: Optional[str] = None
    problemType: Optional[str] = None
    problemName: Optional[str] = None
    systemPrompt: Optional[str] = None
    problemPrompt: Optional[str] = None
    propertyType: Optional[str] = None
    customPropertyName: Optional[str] = None
    customPropertyDesc: Optional[str] = None
    customPropertyAscending: Optional[bool] = None
    customization: Optional[JsonObject] = None
    treeNodes: Optional[JsonArray] = None
    edges: Optional[JsonArray] = None
    metricsHistory: Optional[JsonArray] = None
    visibleMetrics: Optional[dict[str, bool]] = None
    graphState: Optional[JsonObject] = None
    autoZoom: Optional[bool] = None
    sidebarState: Optional[JsonObject] = None
    experimentContext: Optional[JsonValue] = None


class ProjectData(BaseModel):
    id: str
    name: str
    createdAt: str
    lastModified: str
    experiments: List[ExperimentData]


class CreateProjectRequest(BaseModel):
    name: str


class CreateExperimentRequest(BaseModel):
    name: str


class UpdateExperimentRequest(BaseModel):
    experiment: ExperimentData
