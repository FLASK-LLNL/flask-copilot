################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
Shared session state service for WebSocket and API integrations.

This module contains the session persistence models and business logic,
without binding to FastAPI router/dependency wiring.
"""

from datetime import datetime
from typing import Any, List, Optional
import asyncio
import random
import time

from fastapi import HTTPException
from loguru import logger
from pydantic import BaseModel
from sqlalchemy import and_, select
from sqlalchemy.exc import OperationalError
from sqlalchemy.ext.asyncio import AsyncSession

from db_backend.database import models

# Max retries for concurrent modification errors (MariaDB error 1020)
_MAX_RETRIES = 3


def generate_id(prefix: str) -> str:
    """Generate unique ID."""
    timestamp = int(time.time() * 1000)
    random_str = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=9))
    return f"{prefix}_{timestamp}_{random_str}"


class SessionState(BaseModel):
    """The state to save/restore for a session."""

    projectId: Optional[str] = None
    projectName: Optional[str] = None
    experimentId: Optional[str] = None
    experimentName: Optional[str] = None

    smiles: Optional[str] = None
    problemType: Optional[str] = None
    systemPrompt: Optional[str] = None
    problemPrompt: Optional[str] = None
    promptsModified: Optional[bool] = None
    autoZoom: Optional[bool] = None
    treeNodes: Optional[Any] = None
    edges: Optional[Any] = None
    offset: Optional[Any] = None
    zoom: Optional[float] = None
    metricsHistory: Optional[Any] = None
    visibleMetrics: Optional[Any] = None
    isComputing: Optional[bool] = None
    serverSessionId: Optional[str] = None

    propertyType: Optional[str] = None
    customPropertyName: Optional[str] = None
    customPropertyDesc: Optional[str] = None
    customPropertyAscending: Optional[bool] = None

    sidebarMessages: Optional[Any] = None
    sidebarSourceFilterOpen: Optional[bool] = None
    sidebarVisibleSources: Optional[Any] = None


class SessionSaveRequest(BaseModel):
    """Request body for saving a session."""

    sessionId: Optional[str] = None
    name: Optional[str] = None
    state: SessionState


class SessionResponse(BaseModel):
    """Response for session operations."""

    sessionId: str
    name: str
    createdAt: datetime
    lastModified: datetime
    isRunning: Optional[bool] = None
    state: SessionState

    class Config:
        from_attributes = True


class SessionListItem(BaseModel):
    """Summary info for listing sessions."""

    sessionId: str
    name: str
    smiles: Optional[str] = None
    problemType: Optional[str] = None
    createdAt: datetime
    lastModified: datetime
    isRunning: Optional[bool] = None

    class Config:
        from_attributes = True


async def save_session(
    request: SessionSaveRequest,
    user: str,
    db: AsyncSession,
) -> SessionResponse:
    """Save or update a session."""
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")

    state = request.state

    if request.sessionId:
        for attempt in range(_MAX_RETRIES):
            try:
                result = await db.execute(
                    select(models.Experiment).where(
                        and_(
                            models.Experiment.id == request.sessionId,
                            models.Experiment.user == user,
                        )
                    )
                )
                experiment = result.scalar_one_or_none()

                if experiment:
                    incoming_nodes = state.treeNodes or []
                    incoming_edges = state.edges or []
                    incoming_smiles = state.smiles or ""

                    if incoming_smiles or not experiment.smiles:
                        experiment.smiles = state.smiles
                    if incoming_nodes or not experiment.tree_nodes:
                        experiment.tree_nodes = state.treeNodes
                    if incoming_edges or not experiment.edges:
                        experiment.edges = state.edges

                    experiment.problem_type = state.problemType
                    experiment.system_prompt = state.systemPrompt
                    experiment.problem_prompt = state.problemPrompt
                    experiment.auto_zoom = state.autoZoom
                    experiment.metrics_history = state.metricsHistory
                    experiment.visible_metrics = state.visibleMetrics
                    experiment.is_running = state.isComputing

                    incoming_sidebar_msgs = state.sidebarMessages or []
                    existing_sidebar_msgs = []
                    if experiment.sidebar_state and isinstance(
                        experiment.sidebar_state, dict
                    ):
                        existing_sidebar_msgs = experiment.sidebar_state.get(
                            "messages", []
                        )

                    if len(incoming_sidebar_msgs) >= len(existing_sidebar_msgs):
                        experiment.sidebar_state = {
                            "messages": state.sidebarMessages or [],
                            "sourceFilterOpen": state.sidebarSourceFilterOpen or False,
                            "visibleSources": state.sidebarVisibleSources or {},
                        }
                    else:
                        logger.info(
                            f"save_session: keeping existing {len(existing_sidebar_msgs)} "
                            f"sidebar messages for {request.sessionId} "
                            f"(incoming had only {len(incoming_sidebar_msgs)})"
                        )

                    existing_gs = experiment.graph_state or {}
                    existing_gs_msgs = existing_gs.get("sidebarMessages") or []
                    experiment.graph_state = {
                        "offset": state.offset,
                        "zoom": state.zoom,
                        "promptsModified": state.promptsModified,
                        "serverSessionId": state.serverSessionId,
                        "propertyType": state.propertyType,
                        "customPropertyName": state.customPropertyName,
                        "customPropertyDesc": state.customPropertyDesc,
                        "customPropertyAscending": state.customPropertyAscending,
                        "sidebarMessages": (
                            incoming_sidebar_msgs
                            if len(incoming_sidebar_msgs) >= len(existing_gs_msgs)
                            else existing_gs_msgs
                        ),
                        "sidebarSourceFilterOpen": state.sidebarSourceFilterOpen,
                        "sidebarVisibleSources": state.sidebarVisibleSources,
                    }
                    experiment.last_modified = datetime.utcnow()

                    if request.name:
                        experiment.name = request.name

                    await db.commit()
                    await db.refresh(experiment)

                    return SessionResponse(
                        sessionId=experiment.id,
                        name=experiment.name,
                        createdAt=experiment.created_at,
                        lastModified=experiment.last_modified,
                        isRunning=experiment.is_running,
                        state=_experiment_to_state(experiment),
                    )
                break
            except OperationalError as exc:
                if "1020" in str(exc) and attempt < _MAX_RETRIES - 1:
                    logger.warning(
                        f"Concurrent modification on experiment {request.sessionId}, retry {attempt + 1}/{_MAX_RETRIES}"
                    )
                    await db.rollback()
                    await asyncio.sleep(0.05 * (attempt + 1))
                    continue
                raise

    project = None
    if state.projectId:
        result = await db.execute(
            select(models.Project).where(
                and_(
                    models.Project.id == state.projectId,
                    models.Project.user == user,
                )
            )
        )
        project = result.scalar_one_or_none()

    if project is None:
        raise HTTPException(
            status_code=422,
            detail="No valid project specified. Create a project first.",
        )

    session_name = (
        request.name or f"Session {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    if state.smiles:
        session_name = (
            f"{state.smiles[:20]}..." if len(state.smiles) > 20 else state.smiles
        )

    new_experiment = models.Experiment(
        id=generate_id("session"),
        project_id=project.id,
        user=user,
        name=session_name,
        created_at=datetime.utcnow(),
        last_modified=datetime.utcnow(),
        smiles=state.smiles,
        problem_type=state.problemType,
        system_prompt=state.systemPrompt,
        problem_prompt=state.problemPrompt,
        auto_zoom=state.autoZoom,
        tree_nodes=state.treeNodes,
        edges=state.edges,
        metrics_history=state.metricsHistory,
        visible_metrics=state.visibleMetrics,
        is_running=state.isComputing,
        graph_state={
            "offset": state.offset,
            "zoom": state.zoom,
            "promptsModified": state.promptsModified,
            "serverSessionId": state.serverSessionId,
            "propertyType": state.propertyType,
            "customPropertyName": state.customPropertyName,
            "customPropertyDesc": state.customPropertyDesc,
            "customPropertyAscending": state.customPropertyAscending,
            "sidebarMessages": state.sidebarMessages,
            "sidebarSourceFilterOpen": state.sidebarSourceFilterOpen,
            "sidebarVisibleSources": state.sidebarVisibleSources,
        },
        sidebar_state={
            "messages": state.sidebarMessages or [],
            "sourceFilterOpen": state.sidebarSourceFilterOpen or False,
            "visibleSources": state.sidebarVisibleSources or {},
        },
    )

    project.last_modified = datetime.utcnow()

    db.add(new_experiment)
    await db.commit()
    await db.refresh(new_experiment)

    return SessionResponse(
        sessionId=new_experiment.id,
        name=new_experiment.name,
        createdAt=new_experiment.created_at,
        lastModified=new_experiment.last_modified,
        isRunning=new_experiment.is_running,
        state=_experiment_to_state(new_experiment),
    )


async def get_latest_session(
    user: str,
    db: AsyncSession,
) -> Optional[SessionResponse]:
    """Get the most recently modified session for the user."""
    if db is None:
        return None

    result = await db.execute(
        select(models.Experiment)
        .where(models.Experiment.user == user)
        .order_by(models.Experiment.last_modified.desc())
        .limit(1)
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        return None

    return SessionResponse(
        sessionId=experiment.id,
        name=experiment.name,
        createdAt=experiment.created_at,
        lastModified=experiment.last_modified,
        isRunning=experiment.is_running,
        state=_experiment_to_state(experiment),
    )


async def get_session(
    session_id: str,
    user: str,
    db: AsyncSession,
) -> SessionResponse:
    """Get a specific session by ID."""
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")

    result = await db.execute(
        select(models.Experiment).where(
            and_(models.Experiment.id == session_id, models.Experiment.user == user)
        )
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        raise HTTPException(status_code=404, detail="Session not found")

    return SessionResponse(
        sessionId=experiment.id,
        name=experiment.name,
        createdAt=experiment.created_at,
        lastModified=experiment.last_modified,
        isRunning=experiment.is_running,
        state=_experiment_to_state(experiment),
    )


async def list_sessions(
    limit: int,
    user: str,
    db: AsyncSession,
) -> List[SessionListItem]:
    """List all sessions for the user, ordered by last modified."""
    if db is None:
        return []

    result = await db.execute(
        select(models.Experiment)
        .where(models.Experiment.user == user)
        .order_by(models.Experiment.last_modified.desc())
        .limit(limit)
    )
    experiments = result.scalars().all()

    return [
        SessionListItem(
            sessionId=exp.id,
            name=exp.name,
            smiles=exp.smiles,
            problemType=exp.problem_type,
            createdAt=exp.created_at,
            lastModified=exp.last_modified,
            isRunning=exp.is_running,
        )
        for exp in experiments
    ]


async def delete_session(
    session_id: str,
    user: str,
    db: AsyncSession,
) -> dict:
    """Delete a session."""
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")

    result = await db.execute(
        select(models.Experiment).where(
            and_(models.Experiment.id == session_id, models.Experiment.user == user)
        )
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        raise HTTPException(status_code=404, detail="Session not found")

    await db.delete(experiment)
    await db.commit()

    return {"success": True}


def _experiment_to_state(experiment: models.Experiment) -> SessionState:
    graph_state = experiment.graph_state or {}

    sidebar = experiment.sidebar_state
    if sidebar:
        sidebar_messages = sidebar.get("messages")
        sidebar_filter_open = sidebar.get("sourceFilterOpen", False)
        sidebar_visible = sidebar.get("visibleSources", {})
    else:
        sidebar_messages = graph_state.get("sidebarMessages")
        sidebar_filter_open = graph_state.get("sidebarSourceFilterOpen")
        sidebar_visible = graph_state.get("sidebarVisibleSources")

    return SessionState(
        projectId=experiment.project_id,
        experimentId=experiment.id,
        experimentName=experiment.name,
        smiles=experiment.smiles,
        problemType=experiment.problem_type,
        systemPrompt=experiment.system_prompt,
        problemPrompt=experiment.problem_prompt,
        promptsModified=graph_state.get("promptsModified"),
        autoZoom=experiment.auto_zoom,
        treeNodes=experiment.tree_nodes,
        edges=experiment.edges,
        offset=graph_state.get("offset"),
        zoom=graph_state.get("zoom"),
        metricsHistory=experiment.metrics_history,
        visibleMetrics=experiment.visible_metrics,
        isComputing=experiment.is_running,
        serverSessionId=graph_state.get("serverSessionId"),
        propertyType=graph_state.get("propertyType"),
        customPropertyName=graph_state.get("customPropertyName"),
        customPropertyDesc=graph_state.get("customPropertyDesc"),
        customPropertyAscending=graph_state.get("customPropertyAscending"),
        sidebarMessages=sidebar_messages,
        sidebarSourceFilterOpen=sidebar_filter_open,
        sidebarVisibleSources=sidebar_visible,
    )
