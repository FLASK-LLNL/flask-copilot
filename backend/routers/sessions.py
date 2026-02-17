################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
Serves routes for auto-save session management (``/api/sessions/*``).
This provides a simpler interface for automatic session persistence without
requiring explicit project/experiment management.
"""
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.exc import OperationalError
from sqlalchemy import select, and_
from pydantic import BaseModel, Field
from typing import Optional, Any, List
from datetime import datetime
from loguru import logger
import time
import random

# Max retries for concurrent modification errors (MariaDB error 1020)
_MAX_RETRIES = 3

from backend.database.engine import get_db
from backend.auth import get_forwarded_user
from backend.database import models

router = APIRouter(prefix="/api/sessions", tags=["sessions"])


def generate_id(prefix: str) -> str:
    """Generate unique ID"""
    timestamp = int(time.time() * 1000)
    random_str = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=9))
    return f"{prefix}_{timestamp}_{random_str}"


# Pydantic models for session API
class SessionState(BaseModel):
    """The state to save/restore for a session"""
    # Project/Experiment identification
    projectId: Optional[str] = None
    projectName: Optional[str] = None
    experimentId: Optional[str] = None
    experimentName: Optional[str] = None
    
    # Core experiment state
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
    
    # Property optimization
    propertyType: Optional[str] = None
    customPropertyName: Optional[str] = None
    customPropertyDesc: Optional[str] = None
    customPropertyAscending: Optional[bool] = None
    
    # Sidebar state
    sidebarMessages: Optional[Any] = None
    sidebarSourceFilterOpen: Optional[bool] = None
    sidebarVisibleSources: Optional[Any] = None


class SessionSaveRequest(BaseModel):
    """Request body for saving a session"""
    sessionId: Optional[str] = None  # If provided, update existing session
    name: Optional[str] = None  # Optional name for the session
    state: SessionState


class SessionResponse(BaseModel):
    """Response for session operations"""
    sessionId: str
    name: str
    createdAt: datetime
    lastModified: datetime
    isRunning: Optional[bool] = None
    state: SessionState

    class Config:
        from_attributes = True


class SessionListItem(BaseModel):
    """Summary info for listing sessions"""
    sessionId: str
    name: str
    smiles: Optional[str] = None
    problemType: Optional[str] = None
    createdAt: datetime
    lastModified: datetime
    isRunning: Optional[bool] = None

    class Config:
        from_attributes = True


# Default project for auto-save sessions
DEFAULT_PROJECT_NAME = "Auto-Save Sessions"


async def get_or_create_default_project(user: str, db: AsyncSession) -> models.Project:
    """Get or create the default project for auto-save sessions"""
    result = await db.execute(
        select(models.Project).where(
            and_(
                models.Project.user == user,
                models.Project.name == DEFAULT_PROJECT_NAME
            )
        ).order_by(models.Project.created_at.asc())
        .limit(1)
    )
    project = result.scalar_one_or_none()
    
    if not project:
        project = models.Project(
            id=generate_id("project"),
            user=user,
            name=DEFAULT_PROJECT_NAME,
            created_at=datetime.utcnow(),
            last_modified=datetime.utcnow(),
        )
        db.add(project)
        await db.flush()  # Get the ID without committing
    
    return project


@router.post("/save", response_model=SessionResponse)
async def save_session(
    request: SessionSaveRequest,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Save or update a session. If sessionId is provided, updates existing session.
    Otherwise creates a new session.
    """
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")
    
    state = request.state
    
    # If sessionId provided, try to update existing
    if request.sessionId:
        for attempt in range(_MAX_RETRIES):
            try:
                result = await db.execute(
                    select(models.Experiment).where(
                        and_(
                            models.Experiment.id == request.sessionId,
                            models.Experiment.user == user
                        )
                    )
                )
                experiment = result.scalar_one_or_none()
                
                if experiment:
                    # Update existing session
                    experiment.smiles = state.smiles
                    experiment.problem_type = state.problemType
                    experiment.system_prompt = state.systemPrompt
                    experiment.problem_prompt = state.problemPrompt
                    experiment.auto_zoom = state.autoZoom
                    experiment.tree_nodes = state.treeNodes
                    experiment.edges = state.edges
                    experiment.metrics_history = state.metricsHistory
                    experiment.visible_metrics = state.visibleMetrics
                    experiment.is_running = state.isComputing
                    experiment.graph_state = {
                        "offset": state.offset,
                        "zoom": state.zoom,
                        "promptsModified": state.promptsModified,
                        "serverSessionId": state.serverSessionId,
                        # Property optimization
                        "propertyType": state.propertyType,
                        "customPropertyName": state.customPropertyName,
                        "customPropertyDesc": state.customPropertyDesc,
                        "customPropertyAscending": state.customPropertyAscending,
                        # Sidebar state
                        "sidebarMessages": state.sidebarMessages,
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
                        state=_experiment_to_state(experiment)
                    )
                else:
                    break  # Experiment not found, fall through to create
            except OperationalError as e:
                if "1020" in str(e) and attempt < _MAX_RETRIES - 1:
                    logger.warning(f"Concurrent modification on experiment {request.sessionId}, retry {attempt + 1}/{_MAX_RETRIES}")
                    await db.rollback()
                    continue
                raise
    
    # Create new session.
    # The frontend must supply a projectId that maps to an existing project.
    # If none is provided (or it doesn't exist), reject the save so that
    # stale checkpoint/unload saves cannot resurrect data after a clear.
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
    
    # Generate session name if not provided
    session_name = request.name or f"Session {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    if state.smiles:
        session_name = f"{state.smiles[:20]}..." if len(state.smiles) > 20 else state.smiles
    
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
            # Property optimization
            "propertyType": state.propertyType,
            "customPropertyName": state.customPropertyName,
            "customPropertyDesc": state.customPropertyDesc,
            "customPropertyAscending": state.customPropertyAscending,
            # Sidebar state
            "sidebarMessages": state.sidebarMessages,
            "sidebarSourceFilterOpen": state.sidebarSourceFilterOpen,
            "sidebarVisibleSources": state.sidebarVisibleSources,
        }
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
        state=_experiment_to_state(new_experiment)
    )


@router.get("/latest", response_model=Optional[SessionResponse])
async def get_latest_session(
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db)
):
    """Get the most recently modified session for the user"""
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
        state=_experiment_to_state(experiment)
    )


@router.get("/{session_id}", response_model=SessionResponse)
async def get_session(
    session_id: str,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db)
):
    """Get a specific session by ID"""
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")
    
    result = await db.execute(
        select(models.Experiment).where(
            and_(
                models.Experiment.id == session_id,
                models.Experiment.user == user
            )
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
        state=_experiment_to_state(experiment)
    )


@router.get("/", response_model=List[SessionListItem])
async def list_sessions(
    limit: int = 20,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db)
):
    """List all sessions for the user, ordered by last modified"""
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
            isRunning=exp.is_running
        )
        for exp in experiments
    ]


@router.delete("/{session_id}")
async def delete_session(
    session_id: str,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db)
):
    """Delete a session"""
    if db is None:
        raise HTTPException(status_code=503, detail="Database not available")
    
    result = await db.execute(
        select(models.Experiment).where(
            and_(
                models.Experiment.id == session_id,
                models.Experiment.user == user
            )
        )
    )
    experiment = result.scalar_one_or_none()
    
    if not experiment:
        raise HTTPException(status_code=404, detail="Session not found")
    
    await db.delete(experiment)
    await db.commit()
    
    return {"success": True}


def _experiment_to_state(experiment: models.Experiment) -> SessionState:
    """Convert an Experiment model to SessionState"""
    graph_state = experiment.graph_state or {}
    
    return SessionState(
        # Project/Experiment identification
        projectId=experiment.project_id,
        experimentId=experiment.id,
        experimentName=experiment.name,
        
        # Core experiment state
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
        
        # Property optimization
        propertyType=graph_state.get("propertyType"),
        customPropertyName=graph_state.get("customPropertyName"),
        customPropertyDesc=graph_state.get("customPropertyDesc"),
        customPropertyAscending=graph_state.get("customPropertyAscending"),
        
        # Sidebar state
        sidebarMessages=graph_state.get("sidebarMessages"),
        sidebarSourceFilterOpen=graph_state.get("sidebarSourceFilterOpen"),
        sidebarVisibleSources=graph_state.get("sidebarVisibleSources"),
    )
