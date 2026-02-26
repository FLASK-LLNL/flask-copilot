################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
API routes for project and experiment management (``/api/projects/*``).
Enables multi-browser session sharing via MariaDB storage.
"""
import asyncio
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.exc import OperationalError
from sqlalchemy import select, delete
from sqlalchemy.orm import selectinload
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime
from loguru import logger
import time
import random

# Max retries for concurrent modification errors (MariaDB error 1020)
_MAX_RETRIES = 3

from db_backend.database.engine import get_db
from db_backend.auth import get_forwarded_user
from db_backend.database import models

router = APIRouter(prefix="/api/projects", tags=["projects"])


def generate_id(prefix: str) -> str:
    """Generate unique ID"""
    timestamp = int(time.time() * 1000)
    random_str = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=9))
    return f"{prefix}_{timestamp}_{random_str}"


# Pydantic models
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
    customization: Optional[Dict[str, Any]] = None
    treeNodes: Optional[List[Dict[str, Any]]] = None
    edges: Optional[List[Dict[str, Any]]] = None
    metricsHistory: Optional[List[Dict[str, Any]]] = None
    visibleMetrics: Optional[Dict[str, bool]] = None
    graphState: Optional[Dict[str, Any]] = None
    autoZoom: Optional[bool] = None
    sidebarState: Optional[Dict[str, Any]] = None
    experimentContext: Optional[str] = None


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


def experiment_to_dict(exp: models.Experiment) -> Dict[str, Any]:
    """Convert Experiment model to dict matching frontend format"""
    # Handle customization field - for now just return None since it's not in the model
    customization = None

    # Prefer the dedicated sidebar_state column.  Fall back to
    # graph_state (where the session/checkpoint path historically
    # stored sidebar data) so that older rows still surface
    # reasoning messages.
    sidebar = exp.sidebar_state
    if not sidebar and exp.graph_state:
        gs = exp.graph_state
        msgs = gs.get("sidebarMessages")
        if msgs:
            sidebar = {
                "messages": msgs,
                "sourceFilterOpen": gs.get("sidebarSourceFilterOpen", False),
                "visibleSources": gs.get("sidebarVisibleSources", {}),
            }

    return {
        "id": exp.id,
        "name": exp.name,
        "createdAt": exp.created_at.isoformat(),
        "lastModified": exp.last_modified.isoformat(),
        "isRunning": exp.is_running,
        "smiles": exp.smiles,
        "problemType": exp.problem_type,
        "problemName": exp.problem_name,
        "systemPrompt": exp.system_prompt,
        "problemPrompt": exp.problem_prompt,
        "propertyType": exp.property_type,
        "customPropertyName": exp.custom_property_name,
        "customPropertyDesc": exp.custom_property_desc,
        "customPropertyAscending": exp.custom_property_ascending,
        "customization": customization,
        "treeNodes": exp.tree_nodes,
        "edges": exp.edges,
        "metricsHistory": exp.metrics_history,
        "visibleMetrics": exp.visible_metrics,
        "graphState": exp.graph_state,
        "autoZoom": exp.auto_zoom,
        "sidebarState": sidebar,
        "experimentContext": exp.experiment_context,
    }


def project_to_dict(proj: models.Project) -> Dict[str, Any]:
    """Convert Project model to dict matching frontend format"""
    return {
        "id": proj.id,
        "name": proj.name,
        "createdAt": proj.created_at.isoformat(),
        "lastModified": proj.last_modified.isoformat(),
        "experiments": [experiment_to_dict(exp) for exp in proj.experiments],
    }


@router.get("/", response_model=List[ProjectData])
async def get_projects(
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user),
):
    """Get all projects for the current user."""
    try:
        logger.info(f"Loading projects for user: {user}")
        
        if db is None:
            logger.warning("Database not available")
            return []
        
        query = (
            select(models.Project)
            .options(selectinload(models.Project.experiments))
            .where(models.Project.user == user)
            .order_by(models.Project.last_modified.desc())
        )

        result = await db.execute(query)
        projects = result.scalars().all()
        
        # Convert to dicts
        project_dicts = []
        for proj in projects:
            project_dicts.append(project_to_dict(proj))
        
        logger.info(f"Returning {len(project_dicts)} projects with total {sum(len(p['experiments']) for p in project_dicts)} experiments")
        return project_dicts
    except Exception as e:
        logger.error(f"Error in get_projects: {e}")
        import traceback
        logger.error(traceback.format_exc())
        raise


@router.post("/", response_model=ProjectData)
async def create_project(
    request: CreateProjectRequest,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Create a new project"""
    project_id = generate_id("proj")
    now = datetime.utcnow()
    
    project = models.Project(
        id=project_id,
        user=user,
        name=request.name,
        created_at=now,
        last_modified=now,
    )
    
    db.add(project)
    await db.commit()
    await db.refresh(project)
    
    return {
        "id": project.id,
        "name": project.name,
        "createdAt": project.created_at.isoformat(),
        "lastModified": project.last_modified.isoformat(),
        "experiments": [],
    }


@router.delete("/all")
async def delete_all_projects(
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user),
):
    """Delete all projects and experiments for the current user.

    Uses bulk SQL DELETEs to avoid ORM cascade/autoflush conflicts
    with concurrent readers (MariaDB error 1020).
    """
    # Collect project IDs for this user
    result = await db.execute(
        select(models.Project.id).where(models.Project.user == user)
    )
    project_ids = [row[0] for row in result.all()]
    count = len(project_ids)

    if project_ids:
        # Bulk-delete experiments first, then projects (no ORM cascade needed)
        await db.execute(
            delete(models.Experiment).where(
                models.Experiment.project_id.in_(project_ids)
            )
        )
        await db.execute(
            delete(models.Project).where(
                models.Project.id.in_(project_ids)
            )
        )
    await db.commit()
    logger.info(f"Deleted {count} projects for user {user}")
    return {"status": "deleted", "count": count}


@router.delete("/{project_id}")
async def delete_project(
    project_id: str,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Delete a project and all its experiments"""
    result = await db.execute(
        select(models.Project).where(
            models.Project.id == project_id,
            models.Project.user == user
        )
    )
    project = result.scalar_one_or_none()
    
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    
    await db.delete(project)
    await db.commit()
    
    return {"status": "deleted"}


@router.put("/{project_id}")
async def update_project(
    project_id: str,
    request: CreateProjectRequest,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Update project name"""
    result = await db.execute(
        select(models.Project).where(
            models.Project.id == project_id,
            models.Project.user == user
        )
    )
    project = result.scalar_one_or_none()
    
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    
    project.name = request.name
    project.last_modified = datetime.utcnow()
    
    await db.commit()
    await db.refresh(project)
    
    # Load experiments
    exp_result = await db.execute(
        select(models.Experiment)
        .where(models.Experiment.project_id == project.id)
        .order_by(models.Experiment.created_at)
    )
    project.experiments = exp_result.scalars().all()
    
    return project_to_dict(project)


@router.post("/{project_id}/experiments", response_model=ExperimentData)
async def create_experiment(
    project_id: str,
    request: CreateExperimentRequest,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Create a new experiment in a project"""
    logger.info(f"Creating experiment '{request.name}' in project {project_id} for user {user}")
    
    # Verify project exists
    result = await db.execute(
        select(models.Project).where(
            models.Project.id == project_id,
            models.Project.user == user
        )
    )
    project = result.scalar_one_or_none()
    
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    
    experiment_id = generate_id("exp")
    now = datetime.utcnow()
    
    experiment = models.Experiment(
        id=experiment_id,
        project_id=project_id,
        user=user,
        name=request.name,
        created_at=now,
        last_modified=now,
    )
    
    db.add(experiment)
    
    # Note: project.last_modified is intentionally NOT updated here.
    # Updating the shared project row concurrently from multiple experiment
    # creates/updates causes MariaDB error 1020.  The 2-second sidebar poll
    # keeps the project list current without relying on this timestamp.
    
    for attempt in range(_MAX_RETRIES):
        try:
            await db.commit()
            await db.refresh(experiment)
            logger.info(f"Created experiment {experiment_id}")
            return experiment_to_dict(experiment)
        except OperationalError as e:
            if "1020" in str(e) and attempt < _MAX_RETRIES - 1:
                logger.warning(f"Concurrent modification creating experiment {experiment_id}, retry {attempt + 1}/{_MAX_RETRIES}")
                await db.rollback()
                await asyncio.sleep(0.05 * (attempt + 1))
                # Re-add the experiment after rollback clears it from the session
                experiment = models.Experiment(
                    id=experiment_id,
                    project_id=project_id,
                    user=user,
                    name=request.name,
                    created_at=now,
                    last_modified=now,
                )
                db.add(experiment)
                continue
            raise


@router.put("/{project_id}/experiments/{experiment_id}")
async def update_experiment(
    project_id: str,
    experiment_id: str,
    request: UpdateExperimentRequest,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Update an experiment's full state"""
    exp_data = request.experiment

    for attempt in range(_MAX_RETRIES):
        try:
            result = await db.execute(
                select(models.Experiment).where(
                    models.Experiment.id == experiment_id,
                    models.Experiment.project_id == project_id,
                    models.Experiment.user == user
                )
            )
            experiment = result.scalar_one_or_none()
            
            if not experiment:
                raise HTTPException(status_code=404, detail="Experiment not found")
            
            # Update all fields from request.  For tree_nodes,
            # edges, and smiles: protect non-empty DB data from being
            # overwritten by an empty payload (e.g. a stale frontend
            # save after session restore when the backend already
            # persisted computation results via save_session_to_db).
            experiment.name = exp_data.name
            experiment.last_modified = datetime.utcnow()
            experiment.is_running = exp_data.isRunning
            experiment.problem_type = exp_data.problemType
            experiment.problem_name = exp_data.problemName
            experiment.system_prompt = exp_data.systemPrompt
            experiment.problem_prompt = exp_data.problemPrompt
            experiment.property_type = exp_data.propertyType
            experiment.custom_property_name = exp_data.customPropertyName
            experiment.custom_property_desc = exp_data.customPropertyDesc
            experiment.custom_property_ascending = exp_data.customPropertyAscending
            experiment.metrics_history = exp_data.metricsHistory
            experiment.visible_metrics = exp_data.visibleMetrics
            experiment.auto_zoom = exp_data.autoZoom
            experiment.experiment_context = exp_data.experimentContext

            # Guard: don't overwrite non-empty sidebar_state / graph_state
            # with empty data — same rationale as tree_nodes/edges below.
            incoming_sidebar_msgs = []
            if exp_data.sidebarState and isinstance(exp_data.sidebarState, dict):
                incoming_sidebar_msgs = exp_data.sidebarState.get("messages", [])
            existing_sidebar_msgs = []
            if experiment.sidebar_state and isinstance(experiment.sidebar_state, dict):
                existing_sidebar_msgs = experiment.sidebar_state.get("messages", [])

            if incoming_sidebar_msgs or not existing_sidebar_msgs:
                experiment.sidebar_state = exp_data.sidebarState
            else:
                logger.info(
                    f"update_experiment: keeping existing {len(existing_sidebar_msgs)} "
                    f"sidebar messages for {experiment_id} (incoming was empty)"
                )

            # Guard graph_state: preserve sidebarMessages if incoming
            # graph_state has none but existing one does.
            incoming_gs = exp_data.graphState or {}
            existing_gs = experiment.graph_state or {}
            incoming_gs_msgs = incoming_gs.get("sidebarMessages") or []
            existing_gs_msgs = existing_gs.get("sidebarMessages") or []
            if not incoming_gs_msgs and existing_gs_msgs:
                merged_gs = {**incoming_gs, "sidebarMessages": existing_gs_msgs}
                experiment.graph_state = merged_gs
                logger.info(
                    f"update_experiment: keeping existing {len(existing_gs_msgs)} "
                    f"graph_state.sidebarMessages for {experiment_id}"
                )
            else:
                experiment.graph_state = exp_data.graphState

            # Guard: don't overwrite non-empty data with empty data.
            # This prevents a stale frontend save (e.g. from session
            # restore after a browser close) from erasing results that
            # the backend saved via save_session_to_db.
            incoming_nodes = exp_data.treeNodes or []
            incoming_edges = exp_data.edges or []
            incoming_smiles = exp_data.smiles or ""

            if incoming_nodes or not experiment.tree_nodes:
                experiment.tree_nodes = exp_data.treeNodes
            else:
                logger.info(
                    f"update_experiment: keeping existing {len(experiment.tree_nodes)} "
                    f"tree_nodes for {experiment_id} (incoming was empty)"
                )
            if incoming_edges or not experiment.edges:
                experiment.edges = exp_data.edges
            else:
                logger.info(
                    f"update_experiment: keeping existing {len(experiment.edges)} "
                    f"edges for {experiment_id} (incoming was empty)"
                )
            if incoming_smiles or not experiment.smiles:
                experiment.smiles = exp_data.smiles
            else:
                logger.info(
                    f"update_experiment: keeping existing smiles for "
                    f"{experiment_id} (incoming was empty)"
                )
            
            # Note: project.last_modified is intentionally NOT updated here
            # to avoid contention on the shared project row when multiple
            # experiments in the same project are updated concurrently.
            
            await db.commit()
            await db.refresh(experiment)
            
            return experiment_to_dict(experiment)
        except OperationalError as e:
            if "1020" in str(e) and attempt < _MAX_RETRIES - 1:
                logger.warning(f"Concurrent modification on experiment {experiment_id}, retry {attempt + 1}/{_MAX_RETRIES}")
                await db.rollback()
                await asyncio.sleep(0.05 * (attempt + 1))  # brief backoff
                continue
            raise


@router.delete("/{project_id}/experiments/{experiment_id}")
async def delete_experiment(
    project_id: str,
    experiment_id: str,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Delete an experiment"""
    result = await db.execute(
        select(models.Experiment).where(
            models.Experiment.id == experiment_id,
            models.Experiment.project_id == project_id,
            models.Experiment.user == user
        )
    )
    experiment = result.scalar_one_or_none()
    
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")
    
    await db.delete(experiment)
    
    # Update project's last_modified
    proj_result = await db.execute(
        select(models.Project).where(models.Project.id == project_id)
    )
    project = proj_result.scalar_one_or_none()
    if project:
        project.last_modified = datetime.utcnow()
    
    await db.commit()
    
    return {"status": "deleted"}


@router.put("/{project_id}/experiments/{experiment_id}/running")
async def set_experiment_running(
    project_id: str,
    experiment_id: str,
    is_running: bool,
    db: AsyncSession = Depends(get_db),
    user: str = Depends(get_forwarded_user)
):
    """Update experiment's running status"""
    for attempt in range(_MAX_RETRIES):
        try:
            result = await db.execute(
                select(models.Experiment).where(
                    models.Experiment.id == experiment_id,
                    models.Experiment.project_id == project_id,
                    models.Experiment.user == user
                )
            )
            experiment = result.scalar_one_or_none()
            
            if not experiment:
                raise HTTPException(status_code=404, detail="Experiment not found")
            
            experiment.is_running = is_running
            experiment.last_modified = datetime.utcnow()
            
            await db.commit()
            
            return {"status": "updated"}
        except OperationalError as e:
            if "1020" in str(e) and attempt < _MAX_RETRIES - 1:
                logger.warning(f"Concurrent modification on experiment {experiment_id} (running), retry {attempt + 1}/{_MAX_RETRIES}")
                await db.rollback()
                await asyncio.sleep(0.05 * (attempt + 1))
                continue
            raise
