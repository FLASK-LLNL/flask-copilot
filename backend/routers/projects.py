################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################
"""
Serves routes for database interaction for projects and experiments (``/api/projects/*``).
"""
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select
from sqlalchemy.orm import selectinload
from typing import List
import time
import random
from datetime import datetime

from backend.database.engine import get_db
from backend.auth import get_forwarded_user
from backend.database import models, schemas

router = APIRouter(prefix="/api/projects", tags=["projects"])


def generate_id(prefix: str) -> str:
    """Generate unique ID"""
    timestamp = int(time.time() * 1000)
    random_str = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=9))
    return f"{prefix}_{timestamp}_{random_str}"


@router.get("/", response_model=List[schemas.Project])
async def get_projects(user: str = Depends(get_forwarded_user), db: AsyncSession = Depends(get_db)):
    """Load all projects for the authenticated user"""
    result = await db.execute(
        select(models.Project)
        .where(models.Project.user == user)
        .options(selectinload(models.Project.experiments))
        .order_by(models.Project.last_modified.desc())
    )
    projects = result.scalars().all()
    return projects


@router.post("/", response_model=schemas.Project)
async def create_project(
    project_data: schemas.ProjectCreate, user: str = Depends(get_forwarded_user), db: AsyncSession = Depends(get_db)
):
    """Create a new project for the authenticated user"""
    new_project = models.Project(
        id=generate_id("project"),
        user=user,
        name=project_data.name,
        created_at=datetime.utcnow(),
        last_modified=datetime.utcnow(),
    )
    db.add(new_project)
    await db.commit()
    await db.refresh(new_project)
    return new_project


@router.put("/{project_id}", response_model=schemas.Project)
async def update_project(
    project_id: str,
    project_data: schemas.ProjectUpdate,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db),
):
    """Update a project (only if owned by user)"""
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()

    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    project.name = project_data.name
    project.last_modified = datetime.utcnow()

    await db.commit()
    await db.refresh(project)
    return project


@router.delete("/{project_id}")
async def delete_project(project_id: str, user: str = Depends(get_forwarded_user), db: AsyncSession = Depends(get_db)):
    """Delete a project (only if owned by user)"""
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()

    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    await db.delete(project)
    await db.commit()
    return {"success": True}


@router.post("/{project_id}/experiments", response_model=schemas.Experiment)
async def create_experiment(
    project_id: str,
    experiment_data: schemas.ExperimentCreate,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db),
):
    """Create a new experiment (only if user owns project)"""
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()

    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    new_experiment = models.Experiment(
        id=generate_id("exp"),
        project_id=project_id,
        user=user,
        name=experiment_data.name,
        created_at=datetime.utcnow(),
        last_modified=datetime.utcnow(),
    )

    project.last_modified = datetime.utcnow()

    db.add(new_experiment)
    await db.commit()
    await db.refresh(new_experiment)
    return new_experiment


@router.put("/{project_id}/experiments/{experiment_id}", response_model=schemas.Experiment)
async def update_experiment(
    project_id: str,
    experiment_id: str,
    experiment_data: schemas.ExperimentUpdate,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db),
):
    """Update an experiment (only if user owns it)"""
    result = await db.execute(
        select(models.Experiment).where(
            models.Experiment.id == experiment_id,
            models.Experiment.project_id == project_id,
            models.Experiment.user == user,
        )
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    # Update all provided fields
    update_data = experiment_data.model_dump(exclude_unset=True, by_alias=False)
    for field, value in update_data.items():
        if hasattr(experiment, field):
            setattr(experiment, field, value)

    experiment.last_modified = datetime.utcnow()

    # Update project's last_modified
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()
    if project:
        project.last_modified = datetime.utcnow()

    await db.commit()
    await db.refresh(experiment)
    return experiment


@router.delete("/{project_id}/experiments/{experiment_id}")
async def delete_experiment(
    project_id: str, experiment_id: str, user: str = Depends(get_forwarded_user), db: AsyncSession = Depends(get_db)
):
    """Delete an experiment (only if user owns it)"""
    result = await db.execute(
        select(models.Experiment).where(
            models.Experiment.id == experiment_id,
            models.Experiment.project_id == project_id,
            models.Experiment.user == user,
        )
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    await db.delete(experiment)

    # Update project's last_modified
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()
    if project:
        project.last_modified = datetime.utcnow()

    await db.commit()
    return {"success": True}


@router.patch("/{project_id}/experiments/{experiment_id}/running")
async def set_experiment_running(
    project_id: str,
    experiment_id: str,
    is_running: bool,
    user: str = Depends(get_forwarded_user),
    db: AsyncSession = Depends(get_db),
):
    """Set experiment running status (only if user owns it)"""
    result = await db.execute(
        select(models.Experiment).where(
            models.Experiment.id == experiment_id,
            models.Experiment.project_id == project_id,
            models.Experiment.user == user,
        )
    )
    experiment = result.scalar_one_or_none()

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    experiment.is_running = is_running
    experiment.last_modified = datetime.utcnow()

    # Update project's last_modified
    result = await db.execute(
        select(models.Project).where(models.Project.id == project_id, models.Project.user == user)
    )
    project = result.scalar_one_or_none()
    if project:
        project.last_modified = datetime.utcnow()

    await db.commit()
    return {"success": True, "is_running": is_running}
