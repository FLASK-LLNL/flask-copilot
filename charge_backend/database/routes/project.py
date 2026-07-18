from fastapi import APIRouter
from sqlalchemy import select
from typing import List

from charge_backend.database.deps import GetSession, CurrentUser, ValidatedProject
from charge_backend.database.models import (
    Project,
    ProjectCreate,
    ProjectResponse,
    ProjectResponseWithExperiments,
    ProjectUpdate,
)


router = APIRouter(prefix="/projects", tags=["projects"])


@router.post("", response_model=ProjectResponse)
async def create_project(
    *, session: GetSession, project: ProjectCreate, user: CurrentUser
):
    proj_dict = project.model_dump()
    proj_dict["user_id"] = user.id
    db_proj = Project(**proj_dict)
    session.add(db_proj)
    await session.commit()
    await session.refresh(db_proj)
    return db_proj


@router.get("", response_model=List[ProjectResponseWithExperiments])
async def get_projects(
    *,
    session: GetSession,
    current_user: CurrentUser,
):
    await session.refresh(current_user)
    return current_user.projects


@router.get("/debug-all", response_model=List[ProjectResponseWithExperiments])
async def get_all_projects_debug(
    *,
    session: GetSession,
):
    return (await session.scalars(select(Project))).all()


@router.get("/{project_id}", response_model=ProjectResponseWithExperiments)
async def get_project(
    *,
    session: GetSession,
    db_project: ValidatedProject,
):
    return db_project


@router.put("/{project_id}", response_model=ProjectResponse)
async def update_project(
    *,
    session: GetSession,
    db_project: ValidatedProject,
    project_update: ProjectUpdate,
):
    update_dict = project_update.model_dump(exclude_unset=True)
    for key, value in update_dict.items():
        setattr(db_project, key, value)
    await session.commit()
    await session.refresh(db_project)
    return db_project


@router.delete("/{project_id}")
async def delete_project(
    *,
    session: GetSession,
    db_project: ValidatedProject,
):
    await session.delete(db_project)
    await session.commit()
    return {"ok": True}
