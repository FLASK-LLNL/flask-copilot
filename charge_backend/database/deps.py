import uuid
from fastapi import Depends, HTTPException
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from typing import Annotated

from .database import get_session
from .models import Experiment, Project, User

GetSession = Annotated[AsyncSession, Depends(get_session)]


def get_current_user_name():
    return "nobody"


async def get_current_user(session: GetSession) -> User:
    db_user = await session.scalar(
        select(User).where(User.name == get_current_user_name())
    )
    return db_user


CurrentUser = Annotated[User, Depends(get_current_user)]


async def validate_project_id(
    session: GetSession, user: CurrentUser, project_id: uuid.UUID
) -> Project:
    """A helper function to verify the project id is valid and
    writeable by this user.

    """
    db_proj = await session.get(Project, project_id)
    if db_proj is None:
        raise HTTPException(status_code=404, detail="Invalid project id")
    if db_proj.user_id != user.id:
        raise HTTPException(status_code=403, detail="Bad project permissions")
    return db_proj


ValidatedProject = Annotated[Project, Depends(validate_project_id)]


async def validate_experiment_id(
    session: GetSession, project: ValidatedProject, experiment_id: uuid.UUID
) -> Experiment:
    """A helper function to verify the experiment id relates to a
    valid experiment that is owned by the specified project, which is
    validated to belong to the current user.

    """
    db_experiment = await session.get(Experiment, experiment_id)
    if db_experiment is None:
        raise HTTPException(status_code=404, detail="Invalid experiment id")
    if db_experiment.project_id != project.id:
        raise HTTPException(status_code=409, detail="Bad experiment")
    # FIXME (trb): I *think* 409 is the right code here? The user is
    # valid, and the user has permissions for this project. The issue
    # is that there's a mismatch between the experiment id and the
    # project that owns it.
    return db_experiment


ValidatedExperiment = Annotated[Experiment, Depends(validate_experiment_id)]
