from datetime import datetime, timezone
from sqlalchemy import func, Column, DateTime, ForeignKey, JSON
from sqlalchemy.orm import (
    mapped_column,
    relationship,
    DeclarativeBase,
    Mapped,
    MappedAsDataclass,
    Session,
)
from pydantic import BaseModel, ConfigDict, Field
from typing import Any, List, Optional
import uuid


# The silly "one-true-base" base class to make all the SQLA stuff
# work correctly. Everyone calls it "Base", so we do too.
class Base(DeclarativeBase):
    pass


class TimestampSQLMixin:
    """Add columns tracking creation and modification times for objects."""

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        server_default=func.now(),
        nullable=False,
    )
    last_modified: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        onupdate=lambda: datetime.now(timezone.utc),
        server_default=func.now(),
        server_onupdate=func.now(),
        nullable=False,
    )


# We need User notions now. Yippee. Starting SUPER simple. I don't
# care about name/email or anything like that. Just a UUID id and a
# list of owned projects.
class User(Base):
    __tablename__ = "users"
    id: Mapped[uuid.UUID] = mapped_column(primary_key=True, default=uuid.uuid4)
    name: Mapped[str]

    # Make life easier. Deleting a user deletes their projects.
    projects: Mapped[List["Project"]] = relationship(
        cascade="all, delete-orphan",
        order_by="Project.created_at",
        lazy="selectin",
    )


# "Project" is the top-level operational concept. It holds experiments
# grouped under a user-defined "name".
#
# FIXME: Add event listeners so that when experiments get
# inserted/deleted/updated, Project 'last_modified' get updated
# appropriately.
class Project(Base, TimestampSQLMixin):
    __tablename__ = "projects"

    id: Mapped[uuid.UUID] = mapped_column(primary_key=True, default=uuid.uuid4)
    user_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("users.id"))

    name: Mapped[str]

    experiments: Mapped[List["Experiment"]] = relationship(
        back_populates="project",
        cascade="all, delete-orphan",
        order_by="Experiment.created_at",
        lazy="selectin",
    )


# An "Experiment" is the primary repository of contextual information.


class Experiment(Base, TimestampSQLMixin):
    __tablename__ = "experiments"

    id: Mapped[uuid.UUID] = mapped_column(primary_key=True, default=uuid.uuid4)
    project_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("projects.id"))

    data: Mapped[dict] = mapped_column(JSON)

    project: Mapped["Project"] = relationship(back_populates="experiments")


##############################
# PYDANTIC MODELS (schemas)
##############################


# Timestamps are database artifacts, not user-settables. So we mark
# these 'required' and will only use them for the "*Response" models
# below.
class TimestampMixin:
    created_at: datetime
    last_modified: datetime


class UserBase(BaseModel):
    name: str


class UserCreate(UserBase):
    pass


class UserResponse(UserBase):
    id: uuid.UUID


# To support the different "modes" with which the front-end does/will
# deal with this data, we present a few forms of these models.
# Projects can include no experiment data, experiment "metadata", or
# full experiment data.


class ProjectBase(BaseModel):
    name: str


class ProjectCreate(ProjectBase):
    pass


# This is basically "ProjectMetadata". If we just need a simple
# reference to the Project ID or Name, this is sufficent.
class ProjectResponse(ProjectBase, TimestampMixin):
    id: uuid.UUID
    user_id: uuid.UUID
    model_config = ConfigDict(from_attributes=True)


# This provides a full set of fully populated experiments, with their
# full contexts fully joined at the database level. This most easily
# shims into the copilot front-end, especially when cached in
# LocalStorage.
class ProjectResponseWithExperiments(ProjectResponse):
    experiments: List["ExperimentResponse"] = []


# The only thing we can change on a project is the name.
class ProjectUpdate(BaseModel):
    name: str | None = Field(default=None)


# NOTE (trb): Currently, the copilot frontend only passes a "name"
# field to the createProject function. So that's all we require here.
# The full context gets added with an "update", so we defer the
# addition of that field until that point.
class ExperimentCreate(BaseModel):
    name: str


# FIXME (trb): Flesh out the context.
#
# For now, the "name" that gets passed in the ExperimentCreate
# internally gets wrapped into the "data" dict, hence it not being
# explicit here.
class ExperimentResponse(BaseModel, TimestampMixin):
    id: uuid.UUID
    project_id: uuid.UUID
    data: dict[str, Any]  # FIXME: Is there a better type for JSON data?

    model_config = ConfigDict(from_attributes=True)


# FIXME (trb): Flesh out the context.
class ExperimentUpdate(BaseModel):
    data: dict[str, Any] | None = Field(default=None)
