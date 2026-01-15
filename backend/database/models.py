################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from sqlalchemy import Column, String, DateTime, Boolean, ForeignKey, Index, Text, JSON
from sqlalchemy.orm import relationship, Mapped, mapped_column
from datetime import datetime
from typing import Optional
from backend.database.engine import Base


class Project(Base):
    __tablename__ = "projects"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)
    user: Mapped[str] = mapped_column(String(255), nullable=False, index=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)
    last_modified: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False
    )

    experiments: Mapped[list["Experiment"]] = relationship(
        "Experiment", back_populates="project", cascade="all, delete-orphan"
    )

    __table_args__ = (Index("idx_user_last_modified", "user", "last_modified"),)


class Experiment(Base):
    __tablename__ = "experiments"

    id: Mapped[str] = mapped_column(String(255), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(255), ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
    user: Mapped[str] = mapped_column(String(255), nullable=False, index=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)
    last_modified: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False
    )
    is_running: Mapped[Optional[bool]] = mapped_column(Boolean, default=False)

    # System state fields
    smiles: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    problem_type: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    problem_name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    system_prompt: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    problem_prompt: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    property_type: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    custom_property_name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    custom_property_desc: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    custom_property_ascending: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)

    # Complex JSON fields for nested data structures
    tree_nodes: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)
    edges: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)
    metrics_history: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)
    visible_metrics: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)
    graph_state: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)
    auto_zoom: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    sidebar_state: Mapped[Optional[str]] = mapped_column(JSON, nullable=True)

    # Experiment context
    experiment_context: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    project: Mapped["Project"] = relationship("Project", back_populates="experiments")

    __table_args__ = (Index("idx_user_project", "user", "project_id"),)
