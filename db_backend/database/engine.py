################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from __future__ import annotations

import os
import ssl
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from dotenv import load_dotenv
from sqlalchemy.exc import OperationalError
from sqlalchemy.ext.asyncio import (
    AsyncEngine,
    AsyncSession,
    async_sessionmaker,
    create_async_engine,
)
from sqlalchemy.orm import declarative_base

# Load credentials from .env file (user-local, gitignored)
# Override with: FLASK_ENV_FILE=.env-local python3 ...
_project_root = Path(__file__).resolve().parent.parent.parent
_env_file = _project_root / os.getenv("FLASK_ENV_FILE", ".env")
load_dotenv(dotenv_path=_env_file)

Base = declarative_base()


@dataclass(frozen=True)
class DatabaseRuntime:
    """Process-wide database runtime resources."""

    engine: AsyncEngine | None
    session_factory: async_sessionmaker[AsyncSession] | None


def _resolve_database_url() -> tuple[bool, str]:
    """Resolve the async database URL from environment settings."""
    use_sqlite_fallback = os.getenv("USE_SQLITE_FALLBACK", "0") == "1"
    if use_sqlite_fallback:
        sqlite_path = os.getenv("SQLITE_PATH", "sqlite:///flaskcopilot.db")
        return (
            True,
            sqlite_path.replace("sqlite:///", "sqlite+aiosqlite:///"),
        )

    # All credentials loaded from .env file (see .env.example for template)
    db_user = os.getenv("DB_USER")
    db_password = os.getenv("DB_PASSWORD", "")  # Empty password allowed for local dev
    db_host = os.getenv("DB_HOST")
    db_port = os.getenv("DB_PORT")
    db_name = os.getenv("DB_NAME")

    missing = [
        key
        for key, value in {
            "DB_USER": db_user,
            "DB_HOST": db_host,
            "DB_PORT": db_port,
            "DB_NAME": db_name,
        }.items()
        if not value
    ]
    if missing:
        print(f"WARNING: Missing database config in .env: {', '.join(missing)}")
        print("See .env.example for the required variables.")
        print("Falling back to SQLite...")
        sqlite_path = os.getenv("SQLITE_PATH", "sqlite:///flaskcopilot.db")
        return (
            True,
            sqlite_path.replace("sqlite:///", "sqlite+aiosqlite:///"),
        )

    return (
        False,
        f"mysql+aiomysql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}",
    )


def _build_connect_args(use_sqlite_fallback: bool) -> dict[str, Any]:
    """Build connection arguments for the selected backend."""
    if use_sqlite_fallback:
        return {"check_same_thread": False}

    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE
    return {"ssl": ssl_context}


def _build_engine_kwargs(
    use_sqlite_fallback: bool, connect_args: dict[str, Any]
) -> dict[str, Any]:
    """Common engine settings, without applying MariaDB-only pool options to SQLite."""
    kwargs: dict[str, Any] = {
        "echo": False,
        "connect_args": connect_args,
    }
    if not use_sqlite_fallback:
        kwargs.update(
            {
                "pool_size": 10,
                "max_overflow": 20,
                "pool_pre_ping": True,
                "pool_recycle": 3600,
            }
        )
    return kwargs


@lru_cache(maxsize=1)
def _get_database_runtime() -> DatabaseRuntime:
    """Create and cache database runtime resources lazily."""
    try:
        use_sqlite_fallback, async_database_url = _resolve_database_url()
        connect_args = _build_connect_args(use_sqlite_fallback)
        engine_kwargs = _build_engine_kwargs(use_sqlite_fallback, connect_args)

        async_engine = create_async_engine(async_database_url, **engine_kwargs)
        session_factory = async_sessionmaker(
            async_engine, class_=AsyncSession, expire_on_commit=False
        )
        return DatabaseRuntime(
            engine=async_engine,
            session_factory=session_factory,
        )
    except (OperationalError, Exception) as exc:
        print(f"Warning: Could not connect to database: {exc}")
        print("Database features will be disabled.")
        return DatabaseRuntime(engine=None, session_factory=None)


def get_async_engine() -> AsyncEngine | None:
    """Create and cache the async engine lazily.

    The async engine is the primary runtime path for API and WebSocket-backed
    persistence.
    """
    return _get_database_runtime().engine


def get_async_session_factory() -> async_sessionmaker[AsyncSession] | None:
    """Return the lazily initialized async session factory."""
    return _get_database_runtime().session_factory


async def get_db():
    """Dependency for getting async database sessions."""
    session_factory = get_async_session_factory()
    if session_factory is None:
        yield None
        return

    async with session_factory() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
        finally:
            await session.close()
