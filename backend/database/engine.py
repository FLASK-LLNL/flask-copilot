################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession, async_sessionmaker
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import declarative_base
import os

# Check for MariaDB configuration, fallback to local settings
DB_USER = os.getenv("MARIADB_USER", "flask_user")
DB_PASSWORD = os.getenv("MARIADB_PASSWORD", "flask_password")
DB_HOST = os.getenv("MARIADB_HOST", "localhost")
DB_PORT = os.getenv("MARIADB_PORT", "3306")
DB_NAME = os.getenv("MARIADB_DATABASE", "flask_experiments")

DATABASE_URL = f"mysql+aiomysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

try:
    engine = create_async_engine(
        DATABASE_URL,
        echo=False,  # Set to True for SQL debug logging
        pool_size=10,
        max_overflow=20,
        pool_pre_ping=True,
        pool_recycle=3600,
    )

    AsyncSessionLocal = async_sessionmaker(engine, class_=AsyncSession, expire_on_commit=False)
except OperationalError:
    print("Warning: Could not connect to MariaDB. Database features will be disabled.")
    engine = AsyncSessionLocal = None

Base = declarative_base()


async def get_db():
    """Dependency for getting async database sessions"""
    if AsyncSessionLocal is None:
        yield None
        return
    async with AsyncSessionLocal() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
        finally:
            await session.close()
