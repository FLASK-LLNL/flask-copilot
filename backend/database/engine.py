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
import ssl
from pathlib import Path
from dotenv import load_dotenv

# Load credentials from .env file (user-local, gitignored)
load_dotenv(dotenv_path=Path(__file__).resolve().parent.parent.parent / '.env')

# Allow a local SQLite fallback so the app can run even if MariaDB is unreachable.
USE_SQLITE_FALLBACK = os.getenv("USE_SQLITE_FALLBACK", "0") == "1"

if USE_SQLITE_FALLBACK:
    SQLITE_PATH = os.getenv("SQLITE_PATH", "sqlite:///flaskcopilot.db")
    DATABASE_URL = SQLITE_PATH.replace("sqlite:///", "sqlite+pysqlite:///")
    ASYNC_DATABASE_URL = SQLITE_PATH.replace("sqlite:///", "sqlite+aiosqlite:///")
else:
    # Remote LLNL LaunchIT MariaDB configuration (via SSH tunnel)
    # SSH tunnel command: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 <user>@oslic.llnl.gov
    # All credentials loaded from .env file (see .env.example for template)
    DB_USER = os.getenv("DB_USER")
    DB_PASSWORD = os.getenv("DB_PASSWORD")
    DB_HOST = os.getenv("DB_HOST")
    DB_PORT = os.getenv("DB_PORT")
    DB_NAME = os.getenv("DB_NAME")

    # Validate all required DB settings are present in .env
    _missing = [k for k, v in {"DB_USER": DB_USER, "DB_PASSWORD": DB_PASSWORD,
                                "DB_HOST": DB_HOST, "DB_PORT": DB_PORT,
                                "DB_NAME": DB_NAME}.items() if not v]
    if _missing:
        print(f"WARNING: Missing database config in .env: {', '.join(_missing)}")
        print("See .env.example for the required variables.")
        print("Falling back to SQLite...")
        USE_SQLITE_FALLBACK = True
        SQLITE_PATH = os.getenv("SQLITE_PATH", "sqlite:///flaskcopilot.db")
        DATABASE_URL = SQLITE_PATH.replace("sqlite:///", "sqlite+pysqlite:///")
        ASYNC_DATABASE_URL = SQLITE_PATH.replace("sqlite:///", "sqlite+aiosqlite:///")
    else:
        # MariaDB over SSH tunnel with mandatory TLS; use PyMySQL/AioMySQL with lax cert check (self-signed)
        DATABASE_URL = f"mysql+pymysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
        ASYNC_DATABASE_URL = f"mysql+aiomysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"

# Common SSL settings (server uses self-signed cert; tunnel already encrypts)
if not USE_SQLITE_FALLBACK:
    SSL_CONTEXT = ssl.create_default_context()
    SSL_CONTEXT.check_hostname = False
    SSL_CONTEXT.verify_mode = ssl.CERT_NONE
    SSL_KW = {"ssl": SSL_CONTEXT}
else:
    SSL_KW = {}

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Create sync engine (SSL required by server; skip cert verification because tunnel is trusted)
try:
    connect_args = {"check_same_thread": False} if USE_SQLITE_FALLBACK else SSL_KW
    sync_engine = create_engine(
        DATABASE_URL,
        echo=False,  # Set to True for SQL debug logging
        pool_size=10,
        max_overflow=20,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args=connect_args
    )
    SyncSessionLocal = sessionmaker(bind=sync_engine, expire_on_commit=False)
    
    # Also create async engine for compatibility (may have SSL issues)
    engine = create_async_engine(
        ASYNC_DATABASE_URL,
        echo=False,
        pool_size=10,
        max_overflow=20,
        pool_pre_ping=True,
        pool_recycle=3600,
        connect_args=SSL_KW
    )
    AsyncSessionLocal = async_sessionmaker(engine, class_=AsyncSession, expire_on_commit=False)
except OperationalError:
    print("Warning: Could not connect to MariaDB. Database features will be disabled.")
    engine = AsyncSessionLocal = None
    sync_engine = SyncSessionLocal = None

Base = declarative_base()


def get_sync_db():
    """Get a synchronous database session (more reliable with SSL)"""
    if SyncSessionLocal is None:
        return None
    session = SyncSessionLocal()
    try:
        return session
    except Exception:
        session.rollback()
        raise


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
