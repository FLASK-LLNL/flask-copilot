"""
Database configuration for Flask Copilot experiments.

To use in your application:
    from db_config import get_db_engine, get_db_session
    
    engine = get_db_engine()
    session = get_db_session()

Credentials are loaded from .env file (see .env.example for template).
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# Load credentials from .env file (user-local, gitignored)
load_dotenv(dotenv_path=Path(__file__).resolve().parent / '.env')

# Database connection configuration
# Remote LLNL LaunchIT MariaDB (via SSH tunnel)
# SSH tunnel command: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 marathe1@oslic.llnl.gov
DB_CONFIG = {
    'host': os.getenv('DB_HOST', '127.0.0.1'),  # Use localhost when SSH tunnel is active
    'port': int(os.getenv('DB_PORT', '32636')),
    'user': os.getenv('DB_USER', 'marathe1'),
    'password': os.getenv('DB_PASSWORD'),  # Required - loaded from .env
    'database': os.getenv('DB_NAME', 'flaskcopilot'),
}

if not DB_CONFIG['password']:
    print("ERROR: DB_PASSWORD is required. Set it in .env file (see .env.example)")
    sys.exit(1)

# Connection URL
DATABASE_URL = f"mysql+pymysql://{DB_CONFIG['user']}:{DB_CONFIG['password']}@{DB_CONFIG['host']}:{DB_CONFIG['port']}/{DB_CONFIG['database']}"

def get_db_engine(echo=False):
    """
    Create and return a SQLAlchemy engine.
    
    Args:
        echo (bool): If True, log all SQL statements
    
    Returns:
        sqlalchemy.engine.Engine
    """
    import ssl as _ssl
    ssl_context = _ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = _ssl.CERT_NONE
    return create_engine(
        DATABASE_URL, 
        echo=echo, 
        pool_pre_ping=True,
        connect_args={
            'ssl': ssl_context
        }
    )

def get_db_session():
    """
    Create and return a SQLAlchemy session.
    
    Returns:
        sqlalchemy.orm.Session
    """
    engine = get_db_engine()
    Session = sessionmaker(bind=engine)
    return Session()

# For convenience
engine = get_db_engine()
SessionLocal = sessionmaker(bind=engine)
