"""
Database configuration for Flask Copilot experiments.

To use in your application:
    from db_config import get_db_engine, get_db_session
    
    engine = get_db_engine()
    session = get_db_session()
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import os

# Database connection configuration
# Remote LLNL LaunchIT MariaDB (via SSH tunnel)
# SSH tunnel command: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 marathe1@oslic.llnl.gov
DB_CONFIG = {
    'host': os.getenv('MARIADB_HOST', '127.0.0.1'),  # Use localhost when SSH tunnel is active
    'port': int(os.getenv('MARIADB_PORT', '32636')),
    'user': os.getenv('MARIADB_USER', 'marathe1'),
    'password': os.getenv('MARIADB_PASSWORD', 'Eked1c2OWATXtD0YhHKP5CUKh5FGlbIkTaIDGtl1vKMHSB5lrW1FmA8RJB5k0V4x0lgxNSkMAhYbxo4f'),
    'database': os.getenv('MARIADB_DATABASE', 'flaskcopilot'),
}

# Connection URL
DATABASE_URL = f"mysql+pymysql://{DB_CONFIG['user']}:{DB_CONFIG['password']}@{DB_CONFIG['host']}:{DB_CONFIG['port']}/{DB_CONFIG['database']}"

# Alternative: Use environment variables for production
# DATABASE_URL = os.getenv('DATABASE_URL', DATABASE_URL)

def get_db_engine(echo=False):
    """
    Create and return a SQLAlchemy engine.
    
    Args:
        echo (bool): If True, log all SQL statements
    
    Returns:
        sqlalchemy.engine.Engine
    """
    return create_engine(
        DATABASE_URL, 
        echo=echo, 
        pool_pre_ping=True,
        connect_args={
            'ssl': {'fake_flag_to_enable_tls': True}  # Enable SSL/TLS for remote connection
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
