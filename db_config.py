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
DB_CONFIG = {
    'host': 'localhost',
    'port': 3306,
    'user': 'flask_user',
    'password': 'flask_password',
    'database': 'flask_experiments',
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
    return create_engine(DATABASE_URL, echo=echo, pool_pre_ping=True)

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
