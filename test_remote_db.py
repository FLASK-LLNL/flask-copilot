#!/usr/bin/env python3
"""
Test remote LLNL LaunchIT MariaDB connection
"""

import os
import sys
import asyncio
import time
import ssl
from pathlib import Path
from dotenv import load_dotenv
from sqlalchemy.exc import OperationalError
from sqlalchemy import text, create_engine
from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession, async_sessionmaker
import uuid
from datetime import datetime

# Load credentials from .env file (user-local, gitignored)
# Looks in current directory, then parent directories
load_dotenv(dotenv_path=Path(__file__).resolve().parent / '.env')

# Remote LLNL LaunchIT MariaDB connection details
# When using SSH tunnel: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 marathe1@oslic.llnl.gov
# Connect to localhost (127.0.0.1) instead of the remote hostname
#
# Credentials are loaded from .env file (see .env.example for template).
# You can also set them via environment variables:
#   export DB_PASSWORD='your_password_here'
#
# For production with certificate verification:
#   export DB_SSL_CA='/path/to/ca-cert.pem'
#   export DB_SSL_VERIFY='true'
REMOTE_DB_CONFIG = {
    'host': os.getenv('DB_HOST', '127.0.0.1'),  # Use localhost when SSH tunnel is active
    'port': int(os.getenv('DB_PORT', '32636')),
    'user': os.getenv('DB_USER', 'marathe1'),
    'password': os.getenv('DB_PASSWORD'),  # Required - no default for security
    'database': os.getenv('DB_NAME', 'flaskcopilot')
}

# Validate password exists
if not REMOTE_DB_CONFIG['password']:
    print("ERROR: DB_PASSWORD environment variable is required")
    print("Set it with: export DB_PASSWORD='your_password_here'")
    sys.exit(1)


def get_ssl_config():
    """
    Get SSL configuration for database connections.
    Supports optional certificate verification via environment variables.
    
    Environment variables:
      DB_SSL_CA: Path to CA certificate file
      DB_SSL_VERIFY: Set to 'true' to require certificate verification (default: false for dev)
    """
    ssl_ca = os.getenv('DB_SSL_CA')
    ssl_verify = os.getenv('DB_SSL_VERIFY', 'false').lower() == 'true'
    
    # For development/SSH tunnel: basic SSL without verification
    # For production: provide DB_SSL_CA and set DB_SSL_VERIFY=true
    if ssl_ca and os.path.exists(ssl_ca):
        return {
            'ca': ssl_ca,
            'check_hostname': ssl_verify,
            'verify_cert': ssl_verify
        }
    else:
        # Basic SSL encryption without certificate verification
        # Suitable for SSH tunnels or internal networks
        return {'ssl': True}


def connect_with_retry(create_fn, max_retries=3, delay=2):
    """Attempt connection with retries for transient failures"""
    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            return create_fn()
        except Exception as e:
            last_error = e
            if attempt < max_retries:
                print(f"  Attempt {attempt} failed, retrying in {delay}s...")
                time.sleep(delay)
    raise last_error

# Connection URLs
SYNC_URL = f"mysql+pymysql://{REMOTE_DB_CONFIG['user']}:{REMOTE_DB_CONFIG['password']}@{REMOTE_DB_CONFIG['host']}:{REMOTE_DB_CONFIG['port']}/{REMOTE_DB_CONFIG['database']}"
ASYNC_URL = f"mysql+aiomysql://{REMOTE_DB_CONFIG['user']}:{REMOTE_DB_CONFIG['password']}@{REMOTE_DB_CONFIG['host']}:{REMOTE_DB_CONFIG['port']}/{REMOTE_DB_CONFIG['database']}"

async def test_remote_connection():
    """Test remote database connection"""
    print("=== Remote LLNL LaunchIT MariaDB Test ===\n")
    
    print("Connection Details:")
    print(f"Host: {REMOTE_DB_CONFIG['host']}")
    print(f"Port: {REMOTE_DB_CONFIG['port']}")
    print(f"User: {REMOTE_DB_CONFIG['user']}")
    print(f"Database: {REMOTE_DB_CONFIG['database']}")
    print(f"Password: {'*' * len(REMOTE_DB_CONFIG['password'])}")
    print()
    
    sync_success = False
    async_success = False
    sync_engine = None
    async_engine = None
    
    # Test synchronous connection
    print("Testing Sync Connection...")
    try:
        def create_sync_engine():
            ssl_config = get_ssl_config()
            engine = create_engine(
                SYNC_URL, 
                pool_pre_ping=True, 
                connect_args={
                    'connect_timeout': 10,
                    **ssl_config  # Merge SSL config
                }
            )
            # Test connection immediately
            with engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            return engine
        
        sync_engine = connect_with_retry(create_sync_engine)
        with sync_engine.connect() as conn:
            result = conn.execute(text("SELECT 1 as test_connection"))
            row = result.fetchone()
            print(f"Sync connection successful: {row[0]}")
            
            # Get MariaDB version
            result = conn.execute(text("SELECT VERSION() as version"))
            row = result.fetchone()
            print(f"Remote MariaDB version: {row[0]}")
            sync_success = True
            
    except Exception as e:
        print(f"Sync connection failed: {e}")
        if sync_engine:
            sync_engine.dispose()
    
    # Test asynchronous connection
    print("\nTesting Async Connection...")
    try:
        async def create_async_engine_with_test():
            ssl_config = get_ssl_config()
            engine = create_async_engine(
                ASYNC_URL,
                echo=False,
                pool_size=5,
                max_overflow=10,
                pool_pre_ping=True,
                pool_recycle=3600,
                connect_args={
                    **ssl_config  # Merge SSL config
                }
            )
            # Test connection immediately
            async with engine.begin() as conn:
                await conn.execute(text("SELECT 1"))
            return engine
        
        # Retry logic for async (with sync wrapper)
        last_error = None
        for attempt in range(1, 4):
            try:
                async_engine = await create_async_engine_with_test()
                break
            except Exception as e:
                last_error = e
                if attempt < 3:
                    print(f"  Attempt {attempt} failed, retrying in 2s...")
                    await asyncio.sleep(2)
        else:
            raise last_error
        
        async with async_engine.begin() as conn:
            result = await conn.execute(text("SELECT 1 as test_connection"))
            row = result.fetchone()
            print(f"Async connection successful: {row[0]}")
            
            # Get MariaDB version
            result = await conn.execute(text("SELECT VERSION() as version"))
            row = result.fetchone()
            print(f"Remote MariaDB version: {row[0]}")
            async_success = True
            
    except Exception as e:
        print(f"Async connection failed: {e}")
        if async_engine:
            await async_engine.dispose()
    
    return sync_success, async_success, sync_engine, async_engine

async def test_database_setup(sync_engine=None, async_engine=None):
    """Test if database and tables exist"""
    print("\n=== Database Setup Test ===\n")
    
    if not sync_engine and not async_engine:
        print("No working connection available")
        return False
        
    try:
        # Prefer sync engine since it's more reliable with SSL
        if sync_engine:
            print("Using sync connection for setup...")
            with sync_engine.connect() as conn:
                # List all databases
                result = conn.execute(text("SHOW DATABASES"))
                databases = [row[0] for row in result]
                print(f"Available databases: {databases}")
                
                # Check for tables
                result = conn.execute(text("SHOW TABLES"))
                tables = [row[0] for row in result]
                print(f"Existing tables: {tables}")
        elif async_engine:
            print("Using async connection for setup...")
            async with async_engine.begin() as conn:
                # List all databases
                result = await conn.execute(text("SHOW DATABASES"))
                databases = [row[0] for row in result]
                print(f"Available databases: {databases}")
                
                # Check for tables
                result = await conn.execute(text("SHOW TABLES"))
                tables = [row[0] for row in result]
                print(f"Existing tables: {tables}")
                    
        return True
        
    except Exception as e:
        print(f"Database setup check failed: {e}")
        import traceback
        traceback.print_exc()
        return False

async def main():
    """Main test function"""
    # Test connections
    sync_success, async_success, sync_engine, async_engine = await test_remote_connection()
    
    if not sync_success and not async_success:
        print("\nAll remote connections failed!")
        print("\nPossible issues:")
        print("- Not connected to LLNL network/VPN")
        print("- Firewall blocking access")
        print("- Server is down")
        print("- Incorrect credentials")
        return False
    
    # Test database setup
    setup_success = await test_database_setup(sync_engine, async_engine)
    
    # Summary
    print(f"\n=== Final Results ===")
    print(f"Sync Connection: {'SUCCESS' if sync_success else 'FAILED'}")
    print(f"Async Connection: {'SUCCESS' if async_success else 'FAILED'}")
    print(f"Database Setup: {'SUCCESS' if setup_success else 'FAILED'}")
    
    if sync_success or async_success:
        print(f"\nRemote LLNL LaunchIT MariaDB connection is working!")
        print(f"Ready to use: {REMOTE_DB_CONFIG['host']}:{REMOTE_DB_CONFIG['port']}")
    
    # Clean up engines
    if sync_engine:
        sync_engine.dispose()
    if async_engine:
        await async_engine.dispose()
    
    return sync_success or async_success

# Overall timeout in seconds (prevents hanging on network issues)
OVERALL_TIMEOUT = int(os.getenv('DB_TEST_TIMEOUT', '60'))

if __name__ == "__main__":
    try:
        result = asyncio.run(asyncio.wait_for(main(), timeout=OVERALL_TIMEOUT))
    except asyncio.TimeoutError:
        print(f"\n*** TEST TIMED OUT after {OVERALL_TIMEOUT}s ***")
        print("The test took too long — likely a network/firewall issue.")
        print(f"Adjust timeout with: export DB_TEST_TIMEOUT=120")
        sys.exit(2)
    sys.exit(0 if result else 1)
