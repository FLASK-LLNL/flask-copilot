# MariaDB Setup for Flask Copilot

## Installation Summary

✅ MariaDB 12.1.2 has been installed and configured on your system.

## Service Status

MariaDB is running as a background service and will automatically start on system boot.

### Service Commands

```bash
# Check status
brew services info mariadb

# Start MariaDB
brew services start mariadb

# Stop MariaDB
brew services stop mariadb

# Restart MariaDB
brew services restart mariadb
```

## Database Configuration

### Connection Details

- **Host**: localhost
- **Port**: 3306 (default)
- **Database**: `flask_experiments`
- **User**: `flask_user`
- **Password**: `flask_password`

### Python Connection

```python
# Using PyMySQL + SQLAlchemy
from db_config import get_db_engine, get_db_session

# Get engine
engine = get_db_engine()

# Get session
session = get_db_session()
```

Connection URL format:
```
mysql+pymysql://flask_user:flask_password@localhost:3306/flask_experiments
```

## MariaDB CLI Access

```bash
# Connect to MariaDB
mariadb

# Connect to specific database
mariadb flask_experiments

# Execute query directly
mariadb -e "SHOW DATABASES;"

# Connect with username/password
mariadb -u flask_user -pflask_password flask_experiments
```

## Example Usage

See [db_example.py](./db_example.py) for a complete example including:
- Database initialization
- Creating tables with SQLAlchemy ORM
- Saving experiment data
- Querying experiments
- Updating experiment status

Run the example:
```bash
python db_example.py
```

## Database Schema

The example includes an `experiments` table with:
- `id` - Primary key
- `name` - Experiment name
- `description` - Experiment description
- `parameters` - JSON field for experiment parameters
- `results` - JSON field for experiment results
- `status` - Current status (pending/running/completed/failed)
- `created_at` - Creation timestamp
- `updated_at` - Last update timestamp

## Installed Python Packages

- `pymysql` - MySQL database adapter
- `sqlalchemy` - SQL toolkit and ORM

## Security Notes

⚠️ **Important**: The default password is `flask_password`. For production use:

1. Change the password:
```bash
mariadb -e "ALTER USER 'flask_user'@'localhost' IDENTIFIED BY 'your_secure_password';"
```

2. Use environment variables:
```python
import os
DATABASE_URL = os.getenv('DATABASE_URL')
```

3. Store credentials in a `.env` file (add to `.gitignore`):
```
DATABASE_URL=mysql+pymysql://flask_user:your_secure_password@localhost:3306/flask_experiments
```

## Backup and Restore

### Backup
```bash
mariadb-dump -u flask_user -p flask_experiments > backup.sql
```

### Restore
```bash
mariadb -u flask_user -p flask_experiments < backup.sql
```

## Useful Queries

```sql
-- Show all databases
SHOW DATABASES;

-- Use database
USE flask_experiments;

-- Show all tables
SHOW TABLES;

-- Describe table structure
DESCRIBE experiments;

-- View all experiments
SELECT * FROM experiments;

-- Count experiments by status
SELECT status, COUNT(*) FROM experiments GROUP BY status;
```

## Configuration Files

- Database config: [db_config.py](./db_config.py)
- Example usage: [db_example.py](./db_example.py)
- MariaDB config: `/Users/marathe1/local/etc/my.cnf` (if needed)

## Troubleshooting

### Connection refused
```bash
# Check if MariaDB is running
brew services list | grep mariadb

# Check port is listening
lsof -i :3306
```

### Access denied
```bash
# Verify user exists
mariadb -e "SELECT User, Host FROM mysql.user WHERE User='flask_user';"

# Recreate user if needed
mariadb -e "DROP USER IF EXISTS 'flask_user'@'localhost';"
mariadb -e "CREATE USER 'flask_user'@'localhost' IDENTIFIED BY 'flask_password';"
mariadb -e "GRANT ALL PRIVILEGES ON flask_experiments.* TO 'flask_user'@'localhost';"
mariadb -e "FLUSH PRIVILEGES;"
```

## Data Directory

MariaDB data is stored at: `/Users/marathe1/local/var/mysql`

## Next Steps

1. Integrate database saving into your ChARGe backend
2. Add experiment tracking to the web UI
3. Implement experiment history viewing
4. Add export/import functionality
5. Set up automated backups

## Documentation

- MariaDB Documentation: https://mariadb.com/kb/en/
- SQLAlchemy Documentation: https://docs.sqlalchemy.org/
- PyMySQL Documentation: https://pymysql.readthedocs.io/
