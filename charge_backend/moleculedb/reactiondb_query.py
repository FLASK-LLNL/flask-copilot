import sqlite3
import json
from typing import List, Dict, Any


class ReactionDatabaseReader:
    def __init__(self, db_path: str, read_only: bool = True):
        """
        Initialize database reader.

        :param db_path: Path to SQLite database
        :param read_only: Open in read-only mode
        """
        if read_only:
            uri = f"file:{db_path}?mode=ro"
            self.conn = sqlite3.connect(uri, uri=True)

            # Read-only safe optimizations
            self.conn.execute("PRAGMA cache_size=-64000")  # 64MB cache
            self.conn.execute("PRAGMA mmap_size=268435456")  # 256MB memory-map
            self.conn.execute("PRAGMA temp_store=MEMORY")
        else:
            self.conn = sqlite3.connect(db_path)

            # Full optimizations (requires write access)
            self.conn.execute("PRAGMA journal_mode=WAL")
            self.conn.execute("PRAGMA synchronous=NORMAL")
            self.conn.execute("PRAGMA cache_size=-64000")
            self.conn.execute("PRAGMA mmap_size=268435456")
            self.conn.execute("PRAGMA temp_store=MEMORY")

        self.cursor = self.conn.cursor()

    def get(self, key: str) -> List[Dict[str, Any]]:
        """
        Get all entries matching the exact key.

        :param key: The key to query
        :return: List of dictionaries (empty list if no matches)
        """
        self.cursor.execute("SELECT data FROM entries WHERE key = ?", (key,))
        return [json.loads(row[0]) for row in self.cursor.fetchall()]

    def get_many(self, keys: List[str]) -> Dict[str, List[Dict[str, Any]]]:
        """
        Get all entries for multiple keys.

        :param keys: The list of keys to query
        :return: Dictionary mapping each key to its list of entries
        """
        placeholders = ",".join("?" * len(keys))
        self.cursor.execute(
            f"SELECT key, data FROM entries WHERE key IN ({placeholders})", keys
        )

        # Group results by key
        results = {}
        for key, data in self.cursor.fetchall():
            if key not in results:
                results[key] = []
            results[key].append(json.loads(data))

        # Ensure all queried keys are in result (even if empty)
        for key in keys:
            if key not in results:
                results[key] = []

        return results

    def count_by_key(self, key: str) -> int:
        """Count number of entries for a specific key."""
        self.cursor.execute("SELECT COUNT(*) FROM entries WHERE key = ?", (key,))
        return self.cursor.fetchone()[0]

    def exists(self, key: str) -> bool:
        """Check if key has any entries."""
        self.cursor.execute("SELECT 1 FROM entries WHERE key = ? LIMIT 1", (key,))
        return self.cursor.fetchone() is not None

    def total_entries(self) -> int:
        """Get total number of entries."""
        self.cursor.execute("SELECT COUNT(*) FROM entries")
        return self.cursor.fetchone()[0]

    def unique_keys(self) -> int:
        """Get number of unique keys."""
        self.cursor.execute("SELECT COUNT(DISTINCT key) FROM entries")
        return self.cursor.fetchone()[0]

    def close(self):
        """Close database connection."""
        self.conn.close()

    # Support ``with`` syntax
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
