import functools
import sqlite3
import os
import pandas as pd

try:
    from rdkit import Chem
except ImportError:
    Chem = None


_ZINC_STOCK_DATABASE_PATH = os.getenv("FLASK_ZINC_STOCK_DB", "/data/zinc_stock.hdf5")
_MOLECULE_DATABASE_PATH = os.getenv("FLASK_MOLECULE_DB", "/data/db/molecules.db")


if os.path.exists(_ZINC_STOCK_DATABASE_PATH):
    STOCK_DATABASE = pd.read_hdf(_ZINC_STOCK_DATABASE_PATH, "table")
    # The below line does not work in a uvicorn setting because the Pandas
    # dataframe does not survive the fork well
    # STOCK_DATABASE.set_index("inchi_key", inplace=True)  # Optimize InChI key queries
    STOCK_DATABASE = set(STOCK_DATABASE["inchi_key"])
else:
    STOCK_DATABASE = None


def is_purchasable_zinc(smiles: str) -> bool:
    """
    Checks whether a molecule (as SMILES string) is purchasable from the ZINC
    Stock file.

    :param smiles: The SMILES string representing the molecule
    :return: True if purchasable, False otherwise
    """
    if Chem is None or STOCK_DATABASE is None:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Invalid SMILES
        return False
    inchi: str = str(Chem.MolToInchi(mol))
    inchi_key = Chem.InchiToInchiKey(inchi)

    return (
        inchi_key in STOCK_DATABASE
    )  # or `STOCK_DATABASE.index` if using dataframe directly


def moldb_connect(db_path: str) -> sqlite3.Connection:
    """
    Connect to the molecule database

    :param db_path: Path to the database file
    :return: A sqlite3 connection object
    """
    uri = f"file:{db_path}?mode=ro"
    conn = sqlite3.connect(uri, uri=True)

    # Read-only safe optimizations
    conn.execute("PRAGMA cache_size=-64000")  # 64MB cache
    conn.execute("PRAGMA mmap_size=268435456")  # 256MB memory-map
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def is_purchasable_moldb(conn: sqlite3.Connection, inchi: str) -> list[str]:
    """
    Checks whether a molecule is purchasable and returns a list of sources
    that provide it.

    :param conn: A connection to a moleculedb database
    :param inchi: The InChI of the molecule to query
    :return: A list of strings representing the purchase sources
    """
    cursor = conn.cursor()

    cursor.execute(
        """
        SELECT
            b.source
        FROM buyables AS b
        JOIN molecules AS m
            ON b.molecule_id = m.molecule_id
        WHERE m.key = ?
        """,
        (inchi,),
    )

    rows = cursor.fetchall()
    result = []
    for row in rows:
        result.append(row)

    return result


MOLDB_CONNECTION = None


@functools.cache
def is_purchasable(smiles: str) -> list[str]:
    """
    Checks whether a molecule is purchasable and returns a list of sources
    that provide it.

    :param conn: A connection to a moleculedb database
    :param inchi: The InChI of the molecule to query
    :return: A list of strings representing the purchase sources
    """
    result: list[str] = []
    if Chem is None:
        return result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:  # Invalid SMILES
        return result
    inchi: str = str(Chem.MolToInchi(mol))

    global MOLDB_CONNECTION
    if MOLDB_CONNECTION is None:
        if os.path.exists(_MOLECULE_DATABASE_PATH):
            MOLDB_CONNECTION = moldb_connect(_MOLECULE_DATABASE_PATH)
    if MOLDB_CONNECTION is not None:
        result.extend(is_purchasable_moldb(MOLDB_CONNECTION, inchi))

    if STOCK_DATABASE is not None:
        inchi_key = Chem.InchiToInchiKey(inchi)
        if inchi_key in STOCK_DATABASE:
            result.append("ZINC Stock")

    return result
