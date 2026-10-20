# Licensed with the MIT License, see LICENSE for details

"""Script to drop and re-generate the spatial index."""

import sqlite3

from ...exceptions import DatabaseError
from ...table import registry as table_registry


def create_index(db: sqlite3.Connection, table: str) -> bool:
    """Create the spatial index.

    If the index already exists, nothing is done.


    Parameters
    ----------
    db : sqlite3 Connection
        The Calviacat sqlite3 database object.

    table : str, optional
        The table to index.


    Returns
    -------
    created : bool
        ``True`` if the index was created.

    """

    # the table must be in the table registry
    try:
        table_def = table_registry[table]
    except LookupError:
        raise LookupError(
            f"{table} is not one of the registered tables: "
            + ", ".join(table_registry.keys())
        )

    # check for the index and trigger
    c = db.execute(
        f"SELECT TRUE FROM sqlite_master WHERE type='table' AND name='{table}_skytree';"
    )
    index_exists = c.fetchone() is not None

    c = db.execute(
        "SELECT TRUE FROM sqlite_master "
        f"WHERE type='trigger' AND name='{table}_skytree';"
    )
    trigger_exists = c.fetchone() is not None

    # if both exist, then there is nothing to do
    if index_exists and trigger_exists:
        return False

    # if the index exists, but the trigger is missing, then there is a serious
    # issue
    if index_exists and not trigger_exists:
        raise DatabaseError(
            f"The table {table}_skytree exists but the trigger "
            f"{table}_insert does not.  Consider dropping and "
            "re-creating the index."
        )

    # neither exist, so create the (virtual) table
    db.execute(
        f"""
            CREATE VIRTUAL TABLE
            {table}_skytree USING RTREE(
            {table_def.objid_column} INTEGER PRIMARY KEY,
            x0 FLOAT, x1 FLOAT,
            y0 FLOAT, y1 FLOAT,
            z0 FLOAT, z1 FLOAT
        )"""
    )

    # insert the data to index
    cmd = f"""
        INSERT INTO {table}_skytree
            ({table_def.objid_column}, x0, x1, y0, y1, z0, z1)
        SELECT
            {table_def.objid_column},
            COS({table_def.dec_column} * 0.017453292519943295)
            * COS({table_def.ra_column} * 0.017453292519943295),
            COS({table_def.dec_column} * 0.017453292519943295)
            * COS({table_def.ra_column} * 0.017453292519943295),
            COS({table_def.dec_column} * 0.017453292519943295)
            * SIN({table_def.ra_column} * 0.017453292519943295),
            COS({table_def.dec_column} * 0.017453292519943295)
            * SIN({table_def.ra_column} * 0.017453292519943295),
            SIN({table_def.dec_column} * 0.017453292519943295),
            SIN({table_def.dec_column} * 0.017453292519943295)
        FROM {table};
    """
    db.execute(cmd)

    # trigger insert when new data is added
    db.execute(f"""
        CREATE TRIGGER IF NOT EXISTS
        {table}_insert AFTER INSERT ON {table}
        BEGIN
            INSERT INTO {table}_skytree VALUES (
                new.{table_def.objid_column},
                COS(new.{table_def.dec_column} * 0.017453292519943295)
                * COS(new.{table_def.ra_column} * 0.017453292519943295),
                COS(new.{table_def.dec_column} * 0.017453292519943295)
                * COS(new.{table_def.ra_column} * 0.017453292519943295),
                COS(new.{table_def.dec_column} * 0.017453292519943295)
                * SIN(new.{table_def.ra_column} * 0.017453292519943295),
                COS(new.{table_def.dec_column} * 0.017453292519943295)
                * SIN(new.{table_def.ra_column} * 0.017453292519943295),
                SIN(new.{table_def.dec_column} * 0.017453292519943295),
                SIN(new.{table_def.dec_column} * 0.017453292519943295)
            );
        END
    """)

    db.commit()

    return True


def drop_index(db: sqlite3.Connection, table: str):
    """Drop the spatial index and associated insert trigger.

    No error is raised if the index or trigger does not exist.


    Parameters
    ----------
    db : sqlite3 Connection
        The Calviacat sqlite3 database object.

    table : str
        Drop the index and trigger for this table.

    """

    db.execute(f"DROP TABLE IF EXISTS {table}_skytree")
    db.execute(f"DROP TRIGGER IF EXISTS {table}_insert")
    db.commit()
