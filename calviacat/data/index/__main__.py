# Licensed with the MIT License, see LICENSE for details

import argparse

from ...catalog import Catalog
from . import create_index, drop_index

parser = argparse.ArgumentParser(
    "Drop and re-create Calviacat database spatial indices."
)
parser.add_argument("db", help="Calviacat database file")
parser.add_argument("table", help="Recreate the index for this table.")
args = parser.parse_args()

db = Catalog.connect(args.db)
drop_index(db, args.table)
create_index(db, args.table)
