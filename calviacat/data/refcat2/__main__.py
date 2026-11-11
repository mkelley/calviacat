# Licensed with the MIT License, see LICENSE for details

import argparse

from . import add

parser = argparse.ArgumentParser()
parser.add_argument("csv", help="MAST-ATLAS-Refcat2 csv file")
parser.add_argument("db", help="Calviacat database file")
args = parser.parse_args()

add(args.csv, args.db)
