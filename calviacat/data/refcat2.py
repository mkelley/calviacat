# Licensed with the MIT License, see LICENSE for details

"""
Reference for MAST-ATLAS-Refcat2 data format:

https://archive.stsci.edu/hlsps/atlas-refcat2/hlsp_atlas-refcat2_atlas_ccd_all_multi_v1_readme.txt


README file for ATLAS-REFCAT2 (Tonry et al. 2018) - The ATLAS All-Sky Stellar Reference Catalog
MAST webpage: https://archive.stsci.edu/hlsp/atlas-refcat2
Refer to this HLSP with DOI: http://dx.doi.org/10.17909/t9-2p3r-7651

## Introduction
ATLAS-REFCAT2 is an all-sky reference catalog containing nearly one billion stars down to apparent magnitude m ~19.  The catalog includes PanSTARRS DR1, ATLAS Pathfinder, ATLAS re-flattened APASS, SkyMapper DR1, APASS DR9, Tycho-2, and the Yale Bright Star Catalog.  Gaia DR2 serves as the source of the astrometric solution for ATLAS-REFCAT2, with typical systematic errors of < 5 mmag RMS, although this can be as much as 20 mmag near the Galactic plane.  The ATLAS Pathfinder telescope was used to collect g,r,i photometry for stars brighter than the 14th magnitude bright limit of PanSTARRS, and to extend the reference system below -30 declination.

## Data Products
The ATLAS-REFCAT2 catalog is available in .csv files divided by declintion on the MAST project page.  It is also available in MAST CasJobs (http://mastweb.stsci.edu/mcasjobs/) for SQL access and for cross-matching with numerous other catalogs available at MAST.

### Catalog Columns
Index | Name   |   Data Type   |   Unit   |   Description
---------------------------------------------------------
0       objid      bigint          none       Object ID
1       RA 	       float      	   degrees    Right ascension from Gaia DR2, J2000, epoch 2015.5
2       Dec 	   float      	   degrees    Declination from Gaia DR2, J2000, epoch 2015.5
3       plx 	   real       	   mas	      Parallax from Gaia DR2
4       dplx 	   real       	   mas	      Parallax uncertainty from Gaia DR2
5       pmra 	   real       	   mas/yr     Proper motion in right ascension from Gaia DR2
6       dpmra      real       	   mas/yr     Proper motion uncertainty in right ascension
7       pmdec      real       	   mas/yr     Proper motion in declination from Gaia DR2
8       dpmdec     real       	   mas/yr     Proper motion uncertainty in declination
9       Gaia 	   real       	   mag	      Gaia G magnitude
10      dGaia      real       	   mag	      Gaia G magnitude uncertainty
11      BP 	       real       	   mag	      Gaia G_bp magnitude
12      dBP 	   real       	   mag	      Gaia G_bp magnitude uncertainty
13      RP   	   real       	   mag	      Gaia G_rp magnitude
14      dRP 	   real       	   mag	      Gaia G_rp magnitude uncertainty
15      Teff 	   int        	   K	      Gaia stellar effective temperature
16      AGaia      real       	   mag	      Gaia estimate of G-band extinction for this star
17      dupvar     int        	   none	      Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE
18      Ag 	       real       	   mag	      SFD estimate of total g-band extinction
19      rp1 	   real       	   arcsec     Radius where cummulative G flux exceeds 0.1 x this star
20      r1  	   real       	   arcsec     Radius where cummulative G flux exceeds 1.0 x this star
21      r10 	   real       	   arcsec     Radius where cummulative G flux exceeds 10.0 x this star
22      g   	   real       	   mag	      PanSTARRS g magnitude
23      dg  	   real       	   mag	      PanSTARRS g magnitude uncertainty
24      gchi 	   real       	   none	      chi^2 / DOF for contributors
25      gcontrib   int        	   none	      Bitmap of conributing catalogs to g
26      r   	   real       	   mag	      PanSTARRS r magnitude
27      dr   	   real       	   mag	      PanSTARRS r magnitude uncertainty
28      rchi 	   real       	   none	      chi^2 / DOF for contributors
29      rcontrib   int        	   none	      Bitmap of conributing catalogs to r
30      i    	   real       	   mag	      PanSTARRS i magnitude
31      di  	   real       	   mag	      PanSTARRS i magnitude uncertainty
32      ichi 	   real       	   none	      chi^2 / DOF for contributors
33      icontrib   int        	   none	      Bitmap of conributing catalogs to i
34      z   	   real       	   mag	      PanSTARRS z magnitude
35      dz  	   real       	   mag	      PanSTARRS z magnitude uncertainty
36      zchi 	   real       	   none	      chi^2 / DOF for contributors
37      zcontrib   int        	   none	      Bitmap of conributing catalogs to z
38      nstat      int        	   none	      Count of griz outliers rejected
39      J   	   real       	   mag	      2MASS J magnitude
40      dJ  	   real       	   mag	      2MASS J magnitude uncertainty
41      H   	   real       	   mag	      2MASS H magnitude
42      dH  	   real       	   mag	      2MASS H magnitude uncertainty
43      K   	   real       	   mag	      2MASS K magnitude
44      dK  	   real       	   mag	      2MASS K magnitude uncertainty

"""  # noqa: E501

import argparse

from ..refcat2 import RefCat2
from .index import drop_index


def add(filename: str, database: str) -> None:
    """Add a MAST-ATLAS-Refcat2 file to a calviacat database.

    The MAST-formatted ATLAS-Refcat2 data are available at
    https://archive.stsci.edu/hlsp/atlas-refcat2 .


    Parameters
    ----------
    filename : str
        MAST-formatted input CSV data file

    database : str
        Calviacat database file name name.

    drop_index : bool
        If ``True``, then the spatial index will be dropped before adding rows.
        It will not be re-created.

    """

    rc2 = RefCat2(database, offline=True, create_spatial_index=False)
    drop_index(rc2.db, "refcat2")

    # test on 100,000 rows
    # chunk 300 = 4.4 s
    # chunk 1000 = 3.5 s
    # chunk 3000 = 3.3 s
    # chunk 10000 = 3.6 s
    # chunk 30000 = 3.8 s
    chunk = 3000
    total = 0
    with open(filename, "r") as inf:
        while True:
            rows = []
            for count, line in enumerate(inf):
                cells = line.split(",")
                rows.append(
                    [
                        cells[0],
                        cells[1],
                        cells[2],
                        cells[9],
                        cells[10],
                        cells[11],
                        cells[12],
                        cells[13],
                        cells[14],
                        cells[22],
                        cells[23],
                        cells[26],
                        cells[27],
                        cells[30],
                        cells[31],
                        cells[34],
                        cells[35],
                        cells[39],
                        cells[40],
                        cells[41],
                        cells[42],
                        cells[43],
                        cells[44],
                    ]
                )
                if count >= chunk:
                    break

            if len(rows) == 0:
                break

            rc2.db.executemany(
                """
            INSERT OR IGNORE INTO {}
            VALUES({})
            """.format(rc2.table.name, ",".join("?" * len(rc2.table.columns))),
                rows,
            )
            rc2.db.commit()

            total += len(rows)

    print(f"{total} rows added.")
    print(
        "If the database is complete, create the spatial index with"
        " calviacat.data.index --table=refcat2"
    )


if __name__ == "__main__":  # pragma: no cover
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="MAST-ATLAS-Refcat2 csv file")
    parser.add_argument("db", help="Calviacat database file")
    args = parser.parse_args()

    add(args.csv, args.db)
