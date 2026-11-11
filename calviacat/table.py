# Licensed with the MIT License, see LICENSE for details

"""SQLite table definitions and registry."""

from dataclasses import dataclass


@dataclass
class TableDefinition:
    name: str
    column_definitions: tuple[str, str]
    objid_column: str
    ra_column: str
    dec_column: str
    filter2col: dict[str, str]

    def __post_init__(self):
        """Register the table definition with Calviacat."""
        global registry
        registry[self.name] = self

    @property
    def columns(self):
        return tuple((c[0] for c in self.column_definitions))


registry: dict[str, TableDefinition] = {}
