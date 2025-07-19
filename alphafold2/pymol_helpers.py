"""Helper functions for pymol

"""

import pymol
from pymol import cmd


class Item:

    def __init__(self, value, kind):
        self.value = value
        self.kind = kind


def load_pymol_item(ref, item):
    if item.kind == "file":
        return cmd.load(item.value, ref)
    elif item.kind == "id":
        return cmd.fetch(item.value, ref)
    else:
        raise RuntimeError(f"Unknown kind {item.kind}")