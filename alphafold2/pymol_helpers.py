"""Helper functions for pymol

"""

from pymol import cmd


class Item:

    def __init__(self, value, kind, id=None):
        self.value = value
        self.kind = kind
        self.id = id if id else value


def load_pymol_item(ref, item):
    if item.kind == "file":
        return cmd.load(item.value, ref)
    elif item.kind == "id":
        return cmd.fetch(item.value, ref)
    else:
        raise RuntimeError(f"Unknown kind {item.kind}")
