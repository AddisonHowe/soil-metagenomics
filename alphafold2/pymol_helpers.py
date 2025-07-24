"""Helper functions for pymol

"""

from pymol import cmd
from tqdm import tqdm


class Item:

    def __init__(self, value, kind, id=None):
        self.value = value
        self.kind = kind
        self.id = id if id else value


def load_pymol_item(ref, item, fetch_path=None):
    if item.kind == "file":
        return cmd.load(item.value, ref)
    elif item.kind == "id":
        return cmd.fetch(item.value, ref, path=fetch_path, finish=1)
    else:
        raise RuntimeError(f"Unknown kind {item.kind}")


def get_printv(verbosity):
    def printv(s, v=1, pbar=False):
        logger = tqdm.write if pbar else print
        if v <= verbosity:
            logger(s)
    return printv
