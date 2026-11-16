"""
Registry of clustering methods.

Each method module defines:
    NAME        short name used in result files
    PARAM_NAME  name of its main parameter ("k" or "resolution")
    run(data, param, seed, **options) -> {"labels": ..., "runtime_sec": ...}

To add a new method, write a module with the same three things and add it here.
"""

from . import kmeans_method, leiden_method, sc3s_method

METHODS = {m.NAME: m for m in (kmeans_method, leiden_method, sc3s_method)}
