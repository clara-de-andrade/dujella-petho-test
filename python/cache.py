#!/usr/bin/python3

import os
import pickle


def update_cache(new_cache, filename="cache.pkl"):
    if not os.path.exists(filename):
        with open(filename, "x") as f:
            pass

    cache = {}

    with open(filename, "rb") as f:
        try:
            cache = pickle.load(f)
        except EOFError:
            pass
    
    cache.update(new_cache)

    with open(filename, "wb") as f:
        pickle.dump(cache, f)

    return cache


def load_cache(key=None, default_value=None, filename="cache.pkl"): 
    if not os.path.exists(filename):
        with open(filename, "x") as f:
            pass

    cache = {}

    with open(filename, "rb") as f:
        try:
            cache = pickle.load(f)
        except EOFError:
            pass

    if key:
        return cache.get(key, default_value)
    else:
        return cache


def purge_cache(filename="cache.pkl"):
    if not os.path.exists(filename):
        return

    with open(filename, "wb") as f:
        pickle.dump({}, f)

    return {}