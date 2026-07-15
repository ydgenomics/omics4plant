import os
import json
import argparse
import re

import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import Dataset
import pyBigwig
import pyfaidx
import logging
import pandas as pd
from typing import Callable, Dict, List

# 这个import的极致是怎样的？
from src.util import dist_print

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

# 感觉没什么用，可能是提高chr的通用性
def _chrom_aliases(chrom):
    chrom=str(chrom)
    aliases={chrom}
    match=re.search(r"(?:Chr|chr)?0?(1[0-2]|[1-9])(?!\d)", chrom)
    if match:
        idx = int(match.group(1))
        aliases.update({
            f"Chr{idx:02d}",
            f"Chr{idx}",
            f"chr{idx:02d}",
            f"chr{idx}",
            str(idx),
            f"{idx:02d}",
        })
    return aliases

def load_chromosome_features(path, normalize=True):
    """
    load chromosome level numeric features

    supported json formats:
        ...
    """
    with open(path, 'r') as f:
        payload = json.load(f)
    
    if isinstance(payload, dict) and "features" in payload:
        raw_features = payload["features"]
    elif isinstance(payload, dict):
        raw_features = {
            key: value for key, value in payload.items()
            if key not in {"feature_names", "normalize", "metadata"}
        }
    else: 
        raw_features=payload
    if not isinstance(raw_features, dict):
        raise ValueError(f"chromosome features must be a JSON object: {path}")
    
    feature_names = payload.get("feature_names") if isinstance(payload, dict) else None
    if feature_names is None:
        dict_values = [v for v in raw_features.values() if isinstance(v, dict)]
        if dict_values:
            keys = set()
            for value in dict_values:
                keys.update(value.keys())
            feature_names = sorted(keys)

    rows = {}
    for chrom, values in raw_features.items():
        if isinstance(values, dict):
            if feature_names is None:
                raise ValueError(f"Cannot infer feature names for chromosome feature dict: {path}")
            vector = [float(value.get(name, 0.0)) for name in feature_names]
        else:
            vector = [float(v) for v in values]
            if feature_names is None:
                feature_names = [f"feature_{i}" for i in range(len(vector))]
        if len(vector) != len(feature_names):
            raise ValueError(
                f"Chromosome {chrom} feature length mismatch: "
                f"{len(vector)} != {len(feature_names)}"
            )
        rows[str(chrom)] = vector
    
    if not rows:
        raise ValueError(f"No chromosome features found: {path}")
    

    matrix=np.asarray(list(rows.values()), dtype=np.float32)
    if normalize:
        mean = matrix.mean(axis=0, keepdims=True)
        std = matrix.std(axis=0, keepdims=True)
        std = np.where(std < 1e-6, 1.0, std)
        matrix = (matrix - mean)/std

    normalized_rows = {}
    for chrom, vector in zip(row.keys(), matrix):
        for alias in _chrom_aliases(chrom):
            normalized_rows[alias] = vector.astype(np.float32).tolist()

    return {
        "feature_names": list(feature_names),
        "features": normalized_rows,
        "dim": len(feature_names),
    }