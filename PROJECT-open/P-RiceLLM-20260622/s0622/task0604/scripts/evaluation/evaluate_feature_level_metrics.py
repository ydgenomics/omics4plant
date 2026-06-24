#!/usr/bin/env python3
"""
Evaluate prediction agreement at gene and exon level.

The script only needs a prediction run directory and a GFF/GTF root. It scans
all *_predictions.csv files, matches each sample directory to P-number GFF
files, aggregates overlapping base-level predictions over gene/exon intervals,
and writes per-feature metrics plus summary tables.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm


DEFAULT_PREDICT_RUN_DIR = Path(
    "/mnt/rice/default/Workspace/xuxiaolong/RNAprediction/zhongLY/outputs/predict/202605140209-zeroloss"
)
DEFAULT_GFF_ROOT = Path(
    "/mnt/rice/default/Workspace/Rice-Genome/application/RNAseq/riceRNAseqData/18k/ref"
)
DUPLICATE_FEATURE_KEY_COLUMNS = [
    "sample_id",
    "annotation_path",
    "chromosome",
    "feature_type",
    "feature_id",
    "duplicate_count",
    "parent_ids",
    "strands",
    "intervals",
]


@dataclass(frozen=True)
class Feature:
    chrom: str
    start0: int
    end0: int
    eval_start0: int
    eval_end0: int
    feature_type: str
    feature_id: str
    parent_id: str
    strand: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute gene/exon-level Pearson, R2, and related metrics from prediction CSV files."
    )
    parser.add_argument(
        "predict_run_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_PREDICT_RUN_DIR,
        help=f"Prediction run directory. Default: {DEFAULT_PREDICT_RUN_DIR}",
    )
    parser.add_argument(
        "gff_root",
        nargs="?",
        type=Path,
        default=DEFAULT_GFF_ROOT,
        help=f"GFF/GTF root directory. Default: {DEFAULT_GFF_ROOT}",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=None,
        help="Output directory. Default: <predict_run_dir>/feature_level_metrics",
    )
    parser.add_argument(
        "--features",
        nargs="+",
        default=["gene", "exon"],
        help="Feature types to evaluate. Default: gene exon",
    )
    parser.add_argument(
        "--min_overlap_bp",
        type=int,
        default=1,
        help="Minimum overlapped predicted bases required for a feature. Default: 1",
    )
    parser.add_argument(
        "--min_nonzero_bp",
        type=int,
        default=2,
        help="Minimum bases with both pred and true > 0 for Pearson/Spearman. Default: 2",
    )
    parser.add_argument(
        "--min_r2_variance",
        type=float,
        default=1e-8,
        help="Minimum true-signal sum of squares required to report R2. Default: 1e-8",
    )
    parser.add_argument(
        "--feature_flank_bp",
        type=int,
        default=0,
        help="Extend each feature evaluation interval by this many bp on both sides. Default: 0",
    )
    parser.add_argument(
        "--feature_window_bp",
        type=int,
        default=None,
        help=(
            "Evaluate a fixed-width window centered on each feature midpoint. "
            "Applied after flank expansion if set. Default: use feature +/- flank."
        ),
    )
    parser.add_argument(
        "--pearson_delta_normalize",
        choices=["none", "global_mean", "nonzero_mean"],
        default="global_mean",
        help=(
            "Normalization for Pearson-delta-style pred vs true comparison. "
            "Default follows pearson_delta.py: global_mean."
        ),
    )
    parser.add_argument(
        "--pearson_delta_ref_mode",
        choices=["feature_normalized_mean", "true_global_mean", "zero", "pred_true_pooled_mean"],
        default="feature_normalized_mean",
        help=(
            "Reference for Pearson-delta-style feature-level comparison. "
            "feature_normalized_mean uses each feature's mean normalized true expression across training samples; "
            "true_global_mean uses the current tissue/track global mean of true feature values; "
            "zero reproduces the old implicit zero-ref behavior. Default: feature_normalized_mean."
        ),
    )
    parser.add_argument(
        "--strict_full_feature",
        action="store_true",
        help="Keep only features whose evaluation interval is fully covered by prediction windows.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute outputs even if per-CSV feature metric files already exist.",
    )
    parser.add_argument(
        "--csv_glob",
        default="*_multitrack/*/*_predictions.csv",
        help="Glob under predict_run_dir for prediction CSV files. Default: *_multitrack/*/*_predictions.csv",
    )
    parser.add_argument(
        "--max_csv",
        type=int,
        default=None,
        help="Optional maximum number of CSV files to process, useful for smoke tests.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes for per-CSV aggregation. Default: 1",
    )
    parser.add_argument(
        "--summary_only_from_existing",
        action="store_true",
        help=(
            "Skip prediction CSV aggregation and rebuild summary files from existing "
            "*_feature_metrics.csv files under output_dir."
        ),
    )
    return parser.parse_args()


def parse_attr(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def pick_id(attrs: dict[str, str], feature_type: str) -> str:
    for key in ("ID", "gene_id", "transcript_id", "Name"):
        value = attrs.get(key)
        if value:
            return value
    return f"{feature_type}:unknown"


def pick_parent(attrs: dict[str, str]) -> str:
    for key in ("Parent", "gene_id", "transcript_id"):
        value = attrs.get(key)
        if value:
            return value.split(",")[0]
    return ""


def sample_id_from_dir(sample_name: str) -> str | None:
    match = re.search(r"(?:^|_)P(\d+)(?:_|$)", sample_name)
    if not match:
        return None
    return f"P{int(match.group(1))}"


def find_annotation(gff_root: Path, sample_id: str) -> Path:
    candidates = [
        gff_root / f"{sample_id}.IRGSP.liftoff.gff3",
        gff_root / f"{sample_id}.MH63.liftoff.gff3",
        gff_root / f"{sample_id}.liftoff.gff3",
        gff_root / f"{sample_id}_EVM.all.gff3",
        gff_root / f"{sample_id}_EVM.all.gtf",
        gff_root / f"{sample_id}.IRGSP.liftoff.gff",
        gff_root / f"{sample_id}.MH63.liftoff.gff",
        gff_root / sample_id / f"{sample_id}_EVM.all.gff3",
        gff_root / sample_id / f"{sample_id}_EVM.all.gtf",
        gff_root / sample_id / f"{sample_id}.IRGSP.liftoff.gff3",
        gff_root / sample_id / f"{sample_id}.MH63.liftoff.gff3",
        gff_root / sample_id / f"{sample_id}.liftoff.gff3",
    ]
    for path in candidates:
        if path.is_file():
            return path
    globbed = (
        sorted(gff_root.glob(f"**/{sample_id}*liftoff*.gff3"))
        + sorted(gff_root.glob(f"**/{sample_id}*liftoff*.gff"))
        + sorted(gff_root.glob(f"**/{sample_id}*EVM*.gff3"))
        + sorted(gff_root.glob(f"**/{sample_id}*EVM*.gtf"))
    )
    if globbed:
        return globbed[0]
    raise FileNotFoundError(f"No GFF/GTF found for {sample_id} under {gff_root}")


def build_eval_interval(start0: int, end0: int, flank_bp: int, window_bp: int | None) -> tuple[int, int]:
    eval_start0 = max(0, start0 - flank_bp)
    eval_end0 = end0 + flank_bp
    if window_bp is not None:
        if window_bp <= 0:
            raise ValueError("--feature_window_bp must be a positive integer")
        midpoint = (start0 + end0) // 2
        half = window_bp // 2
        eval_start0 = max(0, midpoint - half)
        eval_end0 = eval_start0 + window_bp
    if eval_end0 <= eval_start0:
        eval_end0 = eval_start0 + 1
    return eval_start0, eval_end0


def load_features(
    annotation_path: Path,
    feature_types: set[str],
    flank_bp: int = 0,
    window_bp: int | None = None,
) -> dict[str, list[Feature]]:
    features_by_chrom: dict[str, list[Feature]] = defaultdict(list)
    with annotation_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, attr_text = parts[:9]
            if feature_type not in feature_types:
                continue
            attrs = parse_attr(attr_text)
            start0 = int(start) - 1
            end0 = int(end)
            if start0 < 0 or end0 <= start0:
                continue
            eval_start0, eval_end0 = build_eval_interval(start0, end0, flank_bp, window_bp)
            feature = Feature(
                chrom=chrom,
                start0=start0,
                end0=end0,
                eval_start0=eval_start0,
                eval_end0=eval_end0,
                feature_type=feature_type,
                feature_id=pick_id(attrs, feature_type),
                parent_id=pick_parent(attrs),
                strand=strand,
            )
            features_by_chrom[chrom].append(feature)

    for chrom in features_by_chrom:
        features_by_chrom[chrom].sort(key=lambda feature: (feature.eval_start0, feature.eval_end0))
    return dict(features_by_chrom)


def find_duplicate_feature_keys(
    sample_id: str,
    annotation_path: Path,
    features_by_chrom: dict[str, list[Feature]],
) -> list[dict[str, object]]:
    duplicate_rows = []
    for chrom, features in features_by_chrom.items():
        by_key: dict[tuple[str, str], list[Feature]] = defaultdict(list)
        for feature in features:
            by_key[(feature.feature_type, feature.feature_id)].append(feature)
        for (feature_type, feature_id), duplicates in by_key.items():
            if len(duplicates) <= 1:
                continue
            duplicate_rows.append(
                {
                    "sample_id": sample_id,
                    "annotation_path": str(annotation_path),
                    "chromosome": chrom,
                    "feature_type": feature_type,
                    "feature_id": feature_id,
                    "duplicate_count": len(duplicates),
                    "parent_ids": ";".join(sorted({feature.parent_id for feature in duplicates})),
                    "strands": ";".join(sorted({feature.strand for feature in duplicates})),
                    "intervals": ";".join(
                        f"{feature.start0}-{feature.end0}" for feature in duplicates
                    ),
                }
            )
    return duplicate_rows


def check_duplicate_feature_keys(
    csv_paths: list[Path],
    gff_root: Path,
    feature_types: set[str],
    flank_bp: int,
    window_bp: int | None,
    output_dir: Path,
    annotation_cache: dict[str, dict[str, list[Feature]]] | None = None,
) -> dict[str, dict[str, list[Feature]]]:
    annotation_cache = annotation_cache if annotation_cache is not None else {}
    sample_ids = sorted(
        {
            context["sample_id"]
            for csv_path in csv_paths
            for context in [parse_csv_context(csv_path)]
            if context["sample_id"]
        }
    )
    duplicate_rows = []
    for sample_id in sample_ids:
        annotation_path = find_annotation(gff_root, sample_id)
        if sample_id not in annotation_cache:
            annotation_cache[sample_id] = load_features(
                annotation_path,
                feature_types,
                flank_bp=flank_bp,
                window_bp=window_bp,
            )
        duplicate_rows.extend(
            find_duplicate_feature_keys(sample_id, annotation_path, annotation_cache[sample_id])
        )

    report_path = output_dir / "feature_level_duplicate_feature_keys.csv"
    pd.DataFrame(duplicate_rows, columns=DUPLICATE_FEATURE_KEY_COLUMNS).to_csv(
        report_path,
        index=False,
    )
    if duplicate_rows:
        print(
            f"WARNING: Found {len(duplicate_rows)} duplicate feature keys that will be merged "
            f"by aggregation. See: {report_path}"
        )
    else:
        print(f"No duplicate feature keys found. Report: {report_path}")
    return annotation_cache


def parse_array(value: object) -> np.ndarray:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return np.array([], dtype=np.float32)
    text = str(value).strip()
    if not text:
        return np.array([], dtype=np.float32)
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        parsed = ast.literal_eval(text)
    return np.asarray(parsed, dtype=np.float32)


def iter_prediction_rows(csv_path: Path) -> Iterable[tuple[str, int, int, np.ndarray, np.ndarray]]:
    for chunk in pd.read_csv(csv_path, chunksize=16):
        required = {"chromosome", "start", "end", "predicted_expression", "true_expression"}
        missing = required - set(chunk.columns)
        if missing:
            raise ValueError(f"{csv_path} missing required columns: {sorted(missing)}")
        for row in chunk.itertuples(index=False):
            chrom = str(getattr(row, "chromosome"))
            start = int(getattr(row, "start"))
            end = int(getattr(row, "end"))
            pred = parse_array(getattr(row, "predicted_expression"))
            true = parse_array(getattr(row, "true_expression"))
            length = end - start
            if length <= 0 or len(pred) != length or len(true) != length:
                continue
            yield chrom, start, end, pred, true


def discover_chromosomes(csv_path: Path) -> set[str]:
    chroms: set[str] = set()
    for chunk in pd.read_csv(csv_path, usecols=["chromosome"], chunksize=1000):
        chroms.update(chunk["chromosome"].astype(str).unique())
    return chroms


def features_for_csv(
    features_by_chrom: dict[str, list[Feature]],
    csv_path: Path,
) -> dict[str, list[Feature]]:
    chroms = discover_chromosomes(csv_path)
    return {chrom: features_by_chrom.get(chrom, []) for chrom in chroms}


def compute_basic_metrics(
    pred: np.ndarray,
    true: np.ndarray,
    min_nonzero_bp: int,
    min_r2_variance: float,
) -> dict[str, float]:
    metrics: dict[str, float] = {}
    if len(pred) == 0 or len(true) == 0:
        return {
            "pearson": np.nan,
            "spearman": np.nan,
            "r2": np.nan,
            "mse": np.nan,
            "mae": np.nan,
        }

    nonzero_mask = (pred > 0) & (true > 0)
    pred_nonzero = pred[nonzero_mask]
    true_nonzero = true[nonzero_mask]
    if (
        len(pred_nonzero) >= min_nonzero_bp
        and len(np.unique(pred_nonzero)) > 1
        and len(np.unique(true_nonzero)) > 1
    ):
        metrics["pearson"] = float(stats.pearsonr(pred_nonzero, true_nonzero).statistic)
        metrics["spearman"] = float(stats.spearmanr(pred_nonzero, true_nonzero).statistic)
    else:
        metrics["pearson"] = np.nan
        metrics["spearman"] = np.nan

    diff = true - pred
    metrics["mse"] = float(np.mean(diff**2))
    metrics["mae"] = float(np.mean(np.abs(diff)))
    ss_res = float(np.sum(diff**2))
    ss_tot = float(np.sum((true - np.mean(true)) ** 2))
    metrics["r2"] = float(1.0 - ss_res / ss_tot) if ss_tot > min_r2_variance else np.nan
    return metrics


def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2:
        return float("nan")
    if float(np.std(a)) <= 1e-8 or float(np.std(b)) <= 1e-8:
        return float("nan")
    return float(stats.pearsonr(a, b).statistic)


def pearson_delta_scale(arr: np.ndarray, method: str) -> float:
    values = np.asarray(arr, dtype=np.float64)
    if method == "none":
        return 1.0
    if method == "global_mean":
        mean_value = float(np.mean(values))
    elif method == "nonzero_mean":
        nonzero = values[values != 0.0]
        if nonzero.size == 0:
            return float("nan")
        mean_value = float(np.mean(nonzero))
    else:
        raise ValueError(f"Unknown pearson_delta normalization: {method}")
    if mean_value <= 0.0 or not np.isfinite(mean_value):
        return float("nan")
    return mean_value


def normalize_for_pearson_delta(arr: np.ndarray, method: str) -> tuple[np.ndarray, float]:
    scale = pearson_delta_scale(arr, method)
    values = np.asarray(arr, dtype=np.float64)
    if not np.isfinite(scale):
        return np.full_like(values, np.nan, dtype=np.float64), scale
    return values / scale, scale


def compute_pearson_delta_metrics(
    pred: np.ndarray,
    true: np.ndarray,
    normalize: str,
    ref_value: float | np.ndarray = 0.0,
) -> dict[str, float]:
    """Pearson-delta-compatible pred/true comparison with a scalar or feature-wise ref."""
    if len(pred) == 0 or len(true) == 0:
        return {
            "pearson_delta": np.nan,
            "delta_rmse": np.nan,
            "pearson_delta_ref_value": (
                float(ref_value) if np.isscalar(ref_value) else float("nan")
            ),
            "pearson_delta_scale_pred": np.nan,
            "pearson_delta_scale_true": np.nan,
            "pearson_delta_scale_ref": np.nan,
        }
    pred_norm, scale_pred = normalize_for_pearson_delta(pred, normalize)
    true_norm, scale_true = normalize_for_pearson_delta(true, normalize)
    if np.isscalar(ref_value):
        ref_value_for_report = float(ref_value)
        ref = np.full_like(np.asarray(true, dtype=np.float64), ref_value_for_report, dtype=np.float64)
        if ref_value_for_report == 0.0:
            ref_norm = ref
            scale_ref = 1.0
        else:
            ref_norm, scale_ref = normalize_for_pearson_delta(ref, normalize)
    else:
        ref = np.asarray(ref_value, dtype=np.float64)
        if ref.shape != np.asarray(true).shape:
            raise ValueError(f"ref shape mismatch: ref={ref.shape}, true={np.asarray(true).shape}")
        if normalize == "none":
            ref_norm = ref
            scale_ref = 1.0
        else:
            # A feature_normalized_mean ref is already in the same normalized space
            # as pred_norm/true_norm, so do not divide it by another track-level scale.
            ref_norm = ref
            scale_ref = 1.0
        finite_ref = ref[np.isfinite(ref)]
        ref_value_for_report = float(np.mean(finite_ref)) if finite_ref.size else float("nan")
    if np.isnan(pred_norm).any() or np.isnan(true_norm).any() or np.isnan(ref_norm).any():
        return {
            "pearson_delta": np.nan,
            "delta_rmse": np.nan,
            "pearson_delta_ref_value": ref_value_for_report,
            "pearson_delta_scale_pred": scale_pred,
            "pearson_delta_scale_true": scale_true,
            "pearson_delta_scale_ref": scale_ref,
        }
    delta_pred = pred_norm - ref_norm
    delta_true = true_norm - ref_norm
    diff = pred_norm - true_norm
    return {
        "pearson_delta": safe_pearson(delta_pred, delta_true),
        "delta_rmse": float(np.sqrt(np.mean(diff**2))),
        "pearson_delta_ref_value": ref_value_for_report,
        "pearson_delta_scale_pred": scale_pred,
        "pearson_delta_scale_true": scale_true,
        "pearson_delta_scale_ref": scale_ref,
    }


def resolve_pearson_delta_ref_value(
    pred: np.ndarray,
    true: np.ndarray,
    ref_mode: str,
) -> float:
    pred = np.asarray(pred, dtype=np.float64)
    true = np.asarray(true, dtype=np.float64)
    if ref_mode == "zero":
        return 0.0
    if ref_mode == "true_global_mean":
        return float(np.mean(true)) if true.size else float("nan")
    if ref_mode == "pred_true_pooled_mean":
        pooled = np.concatenate([pred, true]) if pred.size or true.size else np.array([], dtype=np.float64)
        return float(np.mean(pooled)) if pooled.size else float("nan")
    raise ValueError(f"Unknown pearson_delta ref mode: {ref_mode}")


def add_normalized_feature_values(
    df: pd.DataFrame,
    normalize: str,
) -> pd.DataFrame:
    """Add pred/true values normalized within each sample-track-feature_type group."""
    if df.empty:
        return df
    df = df.copy()
    group_cols = ["sample_dir", "chrom_unit", "biosample", "modality", "feature_type"]
    df["pred_mean_norm"] = np.nan
    df["true_mean_norm"] = np.nan
    df["pearson_delta_scale_pred"] = np.nan
    df["pearson_delta_scale_true"] = np.nan

    for _, idx in df.groupby(group_cols, dropna=False).groups.items():
        idx = list(idx)
        pred = df.loc[idx, "pred_mean"].to_numpy(dtype=float)
        true = df.loc[idx, "true_mean"].to_numpy(dtype=float)
        pred_norm, scale_pred = normalize_for_pearson_delta(pred, normalize)
        true_norm, scale_true = normalize_for_pearson_delta(true, normalize)
        df.loc[idx, "pred_mean_norm"] = pred_norm
        df.loc[idx, "true_mean_norm"] = true_norm
        df.loc[idx, "pearson_delta_scale_pred"] = scale_pred
        df.loc[idx, "pearson_delta_scale_true"] = scale_true
    return df


def add_feature_normalized_refs(df: pd.DataFrame) -> pd.DataFrame:
    """Use each feature's mean normalized true expression across training samples as ref."""
    if df.empty:
        return df
    df = df.copy()
    ref_group_cols = [
        "biosample",
        "modality",
        "feature_type",
        "chrom_unit",
        "parent_id",
        "feature_id",
    ]
    ref_source = df[df["split"] == "train"]
    ref = (
        ref_source.groupby(ref_group_cols, dropna=False)["true_mean_norm"]
        .mean()
        .rename("pearson_delta_ref_value")
        .reset_index()
    )
    df = df.drop(columns=["pearson_delta_ref_value"], errors="ignore")
    df = df.merge(ref, on=ref_group_cols, how="left")
    df["pearson_delta_ref_mode"] = "feature_normalized_mean"
    df["pearson_delta_scale_ref"] = 1.0
    return df


def finalize_pearson_delta_columns(
    df: pd.DataFrame,
    pearson_delta_ref_mode: str,
    pearson_delta_normalize: str,
) -> pd.DataFrame:
    if df.empty:
        return df
    df = add_normalized_feature_values(df, pearson_delta_normalize)
    if pearson_delta_ref_mode == "feature_normalized_mean":
        return add_feature_normalized_refs(df)

    df = df.copy()
    group_cols = ["sample_dir", "chrom_unit", "biosample", "modality", "feature_type"]
    df["pearson_delta_ref_mode"] = pearson_delta_ref_mode
    df["pearson_delta_ref_value"] = np.nan
    df["pearson_delta_scale_ref"] = np.nan
    for _, idx in df.groupby(group_cols, dropna=False).groups.items():
        idx = list(idx)
        pred = df.loc[idx, "pred_mean"].to_numpy(dtype=float)
        true = df.loc[idx, "true_mean"].to_numpy(dtype=float)
        ref_value = resolve_pearson_delta_ref_value(pred, true, pearson_delta_ref_mode)
        df.loc[idx, "pearson_delta_ref_value"] = ref_value
        if np.isfinite(ref_value):
            ref_for_scale = np.full(len(idx), ref_value, dtype=np.float64)
            df.loc[idx, "pearson_delta_scale_ref"] = pearson_delta_scale(
                ref_for_scale,
                pearson_delta_normalize,
            )
    return df


def summarize_across_features(
    df: pd.DataFrame,
    group_cols: list[str],
    min_nonzero_bp: int,
    pearson_delta_normalize: str,
    pearson_delta_ref_mode: str,
) -> pd.DataFrame:
    rows = []
    if df.empty:
        return pd.DataFrame()

    for keys, part in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        valid = part.dropna(subset=["pred_mean", "true_mean"])
        row["n_features"] = int(len(part))
        row["n_valid_features"] = int(len(valid))
        for col in (
            "overlap_bp",
            "feature_length",
            "eval_length",
            "pred_mean",
            "true_mean",
            "pred_sum",
            "true_sum",
            "pearson",
            "spearman",
            "r2",
            "mse",
            "mae",
        ):
            row[f"{col}_mean"] = float(valid[col].mean()) if col in valid and len(valid) else np.nan
            row[f"{col}_median"] = float(valid[col].median()) if col in valid and len(valid) else np.nan

        if len(valid) >= min_nonzero_bp:
            pred_mean = valid["pred_mean"].to_numpy(dtype=float)
            true_mean = valid["true_mean"].to_numpy(dtype=float)
            if len(np.unique(pred_mean)) > 1 and len(np.unique(true_mean)) > 1:
                row["feature_mean_pearson"] = float(stats.pearsonr(pred_mean, true_mean).statistic)
                row["feature_mean_spearman"] = float(stats.spearmanr(pred_mean, true_mean).statistic)
            else:
                row["feature_mean_pearson"] = np.nan
                row["feature_mean_spearman"] = np.nan
            ss_res = float(np.sum((true_mean - pred_mean) ** 2))
            ss_tot = float(np.sum((true_mean - np.mean(true_mean)) ** 2))
            row["feature_mean_r2"] = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
            if pearson_delta_ref_mode == "feature_normalized_mean":
                delta_valid = valid.dropna(
                    subset=["pred_mean_norm", "true_mean_norm", "pearson_delta_ref_value"]
                )
                finite_delta_mask = np.isfinite(
                    delta_valid[["pred_mean_norm", "true_mean_norm", "pearson_delta_ref_value"]]
                    .to_numpy(dtype=float)
                ).all(axis=1)
                delta_valid = delta_valid.loc[finite_delta_mask]
                if len(delta_valid) >= min_nonzero_bp:
                    pred_for_delta = delta_valid["pred_mean_norm"].to_numpy(dtype=float)
                    true_for_delta = delta_valid["true_mean_norm"].to_numpy(dtype=float)
                    ref_for_delta = delta_valid["pearson_delta_ref_value"].to_numpy(dtype=float)
                    delta_metrics = compute_pearson_delta_metrics(
                        pred_for_delta,
                        true_for_delta,
                        "none",
                        ref_value=ref_for_delta,
                    )
                else:
                    delta_metrics = {
                        "pearson_delta": np.nan,
                        "delta_rmse": np.nan,
                        "pearson_delta_ref_value": np.nan,
                        "pearson_delta_scale_pred": np.nan,
                        "pearson_delta_scale_true": np.nan,
                        "pearson_delta_scale_ref": np.nan,
                    }
            else:
                delta_metrics = compute_pearson_delta_metrics(
                    pred_mean,
                    true_mean,
                    pearson_delta_normalize,
                    ref_value=resolve_pearson_delta_ref_value(
                        pred_mean,
                        true_mean,
                        pearson_delta_ref_mode,
                    ),
                )
            row["feature_mean_pearson_delta"] = delta_metrics["pearson_delta"]
            row["feature_mean_delta_rmse"] = delta_metrics["delta_rmse"]
            row["feature_mean_pearson_delta_ref_mode"] = pearson_delta_ref_mode
            row["feature_mean_pearson_delta_ref_value"] = delta_metrics["pearson_delta_ref_value"]
            row["feature_mean_pearson_delta_scale_pred"] = delta_metrics["pearson_delta_scale_pred"]
            row["feature_mean_pearson_delta_scale_true"] = delta_metrics["pearson_delta_scale_true"]
            row["feature_mean_pearson_delta_scale_ref"] = delta_metrics["pearson_delta_scale_ref"]
        else:
            row["feature_mean_pearson"] = np.nan
            row["feature_mean_spearman"] = np.nan
            row["feature_mean_r2"] = np.nan
            row["feature_mean_pearson_delta"] = np.nan
            row["feature_mean_delta_rmse"] = np.nan
            row["feature_mean_pearson_delta_ref_mode"] = pearson_delta_ref_mode
            row["feature_mean_pearson_delta_ref_value"] = np.nan
            row["feature_mean_pearson_delta_scale_pred"] = np.nan
            row["feature_mean_pearson_delta_scale_true"] = np.nan
            row["feature_mean_pearson_delta_scale_ref"] = np.nan

        rows.append(row)
    return pd.DataFrame(rows)


def summarize_summary_rows(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    if df.empty:
        return pd.DataFrame()
    metric_cols = [
        col for col in df.columns
        if col not in set(group_cols)
        and pd.api.types.is_numeric_dtype(df[col])
    ]
    for keys, part in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["n_rows"] = int(len(part))
        for col in metric_cols:
            row[f"{col}_mean"] = float(part[col].mean()) if len(part) else np.nan
            row[f"{col}_median"] = float(part[col].median()) if len(part) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_csv_over_features(
    csv_path: Path,
    features_by_chrom: dict[str, list[Feature]],
    min_overlap_bp: int,
    min_nonzero_bp: int,
    min_r2_variance: float,
    pearson_delta_normalize: str,
    strict_full_feature: bool,
    show_progress: bool = True,
) -> pd.DataFrame:
    feature_values: dict[tuple[str, str, str], dict[str, object]] = {}
    relevant_features = features_for_csv(features_by_chrom, csv_path)

    for chrom, window_start, window_end, pred, true in tqdm(
        iter_prediction_rows(csv_path),
        desc=f"Aggregate {csv_path.parent.parent.name}/{csv_path.parent.name}/{csv_path.name}",
        unit="window",
        disable=not show_progress,
    ):
        chrom_features = relevant_features.get(chrom)
        if not chrom_features:
            continue
        for feature in chrom_features:
            if feature.eval_end0 <= window_start:
                continue
            if feature.eval_start0 >= window_end:
                break
            overlap_start = max(window_start, feature.eval_start0)
            overlap_end = min(window_end, feature.eval_end0)
            if overlap_end <= overlap_start:
                continue
            src_start = overlap_start - window_start
            src_end = overlap_end - window_start
            key = (feature.feature_type, feature.feature_id, feature.chrom)
            entry = feature_values.setdefault(
                key,
                {
                    "feature": feature,
                    "pred_sum": defaultdict(float),
                    "true_sum": defaultdict(float),
                    "counts": defaultdict(int),
                },
            )
            pred_slice = pred[src_start:src_end]
            true_slice = true[src_start:src_end]
            for offset, pos in enumerate(range(overlap_start, overlap_end)):
                entry["pred_sum"][pos] += float(pred_slice[offset])
                entry["true_sum"][pos] += float(true_slice[offset])
                entry["counts"][pos] += 1

    rows = []
    for entry in feature_values.values():
        feature = entry["feature"]
        covered_positions = sorted(entry["counts"].keys())
        overlap_bp = len(covered_positions)
        feature_length = feature.end0 - feature.start0
        eval_length = feature.eval_end0 - feature.eval_start0
        if overlap_bp < min_overlap_bp:
            continue
        if strict_full_feature and overlap_bp < eval_length:
            continue

        pred_values = np.asarray(
            [entry["pred_sum"][pos] / entry["counts"][pos] for pos in covered_positions],
            dtype=np.float32,
        )
        true_values = np.asarray(
            [entry["true_sum"][pos] / entry["counts"][pos] for pos in covered_positions],
            dtype=np.float32,
        )
        metrics = compute_basic_metrics(pred_values, true_values, min_nonzero_bp, min_r2_variance)
        nonzero_mask = (pred_values > 0) & (true_values > 0)
        rows.append(
            {
                "feature_type": feature.feature_type,
                "feature_id": feature.feature_id,
                "parent_id": feature.parent_id,
                "chromosome": feature.chrom,
                "start": feature.start0,
                "end": feature.end0,
                "eval_start": feature.eval_start0,
                "eval_end": feature.eval_end0,
                "strand": feature.strand,
                "feature_length": feature_length,
                "eval_length": eval_length,
                "overlap_bp": overlap_bp,
                "coverage_fraction": overlap_bp / eval_length if eval_length else np.nan,
                "pred_sum": float(np.sum(pred_values)),
                "true_sum": float(np.sum(true_values)),
                "pred_mean": float(np.mean(pred_values)),
                "true_mean": float(np.mean(true_values)),
                "pred_max": float(np.max(pred_values)),
                "true_max": float(np.max(true_values)),
                "pred_zero_ratio": float(np.mean(pred_values == 0)),
                "true_zero_ratio": float(np.mean(true_values == 0)),
                "nonzero_bp": int(np.sum(nonzero_mask)),
                **metrics,
            }
        )

    return pd.DataFrame(rows)


def process_csv_for_feature_metrics(
    csv_path: Path,
    gff_root: Path,
    feature_types: set[str],
    feature_flank_bp: int,
    feature_window_bp: int | None,
    min_overlap_bp: int,
    min_nonzero_bp: int,
    min_r2_variance: float,
    pearson_delta_normalize: str,
    pearson_delta_ref_mode: str,
    strict_full_feature: bool,
    output_dir: Path,
    show_progress: bool = False,
) -> tuple[pd.DataFrame, dict[str, object]]:
    context = parse_csv_context(csv_path)
    sample_id = context["sample_id"]
    if not sample_id:
        raise ValueError(f"cannot parse sample id from {context['sample_dir']}: {csv_path}")

    annotation_path = find_annotation(gff_root, sample_id)
    features_by_chrom = load_features(
        annotation_path,
        feature_types,
        flank_bp=feature_flank_bp,
        window_bp=feature_window_bp,
    )
    feature_metrics = aggregate_csv_over_features(
        csv_path=csv_path,
        features_by_chrom=features_by_chrom,
        min_overlap_bp=min_overlap_bp,
        min_nonzero_bp=min_nonzero_bp,
        min_r2_variance=min_r2_variance,
        pearson_delta_normalize=pearson_delta_normalize,
        strict_full_feature=strict_full_feature,
        show_progress=show_progress,
    )
    if not feature_metrics.empty:
        for key, value in context.items():
            feature_metrics.insert(0, key, value)

    detail_path = output_dir / context["sample_dir"] / context["chrom_unit"] / f"{safe_stem(csv_path)}_feature_metrics.csv"
    manifest_row = {
        **context,
        "csv_path": str(csv_path),
        "annotation_path": str(annotation_path),
        "feature_flank_bp": feature_flank_bp,
        "feature_window_bp": feature_window_bp,
        "pearson_delta_normalize": pearson_delta_normalize,
        "pearson_delta_ref_mode": pearson_delta_ref_mode,
        "n_feature_rows": int(len(feature_metrics)),
        "detail_path": str(detail_path),
    }
    return feature_metrics, manifest_row


def parse_csv_context(csv_path: Path) -> dict[str, str]:
    sample_name = csv_path.parent.parent.name
    chrom_unit = csv_path.parent.name
    match = re.match(r"(.+)__(.+)_predictions\.csv$", csv_path.name)
    biosample = match.group(1) if match else ""
    modality = match.group(2) if match else ""
    sample_id = sample_id_from_dir(sample_name) or ""
    split = sample_name.split("_", 1)[0] if "_" in sample_name else sample_name
    return {
        "split": split,
        "sample_dir": sample_name,
        "sample_id": sample_id,
        "chrom_unit": chrom_unit,
        "biosample": biosample,
        "modality": modality,
    }


def safe_stem(path: Path) -> str:
    return re.sub(r"[^0-9A-Za-z_.+-]+", "_", path.stem)


def write_outputs_for_csv(
    csv_path: Path,
    feature_metrics: pd.DataFrame,
    output_dir: Path,
    min_nonzero_bp: int,
    pearson_delta_normalize: str,
    pearson_delta_ref_mode: str,
    context: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    context = context or parse_csv_context(csv_path)
    csv_out_dir = output_dir / context["sample_dir"] / context["chrom_unit"]
    csv_out_dir.mkdir(parents=True, exist_ok=True)

    if feature_metrics.empty:
        detail_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_metrics.csv"
        pd.DataFrame().to_csv(detail_path, index=False)
        return pd.DataFrame(), pd.DataFrame()

    for key, value in reversed(list(context.items())):
        if key not in feature_metrics.columns:
            feature_metrics.insert(0, key, value)

    detail_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_metrics.csv"
    feature_metrics.to_csv(detail_path, index=False)
    for feature_type, part in feature_metrics.groupby("feature_type", dropna=False):
        type_path = csv_out_dir / f"{safe_stem(csv_path)}_{feature_type}_level_metrics.csv"
        part.to_csv(type_path, index=False)

    per_type = summarize_across_features(
        feature_metrics,
        ["split", "sample_dir", "sample_id", "chrom_unit", "biosample", "modality", "feature_type"],
        min_nonzero_bp,
        pearson_delta_normalize,
        pearson_delta_ref_mode,
    )
    per_csv = summarize_across_features(
        feature_metrics,
        ["split", "sample_dir", "sample_id", "chrom_unit", "biosample", "modality"],
        min_nonzero_bp,
        pearson_delta_normalize,
        pearson_delta_ref_mode,
    )
    per_type_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_summary_by_type.csv"
    per_csv_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_summary.csv"
    per_type.to_csv(per_type_path, index=False)
    per_csv.to_csv(per_csv_path, index=False)
    return feature_metrics, per_type


def feature_metrics_path_matches_output(path: Path) -> bool:
    name = path.name
    return (
        name.endswith("_feature_metrics.csv")
        and not name.endswith("_level_metrics.csv")
        and "_feature_summary" not in name
    )


def load_existing_feature_metrics(output_dir: Path) -> list[pd.DataFrame]:
    paths = sorted(
        path
        for path in output_dir.glob("*/*/*_feature_metrics.csv")
        if feature_metrics_path_matches_output(path)
    )
    if not paths:
        raise FileNotFoundError(f"No existing *_feature_metrics.csv files found under {output_dir}")

    frames = []
    for path in tqdm(paths, desc="Load feature metrics", unit="file"):
        df = pd.read_csv(path)
        if not df.empty:
            frames.append(df)
    if not frames:
        raise ValueError(f"Existing feature metric files under {output_dir} are all empty")
    return frames


def write_summary_outputs(
    all_feature_metrics: list[pd.DataFrame],
    csv_paths: list[Path],
    output_dir: Path,
    min_nonzero_bp: int,
    pearson_delta_ref_mode: str,
    pearson_delta_normalize: str,
) -> None:
    if not all_feature_metrics:
        return

    all_features_df = pd.concat(all_feature_metrics, ignore_index=True)
    all_features_df = finalize_pearson_delta_columns(
        all_features_df,
        pearson_delta_ref_mode,
        pearson_delta_normalize,
    )
    csv_context_by_key = {
        (
            context["sample_dir"],
            context["chrom_unit"],
            context["biosample"],
            context["modality"],
        ): (csv_path, context)
        for csv_path in csv_paths
        for context in [parse_csv_context(csv_path)]
    }

    all_summaries = []
    for csv_key, feature_metrics in all_features_df.groupby(
        ["sample_dir", "chrom_unit", "biosample", "modality"],
        dropna=False,
    ):
        sample_dir, chrom_unit, biosample, modality = csv_key
        csv_path, context = csv_context_by_key.get(
            csv_key,
            (
                output_dir / sample_dir / chrom_unit / f"{biosample}__{modality}_predictions.csv",
                {
                    "split": str(feature_metrics["split"].iloc[0]),
                    "sample_dir": sample_dir,
                    "sample_id": str(feature_metrics["sample_id"].iloc[0]),
                    "chrom_unit": chrom_unit,
                    "biosample": biosample,
                    "modality": modality,
                },
            ),
        )
        _, summary = write_outputs_for_csv(
            csv_path=csv_path,
            feature_metrics=feature_metrics.copy(),
            output_dir=output_dir,
            min_nonzero_bp=min_nonzero_bp,
            pearson_delta_normalize=pearson_delta_normalize,
            pearson_delta_ref_mode=pearson_delta_ref_mode,
            context=context,
        )
        if not summary.empty:
            all_summaries.append(summary)

    if not all_summaries:
        return
    summary_df = pd.concat(all_summaries, ignore_index=True)
    summary_df.to_csv(output_dir / "feature_level_summary_by_csv_and_type.csv", index=False)
    run_summary = summarize_summary_rows(summary_df, ["feature_type"])
    run_summary.to_csv(output_dir / "feature_level_run_summary_by_type.csv", index=False)


def main() -> None:
    args = parse_args()
    predict_run_dir = args.predict_run_dir.resolve()
    gff_root = args.gff_root.resolve()
    output_dir = (args.output_dir or predict_run_dir / "feature_level_metrics").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not predict_run_dir.is_dir():
        raise FileNotFoundError(f"Prediction run directory not found: {predict_run_dir}")
    if not gff_root.is_dir():
        raise FileNotFoundError(f"GFF root directory not found: {gff_root}")

    csv_paths = sorted(predict_run_dir.glob(args.csv_glob))
    if args.max_csv is not None:
        csv_paths = csv_paths[: args.max_csv]
    if not csv_paths:
        raise FileNotFoundError(f"No *_predictions.csv found under {predict_run_dir}")

    feature_types = set(args.features)
    annotation_cache: dict[str, dict[str, list[Feature]]] = {}
    all_feature_metrics = []
    manifest_rows = []

    print(f"Prediction run: {predict_run_dir}")
    print(f"GFF root: {gff_root}")
    print(f"Output dir: {output_dir}")
    print(f"Prediction CSV files: {len(csv_paths)}")
    print(f"Feature flank bp: {args.feature_flank_bp}")
    print(f"Feature window bp: {args.feature_window_bp}")
    print(f"Pearson-delta normalize: {args.pearson_delta_normalize}")
    print(f"Pearson-delta ref mode: {args.pearson_delta_ref_mode}")
    print(f"Workers: {args.workers}")

    workers = max(1, int(args.workers))
    check_cache = annotation_cache if workers == 1 else None
    check_duplicate_feature_keys(
        csv_paths,
        gff_root,
        feature_types,
        args.feature_flank_bp,
        args.feature_window_bp,
        output_dir,
        check_cache,
    )

    if args.summary_only_from_existing:
        all_feature_metrics = load_existing_feature_metrics(output_dir)
        write_summary_outputs(
            all_feature_metrics,
            csv_paths,
            output_dir,
            args.min_nonzero_bp,
            args.pearson_delta_ref_mode,
            args.pearson_delta_normalize,
        )
        print("Done.")
        return

    if workers == 1:
        for csv_path in csv_paths:
            context = parse_csv_context(csv_path)
            sample_id = context["sample_id"]
            if not sample_id:
                print(f"Skip {csv_path}: cannot parse sample id from {context['sample_dir']}")
                continue

            annotation_path = find_annotation(gff_root, sample_id)
            if sample_id not in annotation_cache:
                print(f"Load annotation for {sample_id}: {annotation_path}")
                annotation_cache[sample_id] = load_features(
                    annotation_path,
                    feature_types,
                    flank_bp=args.feature_flank_bp,
                    window_bp=args.feature_window_bp,
                )

            csv_out_dir = output_dir / context["sample_dir"] / context["chrom_unit"]
            detail_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_metrics.csv"
            if detail_path.exists() and not args.force:
                print(f"Existing output, skip: {detail_path}")
                if detail_path.exists():
                    all_feature_metrics.append(pd.read_csv(detail_path))
                continue

            feature_metrics = aggregate_csv_over_features(
                csv_path=csv_path,
                features_by_chrom=annotation_cache[sample_id],
                min_overlap_bp=args.min_overlap_bp,
                min_nonzero_bp=args.min_nonzero_bp,
                min_r2_variance=args.min_r2_variance,
                pearson_delta_normalize=args.pearson_delta_normalize,
                strict_full_feature=args.strict_full_feature,
                show_progress=True,
            )
            if not feature_metrics.empty:
                for key, value in context.items():
                    feature_metrics.insert(0, key, value)
                all_feature_metrics.append(feature_metrics)
            manifest_rows.append(
                {
                    **context,
                    "csv_path": str(csv_path),
                    "annotation_path": str(annotation_path),
                    "feature_flank_bp": args.feature_flank_bp,
                    "feature_window_bp": args.feature_window_bp,
                    "pearson_delta_normalize": args.pearson_delta_normalize,
                    "pearson_delta_ref_mode": args.pearson_delta_ref_mode,
                    "n_feature_rows": int(len(feature_metrics)),
                    "detail_path": str(detail_path),
                }
            )
    else:
        pending_csv_paths = []
        for csv_path in csv_paths:
            context = parse_csv_context(csv_path)
            if not context["sample_id"]:
                print(f"Skip {csv_path}: cannot parse sample id from {context['sample_dir']}")
                continue
            csv_out_dir = output_dir / context["sample_dir"] / context["chrom_unit"]
            detail_path = csv_out_dir / f"{safe_stem(csv_path)}_feature_metrics.csv"
            if detail_path.exists() and not args.force:
                print(f"Existing output, skip: {detail_path}")
                all_feature_metrics.append(pd.read_csv(detail_path))
                continue
            pending_csv_paths.append(csv_path)

        with ProcessPoolExecutor(max_workers=min(workers, len(pending_csv_paths) or 1)) as executor:
            futures = {
                executor.submit(
                    process_csv_for_feature_metrics,
                    csv_path,
                    gff_root,
                    feature_types,
                    args.feature_flank_bp,
                    args.feature_window_bp,
                    args.min_overlap_bp,
                    args.min_nonzero_bp,
                    args.min_r2_variance,
                    args.pearson_delta_normalize,
                    args.pearson_delta_ref_mode,
                    args.strict_full_feature,
                    output_dir,
                    False,
                ): csv_path
                for csv_path in pending_csv_paths
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="Aggregate CSVs", unit="csv"):
                feature_metrics, manifest_row = future.result()
                if not feature_metrics.empty:
                    all_feature_metrics.append(feature_metrics)
                manifest_rows.append(manifest_row)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(output_dir / "feature_level_manifest.csv", index=False)

    write_summary_outputs(
        all_feature_metrics,
        csv_paths,
        output_dir,
        args.min_nonzero_bp,
        args.pearson_delta_ref_mode,
        args.pearson_delta_normalize,
    )

    print("Done.")


if __name__ == "__main__":
    main()
