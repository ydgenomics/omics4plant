#!/bin/bash

set -euo pipefail

CHROMOSOMES=(Chr06 Chr09)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${SCRIPT_DIR}/"
PREDICT_PY="${SCRIPT_DIR}/scripts/predict_zeroloss.py"
OUTPUT_ROOT="${SCRIPT_DIR}/outputs/predict"
TRAIN_ROOT="${SCRIPT_DIR}/outputs/train"
DATA_INDICES_DIR="${SCRIPT_DIR}/data/indices"

mkdir -p "${OUTPUT_ROOT}"

if [[ ! -f "${PREDICT_PY}" ]]; then
    echo "Missing predict.py: ${PREDICT_PY}" >&2
    exit 1
fi

if [[ ! -d "${TRAIN_ROOT}" ]]; then
    echo "Missing train output directory: ${TRAIN_ROOT}" >&2
    exit 1
fi

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

export CUDA_DEVICE_MAX_CONNECTIONS=1
export NVTE_DEBUG=1
export NVTE_DEBUG_LEVEL=2
export NVTE_COMM_OVERLAP=0
export NCCL_P2P_DISABLE=1
export NCCL_P2P_DIRECT_DISABLE=1
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export TOKENIZERS_PARALLELISM=false

detect_gpu_count() {
    if [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]]; then
        IFS=',' read -r -a gpu_ids <<< "${CUDA_VISIBLE_DEVICES}"
        echo "${#gpu_ids[@]}"
        return
    fi

    if command -v nvidia-smi >/dev/null 2>&1; then
        nvidia-smi --query-gpu=index --format=csv,noheader | wc -l
        return
    fi

    echo 0
}

GPU_COUNT="$(detect_gpu_count)"
if [[ "${GPU_COUNT}" -lt 1 ]]; then
    echo "No GPU detected for prediction." >&2
    exit 1
fi

if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    CUDA_VISIBLE_DEVICES="$(seq -s, 0 $((GPU_COUNT - 1)))"
    export CUDA_VISIBLE_DEVICES
fi

export NPROC_PER_NODE="${NPROC_PER_NODE:-${GPU_COUNT}}"
if [[ "${NPROC_PER_NODE}" -gt "${GPU_COUNT}" ]]; then
    echo "NPROC_PER_NODE=${NPROC_PER_NODE} is greater than detected GPU count=${GPU_COUNT}" >&2
    exit 1
fi

DISTRIBUTED_ARGS=(
    --nnodes 1
    --nproc_per_node "${NPROC_PER_NODE}"
    --node_rank 0
    --master_addr localhost
    --master_port 29501
)

detect_latest_ckpt_path() {
    local latest_train_name latest_train_dir ckpt_name ckpt_path

    latest_train_name="$(
        find "${TRAIN_ROOT}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
            | sort \
            | tail -n 1
    )"

    if [[ -z "${latest_train_name}" ]]; then
        echo "No training folders found under ${TRAIN_ROOT}" >&2
        return 1
    fi

    latest_train_dir="${TRAIN_ROOT}/${latest_train_name}"
    ckpt_name="$(
        find "${latest_train_dir}" -mindepth 1 -maxdepth 1 -type d -name 'checkpoint-*' -printf '%f\n' \
            | awk -F- '$2 ~ /^[0-9]+$/ {print $2 "\t" $0}' \
            | sort -n \
            | tail -n 1 \
            | cut -f 2-
    )"

    if [[ -z "${ckpt_name}" ]]; then
        echo "No checkpoint folders found under ${latest_train_dir}" >&2
        return 1
    fi

    ckpt_path="${latest_train_dir}/${ckpt_name}/model.safetensors"
    if [[ ! -f "${ckpt_path}" ]]; then
        echo "Missing checkpoint model: ${ckpt_path}" >&2
        return 1
    fi

    echo "${ckpt_path}"
}

CKPT_PATH="$(detect_latest_ckpt_path)"
LATEST_TRAIN_NAME="$(basename "$(dirname "$(dirname "${CKPT_PATH}")")")"

cd "${WORK_ROOT}"

mapfile -t SAMPLE_DIRS < <(find "${DATA_INDICES_DIR}" -mindepth 1 -maxdepth 1 -type d | sort)
# if [[ "${#SAMPLE_DIRS[@]}" -ne 5 ]]; then
#     echo "Expected 5 sample directories under ${WORK_ROOT}/${DATA_INDICES_DIR}, found ${#SAMPLE_DIRS[@]}" >&2
#     exit 1
# fi


echo "Using checkpoint: ${CKPT_PATH}"
for sample_dir in "${SAMPLE_DIRS[@]}"; do
    sample_name="$(basename "${sample_dir}")"

    for required_file in sequence_split_train.csv index_stat.json bigWig_labels_meta.csv; do
        if [[ ! -f "${sample_dir}/${required_file}" ]]; then
            echo "Missing ${required_file} in sample directory: ${sample_dir}" >&2
            exit 1
        fi
    done

    for chr in "${CHROMOSOMES[@]}"; do
        torchrun "${DISTRIBUTED_ARGS[@]}" "${PREDICT_PY}" \
          --model_path /mnt/rice/default/Workspace/xuxiaolong/rice_1B_stage2_8k_hf \
          --tokenizer_path /mnt/rice/default/Workspace/xuxiaolong/rice_1B_stage2_8k_hf \
          --ckpt_path "${CKPT_PATH}" \
          --sequence_split_test "${sample_dir}/sequence_split_train.csv" \
          --index_stat_json "${sample_dir}/index_stat.json" \
          --bigWig_labels_meta "${sample_dir}/bigWig_labels_meta.csv" \
          --max_predict_samples 500 \
          --output_base_dir "${OUTPUT_ROOT}/${LATEST_TRAIN_NAME}/${sample_name}/${chr}" \
          --test_chromosomes "${chr}" \
          --batch_size 3 \
          --num_workers 4 \
          --use_flash_attn
    done
done
