#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${SCRIPT_DIR}/"
TRAIN_PY="${SCRIPT_DIR}/scripts/train_zeroloss.py"
OUTPUT_ROOT="${SCRIPT_DIR}/outputs/train"
TIMESTAMP="$(date +%Y%m%d%H%M)"
CHRS="Chr06"

mkdir -p "${OUTPUT_ROOT}"

if [[ ! -f "${TRAIN_PY}" ]]; then
    echo "Missing train.py: ${TRAIN_PY}" >&2
    exit 1
fi

# 自动检测当前节点上的可用 GPU 数量并设置 CUDA_VISIBLE_DEVICES

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

export WANDB_API_KEY=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
export WANDB_ENTITY="${WANDB_ENTITY:-ydgenomics2-bgi-group}"
export WANDB_PROJECT="${WANDB_PROJECT:-OneGenomeRice}"
export CUDA_DEVICE_MAX_CONNECTIONS=1
export NVTE_DEBUG=1
export NVTE_DEBUG_LEVEL=2
export NVTE_COMM_OVERLAP=0
export NCCL_P2P_DISABLE=1
export NCCL_P2P_DIRECT_DISABLE=1


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
    echo "No GPU detected for training." >&2
    exit 1
fi

if [[ -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
    CUDA_VISIBLE_DEVICES="$(seq -s, 0 $((GPU_COUNT - 1)))"
    export CUDA_VISIBLE_DEVICES
fi

export nproc_per_node="${nproc_per_node:-${GPU_COUNT}}"
if [[ "${nproc_per_node}" -gt "${GPU_COUNT}" ]]; then
    echo "nproc_per_node=${nproc_per_node} is greater than detected GPU count=${GPU_COUNT}" >&2
    exit 1
fi

USE_WANDB="${USE_WANDB:-1}"
WANDB_ARGS=()
if [[ "${USE_WANDB}" == "1" || "${USE_WANDB}" == "true" || "${USE_WANDB}" == "TRUE" ]]; then
    WANDB_ARGS+=(--use_wandb)
fi

DISTRIBUTED_ARGS=(
    --nnodes 1
    --nproc_per_node "${nproc_per_node}"
    --node_rank 0
    --master_addr localhost
    --master_port 29520
)

cd "${WORK_ROOT}"

# --index_stat_val_json

torchrun "${DISTRIBUTED_ARGS[@]}" "${TRAIN_PY}" \
        --model_path /mnt/rice/default/Workspace/xuxiaolong/guoyafei/rice_1B_stage2_8k_hf \
        --tokenizer_dir /mnt/rice/default/Workspace/xuxiaolong/guoyafei/rice_1B_stage2_8k_hf \
        --sequence_split_train_multi ./data/indices/train_*_multitrack/sequence_split_train.csv \
        --index_stat_multi_json ./data/indices/train_*_multitrack/index_stat.json \
        --sequence_split_val ./data/indices/valid_*_multitrack/sequence_split_train.csv \
        --index_stat_val_json ./data/indices/valid_*_multitrack/index_stat.json \
        --train_chromosomes $CHRS \
        --val_chromosomes $CHRS \
        --output_base_dir "${OUTPUT_ROOT}/${TIMESTAMP}" \
        --lr 0.00005 \
        --batch_size_per_device 1 \
        --gradient_accumulation_steps 10 \
        --num_train_epochs 20 \
        --loss_func mse \
        --max_sequence_length 32768 \
        --use_flash_attn \
        --gpus_per_node "${nproc_per_node}" \
        "${WANDB_ARGS[@]}" \
        "$@"



