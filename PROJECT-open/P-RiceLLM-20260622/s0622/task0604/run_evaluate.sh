#!/usr/bin/env bash
set -euo pipefail
# csv2bw2.py
# calc_metrics_for_batch_bw2.py
# csv_seg_eval.py

# PREDICT_DIR=/mnt/rice/default/Workspace/yangdong/gene_expression_prediction/outputs/predict/202606010533
PREDICT_DIR=/mnt/rice/default/Workspace/yangdong/gene_expression_prediction/outputs/predict/tmp
REF_DIR=/mnt/rice/default/Workspace/Rice-Genome/application/RNAseq/riceRNAseqData/18k/ref
EVAL_SCRIPTS_DIR="/mnt/rice/default/Workspace/yangdong/gene_expression_prediction/scripts/evaluation"

shopt -s nullglob
predict_dirs=("$PREDICT_DIR"/*/*/)
printf '%s\n' "${predict_dirs[@]}"

for CSV_DIR in "${predict_dirs[@]}"; do
    echo "CSV_DIR: $CSV_DIR"
    parent_dir=$(basename "$(dirname "$CSV_DIR")")
    SPECIES=$(echo "$parent_dir" | awk -F'_' '{print $(NF-1)}')
    CHROM=$(basename "$CSV_DIR")
    echo "CHROM=$CHROM, SPECIES=$SPECIES"

    FASTA="${REF_DIR}/${SPECIES}-new.fasta"
    CHROM_SIZES="${REF_DIR}/${SPECIES}-new.chrom.sizes"

    shopt -s nullglob

    cd "$CSV_DIR"

    csvs=(./*.csv) && printf '%s\n' "${csvs[@]}"

    if [ ${#csvs[@]} -eq 0 ]; then
        echo "未找到 CSV 文件：$CSV_DIR"
        exit 0
    fi

    for csv in "${csvs[@]}"; do
        out_txt="${csv%.csv}"
        echo "Convert: $csv -> $out_txt"
        python ${EVAL_SCRIPTS_DIR}/csv2bw2.py \
            --csv "$csv" \
            --output "$out_txt".bw \
            --chrom_sizes "$CHROM_SIZES" \
            --expression_col "predicted_expression"
        python ${EVAL_SCRIPTS_DIR}/calc_metrics_for_batch_bw2.py \
            --pred_files "$out_txt".bw.npy \
            --raw_files "$out_txt".bw_true.npy \
            --output "$out_txt".bw_track-level_stats.txt \
            --fasta "$FASTA" \
            --chrom "$CHROM"
        python ${EVAL_SCRIPTS_DIR}/csv_seg_eval.py \
            --csv "$csv" \
            --output $out_txt --skip_bigwig \
            --chrom_sizes "$CHROM_SIZES" \
            --expression_col "predicted_expression"
    done
done

python ${EVAL_SCRIPTS_DIR}/evaluate_feature_level_metrics.py \
"${PREDICT_DIR}" "$REF_DIR"

rm -f ${PREDICT_DIR}/merged.csv
# shopt -s nullglob; sum_list=(); for d in "${predict_dirs[@]}"; do sum_list+=( "$d"/*chromosome_stats.csv ); done
# head -n1 "${sum_list[0]}" > ${PREDICT_DIR}/merged.csv && for f in "${sum_list[@]}"; do tail -n +2 "$f"; done >> ${PREDICT_DIR}/merged.csv

shopt -s nullglob
sum_list=()
for d in "${predict_dirs[@]}"; do
    sum_list+=( "$d"/*chromosome_stats.csv )
done

head -n1 "${sum_list[0]}" | awk '{print $0",type"}' > "${PREDICT_DIR}/merged.csv"

for f in "${sum_list[@]}"; do
    type="$(basename "$(dirname "$(dirname "$f")")")"
    type="${type%%_*}"
    awk -v t="$type" 'NR>1{print $0","t}' "$f"
done >> "${PREDICT_DIR}/merged.csv"
