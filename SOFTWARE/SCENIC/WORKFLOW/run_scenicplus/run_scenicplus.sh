snakemake -s Snakefile --configfile config.yaml --cores 16 --show-failed-logs -h

mkdir -p /data/work/scenicplus/run_scenicplus/run_scenicplus  && cd /data/work/scenicplus/run_scenicplus/run_scenicplus 
snakemake -s ../Snakefile --configfile ../config.yaml --cores 16 --show-failed-logs

cd /data/work/scenicplus/run_scenicplus/run_scenicplus
scenicplus grn_inference motif_enrichment_cistarget \
--region_set_folder "/data/work/scenicplus/create_cistopic/get_cistopic_model/region_sets" \
--cistarget_db_fname "/data/work/scenicplus/create_cistarget_db/create_cistarget_motif_databases/Os.regions_vs_motifs.rankings.feather" \
--output_fname_cistarget_result "./Output/ctx_results.hdf5" \
--temp_dir "./Temporary" \
--species "Os" \
--fr_overlap_w_ctx_db 0.4 \
--auc_threshold 0.005 \
--nes_threshold 3.0 \
--rank_threshold 0.05 \
--path_to_motif_annotations "/data/work/scenicplus/input/from_planttfdb/process_motif/Os_TF_binding_motifs_information.tbl" \
--annotation_version "os" \
--motif_similarity_fdr 0.001 \
--orthologous_identity_threshold 0.0 \
--annotations_to_use "Direct_annot Orthology_annot" \
--write_html \
--output_fname_cistarget_html "./Output/ctx_results.html" \
--n_cpu 30
                
scenicplus prepare_data prepare_menr \
--paths_to_motif_enrichment_results ./Output/dem_result_fname.h5ad ./Output/ctx_results.hdf5 \
--multiome_mudata_fname ACC_GEX.h5mu \
--out_file_tf_names ./Output/tf_names.txt \
--out_file_direct_annotation ./Output/cistromes_direct.h5ad \
--out_file_extended_annotation ./Output/cistromes_extended.h5ad \
--direct_annotation Direct_annot \
--extended_annotation Orthology_annot

scenicplus grn_inference eGRN \
--TF_to_gene_adj_fname {input.tf_to_gene_adjacencies}\
--region_to_gene_adj_fname {input.region_to_gene_adjacencies} \
--cistromes_fname {input.cistromes_direct}\
--ranking_db_fname {input.ranking_db_fname} \
--eRegulon_out_fname {output.eRegulons_direct} \
--temp_dir {params.temp_dir} \
--order_regions_to_genes_by {params.order_regions_to_genes_by} \
--order_TFs_to_genes_by {params.order_TFs_to_genes_by} \
--gsea_n_perm {params.gsea_n_perm} \
--quantiles {params.quantiles} \
--top_n_regionTogenes_per_gene {params.top_n_regionTogenes_per_gene} \
--top_n_regionTogenes_per_region {params.top_n_regionTogenes_per_region} \
--min_regions_per_gene {params.min_regions_per_gene} \
--rho_threshold {params.rho_threshold} \
--min_target_genes {params.min_target_genes} \
--n_cpu {threads}