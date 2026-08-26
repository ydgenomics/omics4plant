cd /data/work/yita/Integration
input_rds='/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.rds'
out_rds='GM_SCTransform.CCA_integrated.rds'
out_UMAP='GM_SCTransform.CCA_integrated_UMAP.pdf'
batch_key='biosample'
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'

Rscript /data/work/Integration/SCTransform.CCA_integration.R \
--input_rds $input_rds --out_rds $out_rds --out_UMAP $out_UMAP \
--batch_key $batch_key --key_list $key_list --resolution $resolution \
--cluster_name $cluster_name


cd /data/work/yita/Integration
input_rds='/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.rds'
out_rds='GM_SCTransform.harmony_integrated.rds'
out_UMAP='GM_SCTransform.harmony_integrated_UMAP.pdf'
batch_key='biosample'
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'

Rscript /data/work/Integration/SCTransform.harmony_integration.R \
--input_rds $input_rds --out_rds $out_rds --out_UMAP $out_UMAP \
--batch_key $batch_key --key_list $key_list --resolution $resolution \
--cluster_name $cluster_name


cd /data/work/yita/Integration
input_rds='/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.rds'
out_rds='GM_BBKNNR_integrated.rds'
out_UMAP='GM_BBKNNR_integrated_UMAP.pdf'
batch_key='biosample'
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'

Rscript /data/work/Integration/BBKNNR_integration.R \
--input_rds $input_rds --out_rds $out_rds --out_UMAP $out_UMAP \
--batch_key $batch_key --key_list $key_list --resolution $resolution \
--cluster_name $cluster_name


cd /data/work/yita/Integration
input_rds='/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.rds'
out_rds='GM_rliger.INMF_integrated.rds'
out_UMAP='GM_rliger.INMF_integrated_UMAP.pdf'
batch_key='biosample'
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'

Rscript /data/work/Integration/rliger.INMF_integration.R \
--input_rds $input_rds --out_rds $out_rds --out_UMAP $out_UMAP \
--batch_key $batch_key --key_list $key_list --resolution $resolution \
--cluster_name $cluster_name


cd /data/work/yita/Integration
input_h5ad="/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.h5ad"
out_h5ad="GM_harmony_integrated.h5ad"
out_UMAP="GM_harmony_integrated_UMAP.pdf"
batch_key="biosample"
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'
/usr/bin/python /data/work/Integration/harmony_integration.py \
$input_h5ad $out_h5ad $out_UMAP --batch_key biosample --key_list $key_list \
--cluster_name $cluster_name --resolution $resolution


cd /data/work/yita/Integration
input_h5ad="/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.h5ad"
out_h5ad="GM_scVI_integrated.h5ad"
out_UMAP="GM_scVI_integrated_UMAP.pdf"
batch_key="biosample"
key_list='biosample,sample,metaneighbor'
resolution=0.5
cluster_name='celltype'
/usr/bin/python /data/work/Integration/scVI_integration.py \
$input_h5ad $out_h5ad $out_UMAP --batch_key biosample --key_list $key_list \
--cluster_name $cluster_name --resolution $resolution


cd /data/work/yita/Integration/scib-metrics
unintegrated_h5ad="/data/work/yita/MetaNeighbor/GM_meta_metaNeighbor.h5ad"
integrated_file="/data/work/yita/Integration/harmony_SCTransform.harmony_integrated.csv,/data/work/yita/Integration/iNMF_rliger.INMF_integrated.csv,/data/work/yita/Integration/pca_SCTransform.CCA_integrated.csv,/data/work/yita/Integration/X_pca_harmony_harmony_integrated.csv,/data/work/yita/Integration/X_scVI_scVI_integrated.csv"
tests_file="true,true,true,true,true,true,true,true,true,true"
batch_key="biosample"
label_key="metaneighbor"
n_jobs=16
prefix="GM_integrate"
/software/conda/Anaconda/bin/python /data/work/Integration/scIB.py \
--unintegrated_h5ad $unintegrated_h5ad --integrated_file $integrated_file \
--tests_file $tests_file --batch_key $batch_key --label_key $label_key \
--n_jobs $n_jobs --prefix $prefix


