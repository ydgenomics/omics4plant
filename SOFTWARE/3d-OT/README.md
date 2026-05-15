domain identity


```shell
source /opt/software/miniconda3/bin/activate
git clone https://github.com/dbjzs/3d-OT.git && cd 3d-OT
mamba create -n 3d-OT -c conda-forge python==3.10 libopenblas r-base=4.3 r-mclust -y && conda activate 3d-OT
mamba install anndata numpy pandas scanpy torch torch_geometric scikit-learn scipy \
scikit-misc tqdm rpy2 plotly nbformat matplotlib-inline matplotlib ipykernel -y
```