# [3d-OT](https://github.com/dbjzs/3d-OT)

domain identity


```shell
git clone https://github.com/dbjzs/3d-OT.git && cd 3d-OT
source /opt/software/miniconda3/bin/activate
mamba create -n 3d-OT -c conda-forge python==3.10 libopenblas r-base=4.3 r-mclust -y && conda activate 3d-OT
mamba install anndata numpy pandas scanpy scikit-learn scipy \
scikit-misc tqdm rpy2 plotly nbformat matplotlib-inline matplotlib ipykernel -y
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124
pip install torch_geometric
```