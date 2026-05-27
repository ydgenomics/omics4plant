source /opt/software/miniconda3/bin/activate

mamba create -n st python=3.12 -y && conda activate st

mamba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
pip install patchify
pip install fastremap
pip install roifile

pip install spotsweeper

cd / && mkdir -p repo && cd repo

git clone https://github.com/JiazhangCai/SpaDiff.git
mamba install tensorflow numpy matplotlib pandas imageio multiprocess -y
mamba install pathos tqdm scipy ipykernel -y

git clone https://github.com/illuminate6060/SpotGF.git
cd SpotGF
pip install .

git clone https://github.com/danielgchen/SpaceBender.git










# R
conda install bioconda::bioconductor-spotclean -y