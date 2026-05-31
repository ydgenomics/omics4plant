source /opt/software/miniconda3/bin/activate
mamba create -n st python=3.8 -y && conda activate st
mamba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
pip install patchify
pip install fastremap
pip install roifile
pip install ipython
pip install ipykernel
python -m ipykernel install --user --name st --display-name "Python (st)"

pip install descartes
pip install alphashape

pip install spotsweeper