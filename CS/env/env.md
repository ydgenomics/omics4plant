```shell
source /opt/software/miniconda3/bin/activate
conda config --remove-key channels
conda config --add channels conda-forge

conda create -n r4.4 r-base=4.4 -y
conda activate r4.4
conda install conda-forge::r-optparse -y
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes -y
```

```shell
source /opt/software/miniconda3/bin/activate
conda create -n seurat r-base=4.4 -y
conda activate seurat

conda install conda-forge::r-optparse -y
conda install conda-forge::r-biocmanager -y
conda install conda-forge::r-devtools -y
conda install conda-forge::r-remotes -y
# 1_qc
conda install conda-forge::r-seurat -y
conda install conda-forge::r-soupx -y
conda install bioconda::bioconductor-dropletutils -y
conda install bioconda::bioconductor-scdblfinder -y
conda install bioconda::bioconductor-singler -y
conda install bioconda::bioconductor-scater -y
conda install bioconda::r-liger -y


# 02_anno
conda install bioconda::bioconductor-singler -y
conda install conda-forge::r-hgnchelper -y
conda install conda-forge::r-tidyverse -y
conda install conda-forge::r-ggraph -y
conda install conda-forge::r-data.tree -y
conda install bioconda::bioconductor-scater -y

# cluster
conda install conda-forge::r-clustree -y

# 03_integrate
# Rscript -e 'install.packages("bbknnR")'
conda install pwwang::r-seuratdata -y
conda install bioconda::r-liger -y
conda install conda-forge::r-rcppplanc -y
conda install pwwang::r-seuratwrappers -y
conda install bioconda::bioconductor-glmgampoi -y
conda install conda-forge::r-harmony -y

# 04_metaneighbor
conda install bioconda::bioconductor-metaneighbor -y
conda install bioconda::bioconductor-complexheatmap -y
conda install conda-forge::r-circlize -y

# Rscript -e "install.packages('presto')"

# 06_enrich
conda install bioconda::bioconductor-clusterprofiler -y
conda install bioconda::bioconductor-keggrest -y
conda install bioconda::bioconductor-annotationforge -y


source /opt/software/miniconda3/bin/activate
conda create -n scanpy python=3.11 -y
conda activate scanpy

conda install conda-forge::scanpy -y
conda install bioconda::scrublet -y
conda install conda-forge::leidenalg -y
conda install bioconda::harmonypy -y
conda install conda-forge::scvi-tools -y
pip install --use-pep517 bbknn #pip install bbknn
conda install conda-forge::scikit-misc -y
conda install conda-forge::plotly -y


conda install -c conda-forge gcc=12 gxx=12 -y
conda install conda-forge::h5py -y
pip install scib-metrics # Python >=3.11
conda install conda-forge::scanpy -y
conda install conda-forge::ipykernel -y

# 05_dea
pip install memento-de # memento
pip install adjustText

# 07_pseudotime
conda install -c conda-forge cellrank -y
conda install conda-forge::moscot -y
# conda install -c conda-forge -c bioconda palantir -y #python要求3.12
# pip install --user magic-impute
conda install conda-forge::certifi -y
conda install conda-forge::rpy2 -y
# pip install fa2-modified
conda install conda-forge::fa2 -y #不支持安装在大于等于3.12python
conda install anaconda::click -y
# python -m ipykernel install --user --name cellrank2 --display-name "Python (cellrank2)"





# 02_cluster
# conda install bioconda::bioconductor-bluster -y
Rscript -e 'BiocManager::install("bluster")'
# Rscript -e 'devtools::install_github("pengminshi/mrtree")'
# Rscript -e 'devtools::install_github("imbs-hl/ranger")'
# Rscript -e 'remotes::install_github("corceslab/CHOIR", ref="main", repos = BiocManager::repositories(), upgrade = "never")' # CHOIR
source /opt/software/miniconda3/bin/activate
conda activate scline
wget https://github.com/pengminshi/MRtree/archive/refs/heads/master.zip
Rscript -e 'devtools::install_local("master.zip", dependencies = FALSE)'
rm master.zip
Rscript -e 'install.packages("ranger")'
wget https://github.com/corceslab/CHOIR/archive/refs/heads/main.zip
Rscript -e 'remotes::install_local("main.zip", dependencies = FALSE)'
rm main.zip





# conda install bioconda::scib -y # python版本冲突
# pip install scib-metrics






# # 06_enrich
# conda install bioconda::bioconductor-clusterprofiler -y
# conda install bioconda::bioconductor-keggrest -y
# conda install bioconda::bioconductor-annotationforge -y
# Rscript -e 'install.packages("openxlsx", repos="https://cloud.r-project.org")' # install failed with conda

# # recall
# conda install bioconda::bioconductor-scdesign3 -y
# conda install r::r-lamw -y
# Rscript -e "install.packages('knockoff')"
# Rscript -e 'install.packages("countsplit")'
```

# scib
```shell
# scib-metrics
conda create -n scib python=3.11 -y
conda activate scib
conda install -c conda-forge gcc=12 gxx=12 -y
conda install conda-forge::h5py -y
pip install scib-metrics # Python >=3.11
conda install conda-forge::scanpy -y
pip install --use-pep517 bbknn #pip install bbknn
conda install conda-forge::ipykernel -y

# 07_pseudotime
source /opt/software/miniconda3/bin/activate
conda activate scib
conda install -c conda-forge cellrank -y
conda install conda-forge::moscot -y
# conda install -c conda-forge -c bioconda palantir -y #python要求3.12
# pip install --user magic-impute
conda install conda-forge::certifi -y
conda install conda-forge::rpy2 -y
# pip install fa2-modified
conda install conda-forge::fa2 -y #不支持安装在大于等于3.12python
conda install anaconda::click -y
# python -m ipykernel install --user --name cellrank2 --display-name "Python (cellrank2)"

# genes2genes
pip install git+https://github.com/Teichlab/Genes2Genes.git -i https://mirrors.aliyun.com/pypi/simple/ 
pip install optbinning # https://gnpalencia.org/optbinning/installation.html
```

# convert [安装sceasy](https://mp.weixin.qq.com/s/MzMHKjHl-V_XWKYw9DEmCg)
```shell
source /opt/software/miniconda3/bin/activate
conda install conda-forge::mamba -y
mamba create -n convert r-seurat=4.4 r-devtools python=3.11 -y
mamba activate convert
mamba install -c conda-forge scanpy -y
mamba install bioconda::r-sceasy -y
# Rscript -e 'devtools::install_github("cellgeni/schard")' # fail
wget https://github.com/cellgeni/schard/archive/refs/heads/main.zip
Rscript -e 'devtools::install_local("main.zip")'
rm main.zip
conda install conda-forge::r-reticulate -y
conda install conda-forge::loompy -y
conda install conda-forge::r-optparse -y
```

# genes2genes
```shell
source /opt/software/miniconda3/bin/activate
conda create --name g2g python=3.8 -y
conda activate g2g
# pip install genes2genes
# pip install git+https://github.com/Teichlab/Genes2Genes.git
pip install git+https://github.com/Teichlab/Genes2Genes.git -i https://mirrors.aliyun.com/pypi/simple/ 
pip install optbinning # https://gnpalencia.org/optbinning/installation.html
conda install conda-forge::ipykernel -y
conda install conda-forge::scanpy -y
# python -m ipykernel install --user --name g2g --display-name "Python (g2g)"
```