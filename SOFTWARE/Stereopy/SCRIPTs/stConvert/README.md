[gefpy](https://gefpy.readthedocs.io/en/main/Examples.html)

- gef2gem
- gem2gef
- gef2anndata
- anndata2gef
- anndata2rds
- rds2anndata


```python
def gef2gem(gef_path, gem_path=None, binsize=None):
    '''
    将GEF文件转换为GEM文件，支持不同的输入参数组合
    设计：导入需要的库；兼容不同的参数输入方式；最后使用gefpy将gef文件转换为gem文件
    ref: https://gefpy.readthedocs.io/en/main/Examples.html#generate-gem-by-gef
    image: cellpose_035
    '''
    from gefpy.gef_to_gem_cy import gefToGem
    import h5py
    import os
    with h5py.File(gef_path, 'r') as f:
        if 'geneExp' in f:
            print("支持的分辨率（Bin sizes）:", list(f['geneExp'].keys()))
            key_list = list(f['geneExp'].keys())
        else:
            print("文件结构可能为 Cell-bin，包含的顶级组为:", list(f.keys()))
            key_list = list(f.keys())
    if gem_path is None:
        prefix = os.path.splitext(os.path.basename(gef_path))[0]
        if binsize is not None:
            prefix += f'.bin{binsize}'
        else:
            print("未指定binsize，默认使用文件中第一个支持的分辨率。")
            binsize = int(key_list[0].replace('bin', ''))
            prefix += f'.bin{binsize}'
        gem_path = prefix + '.gem'
    else:
        print(f'[gef2gem] gem_path provided: {gem_path}')
        prefix = os.path.splitext(os.path.basename(gem_path))[0]
        if binsize is not None:
            print(f'[gef2gem] binsize provided: {binsize}')
        else:
            print("未指定binsize，默认使用文件中第一个支持的分辨率。")
            binsize = int(key_list[0].replace('bin', ''))
            print(f'[gef2gem] 使用默认binsize: {binsize}')
    print(f'[gef2gem] prefix: {prefix}; gem_path: {gem_path}')
    obj = gefToGem(gem_path, prefix)
    obj.bgef2gem(gef_path, int(binsize))
    

gef2gem(gef_path="/data/work/0511/tissue_seg/wt10.tif_flt.gef")
gef2gem("/data/work/0511/tissue_seg/wt10.tif_flt.gef", "wt10.bin1.gem", 1)
gef2gem("/data/work/0511/tissue_seg/plp15.png_flt.gef", "plp15.bin1.gem", 1)


def gem2h5ad(gem_path, h5ad_path=None, bin_type=None):
    '''
    将GEM文件转换为h5ad文件，支持不同的输入参数组合
    ref: https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_gem.html#scanpy.read_gem
    image: cellpose_035
    '''
    import stereo as st
    import warnings
    warnings.filterwarnings('ignore')
    import os
    if bin_type is None:
        bin_type = 'bins'
    prefix = os.path.splitext(os.path.basename(gem_path))[0]
    data = st.io.read_gem(file_path=gem_path, sep='\t', bin_type=bin_type, is_sparse=True)
    st.io.stereo_to_anndata(data, flavor='seurat',output=prefix + '_seurat.h5ad')
    st.io.stereo_to_anndata(data, flavor='scanpy',output=prefix + '_scanpy.h5ad')
    
gem2h5ad("/data/work/test/wt10.bin1.gem")

# read the gem file
data_path = './SS200000135TL_D1.cellbin.gem'

```


```R
# image: SuperCell
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(ggplot2)
    library(BayesSpace)
    library(Seurat)
    library(cowplot)
    library(RColorBrewer)
    library(stringr)
    library(Seurat)
    library(SeuratDisk)
})
sce <- as.SingleCellExperiment(seurat_spatialObj, assay="SCT")
rds2h5ad <- function(sce, prefix){
    sce.seurat <- as.Seurat(sce, counts="counts", data="logcounts")
    locations <- seurat_spatialObj@meta.data[,c('x','y')]
    sce.seurat[['image']] <- new(Class = 'SlideSeq', assay = "Spatial", coordinates = locations)
    SaveH5Seurat(sce.seurat, filename = paste0(prefix, ".h5Seurat"), assay="Spatial", overwrite=TRUE)
    Convert(paste0(prefix, ".h5Seurat"), dest = "h5ad", overwrite=TRUE)
}
```