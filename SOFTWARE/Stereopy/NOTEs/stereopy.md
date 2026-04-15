
查看`raw.gef`的大小`'width': 23519, 'height': 23519`
```python
import stereo as st
import warnings
warnings.filterwarnings('ignore')
data_path = '/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'
data_info = st.io.read_gef_info(data_path, bin_size=1)
data_info
# {'bin_list': ['bin1'], 'resolution': 500, 'gene_count': 138690, 'offsetX': 0, 'offsetY': 0, 'width': 23519, 'height': 23519, 'maxExp': 99}
```

查看tif像素
```shell
# 查看tif像素
python3 -c "from PIL import Image; Image.MAX_IMAGE_PIXELS = None; img = Image.open('/data/work/Y00710F6/Y00710F6.mRNA.png'); print(f'宽度(Width): {img.width}, 高度(Height): {img.height}')"
# 宽度(Width): 23519, 高度(Height): 23519 # tissue.gef
# 宽度(Width): 23520, 高度(Height): 23520 # bin1_img_tissue_cut.tif
# 宽度(Width): 20021, 高度(Height): 19651 # gem.png
# 26474 27631 .ss
# 
```

查看gem文件内容
```shell
$ head /data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.tissue.gem
#FileFormat=GEMv0.2
#SortedBy=None
#BinType=Bin
#BinSize=1
#Omics=Transcriptomics
#Stereo-seqChip=Y00710F6
#OffsetX=0
#OffsetY=0
geneID  geneName        x       y       MIDCount        ExonCount
CP000001        CP000001        6565    2515    1       1
```
