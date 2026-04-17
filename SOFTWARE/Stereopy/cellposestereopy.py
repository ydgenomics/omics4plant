from stereo.algorithm import cell_pose as cp
image="/data/work/cellpose/Y00710F6.mRNA.tif"  # (输入配准好的图像)
output="/data/work/cellpose/Y00710F6.mRNA_cp.tif"
cp.Cellpose(
    img_path=image, # (输入配准好的图像）
    out_path=output,
    model_type=None, #（按需求选择）
    dmin=30,   # min cell diameter
    dmax=40,   # max cell diameter
    step=10
)
# 模型下载网址：https://www.cellpose.org/models/cyto2torch_0.h5
# cellposestereopy3
# 镜像中已经下载好cellpose官网中能够找到的示例模型，
# 模型使用时需要将模型（多个文件）复制到模型文件路径下：以cyto2为例，标红了需要用到的指令.
# 使用外源的模型需要更改cellpose中/home/stereonote/.cellpose/models.py中的可以被识别的模型名称部分。
# 并将模型复制到模型路径/home/stereonote/.cellpose/models下。

# ls /home/stereonote/.cellpose/models && ls -all /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/cyto2torch_0 /home/stereonote/.cellpose/models
# cp /home/stereonote/model/cellpose/cyto2torch_1 /home/stereonote/.cellpose/models 
# cp /home/stereonote/model/cellpose/cyto2torch_2 /home/stereonote/.cellpose/models 
# cp /home/stereonote/model/cellpose/cyto2torch_3 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/size_cyto2torch_0.npy /home/stereonote/.cellpose/models 

# sudo cp /home/stereonote/model/cellpose/cytotorch_0 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/cytotorch_1 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/cytotorch_2 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/cytotorch_3 /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/size_cytotorch_0.npy /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/demo /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/general /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/gui_models.txt /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/LC1 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/LC2 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/LC3 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/LC4 /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/livecell /home/stereonote/.cellpose/models 
# sudo cp /home/stereonote/model/cellpose/CP /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/CPx /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/CP_20230619_095147 /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/cellpose_residual_on_style_on_concatenation_off_train1_2022_10_28_16_59_59.600570 /home/stereonote/.cellpose/models
# sudo cp /home/stereonote/model/cellpose/cellpose_residual_on_style_on_concatenation_off_train_2023_06_20_15_14_14.218607_epoch_9951 /home/stereonote/.cellpose/models

# Cellpose使用
cp.Cellpose(
    img_path=image,
    out_path=output,
    model_type='cyto2',  #根据需求挑选需要的模型，
    dmin=30,   # min cell diameter
    dmax=40,   # max cell diameter
    step=10
)

# Tif转GEF:
from stereo.tools.cell_cut import CellCut

bgef_path = "/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef" #（使用tissue.gef会报错）
mask_path = "/data/work/Y00710F6/Y00710F6.raw.mRNA_cp.tif" #（上一步的输出文件）

cc = CellCut(cgef_out_dir='.')
out_path = cc.cell_cut(bgef_path=bgef_path, mask_path=mask_path)

# 输出文件为raw.cellbin.gef 
# 这一步可以不使用GPU，根据需求调整，测试中将dmin下调后能够多出一些切片中央的细胞（需要通过测试验证合适参数），
# 上调则效果不太显著并延长细胞切割时长。