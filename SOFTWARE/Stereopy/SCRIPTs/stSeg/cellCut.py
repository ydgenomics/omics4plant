from stereo.tools.cell_cut import CellCut
import os
os.getcwd()

bgef_path = "/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef" #（使用tissue.gef会报错）
mask_path = "/data/work/cell_cut_out/Y00710F6.raw.mRNA_mask.tif" #（上一步的输出文件）

cc = CellCut(cgef_out_dir='.')
out_path = cc.cell_cut(bgef_path=bgef_path, mask_path=mask_path)
# 输出文件路径为当前路径下的Y00710F6.raw.mRNA_cp.gef