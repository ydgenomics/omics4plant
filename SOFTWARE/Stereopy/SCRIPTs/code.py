import stereo as st
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import cv2
import sys
import os

mask_path='/data/work/cell_cut_out/Y00710F6.raw.mRNA_mask.tif'
raw_gef='/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'
st.io.read_gef_info(raw_gef)
data = st.io.read_gef(file_path=raw_gef, bin_type='bins', bin_size=1)

mask = cv2.imread(mask_path, cv2.IMREAD_GRAYSCALE)

# 获取所有点的坐标 (通常是 N x 2 的数组)
coords = data.position.astype(int)

# 提取 x 和 y
x = coords[:, 0]
y = coords[:, 1]

keep_mask = mask[y, x] > 128

# 只保留 mask 内部的点
data.sub_by_index(cell_index=keep_mask)

# 打印过滤后的结果确认
print(f"过滤后的细胞数: {data.cell_names.shape[0]}")
print(f"过滤后的基因数: {data.gene_names.shape[0]}")

st.io.write_mid_gef(data=data,output='/data/work/raw.filtered.gef')