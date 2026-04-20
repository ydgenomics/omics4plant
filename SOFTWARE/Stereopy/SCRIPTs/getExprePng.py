# 260419

import stereo as st
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import cv2
import sys
import os


raw_gef='/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'


st.io.read_gef_info(raw_gef)
data = st.io.read_gef(file_path=raw_gef, bin_type='bins', bin_size=1)
data.tl.cal_qc(); print(data)


# 获取坐标的最大值作为矩阵维度
max_x = data.cells_matrix['spatial'][:, 0].max() + 1
max_y = data.cells_matrix['spatial'][:, 1].max() + 1

# 构建 COO 格式的稀疏矩阵
# 注意：coo_matrix((data, (row, col)), shape=(M, N))
spatial_mtx = coo_matrix(
    (data.cells['total_counts'].values, 
    (data.cells_matrix['spatial'][:, 0], data.cells_matrix['spatial'][:, 1])),
    shape=(max_x, max_y)
); print(spatial_mtx.shape)
spatial_mtx *= 25
mtx = np.clip(spatial_mtx.toarray(), 0, 255)
cv2.imwrite('raw.png', mtx)