# editor: yangdong
# 260419
# input: gef_path

import stereo as st
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import cv2
import sys
import os


gef_path='/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'
prefix = os.path.splitext(os.path.basename(gef_path))[0]  # Returns "Y00710F6.raw"


st.io.read_gef_info(gef_path)
data = st.io.read_gef(file_path=raw_gef, bin_type='bins', bin_size=1)
data.tl.cal_qc(); print(data)


max_x = data.cells_matrix['spatial'][:, 0].max() + 1
max_y = data.cells_matrix['spatial'][:, 1].max() + 1

# coo_matrix((data, (row, col)), shape=(M, N))
spatial_mtx = coo_matrix(
    (data.cells['total_counts'].values, 
    (data.cells_matrix['spatial'][:, 0], data.cells_matrix['spatial'][:, 1])),
    shape=(max_x, max_y)
); print(spatial_mtx.shape)
spatial_mtx *= 25
mtx = np.clip(spatial_mtx.toarray(), 0, 255)
cv2.imwrite(prefix + '.png', mtx)