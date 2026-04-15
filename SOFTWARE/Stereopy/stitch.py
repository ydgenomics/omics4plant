import cv2
import numpy as np
import glob
import os

# 获取所有图片，按文件名排序（保证顺序正确）
tif_files = sorted(glob.glob("/data/work/ss/SS_*.tif"))

# 读取第一张图获取尺寸
sample_img = cv2.imread(tif_files[0], cv2.IMREAD_UNCHANGED)
h, w = sample_img.shape[:2]; print(h, w)

# 假设是 10x10 网格，需要根据实际图片数量调整
cols = 12
rows = len(tif_files) // cols

# 创建画布
canvas = np.zeros((rows * h, cols * w), dtype=sample_img.dtype)

for idx, f in enumerate(tif_files):
    img = cv2.imread(f, cv2.IMREAD_UNCHANGED)
    r = idx // cols
    c = idx % cols
    canvas[r*h:(r+1)*h, c*w:(c+1)*w] = img

cv2.imwrite("stitched_image.tif", canvas)