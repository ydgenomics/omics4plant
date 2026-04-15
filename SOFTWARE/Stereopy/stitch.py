import stereo as st

# 1. 初始化图像处理对象
# 给定小图所在的文件夹路径
img_loader = st.io.ImageLoader(server_path='path/to/ssDNA_tiles_folder/')

# 2. 自动拼接 (Stitching)
# 该函数会读取小图并根据重叠区(overlap)计算位移，生成大图
stitched_img = img_loader.stitch(method='overlap', overlap=0.1) # overlap通常是10%-20%

# 3. 同样的位移逻辑应用到 FB 通道
fb_loader = st.io.ImageLoader(server_path='path/to/FB_tiles_folder/')
# 使用ssDNA计算出的坐标偏移量来拼接FB，确保完全重合
stitched_fb = fb_loader.stitch_by_reference(reference=stitched_img)

# 4. 保存大图
stitched_img.save("ssDNA_full.tif")
stitched_fb.save("FB_full.tif")
