# [SpotGF](https://github.com/illuminate6060/SpotGF)


空间分辨转录组学（SRT）技术结合了高通量基因测序和组织学技术，提供了具有空间背景的基因表达数据，并通过细胞壁钙黄素白（CW）或细胞核4,6-二氨基-2-苯基吲哚二盐酸盐（DAPI）染色等图像进行补充。SRT数据包含了空间位置信息，有助于分析细胞类型，促进细胞功能发现、细胞相互作用研究及其他分析。理想情况下，每个位于特定位置的测序单元（在不同技术中称为珠子或点）应仅捕获原位细胞释放的转录本。然而，由于在液体实验环境中随机扩散，SRT也可能捕获异位转录本。这种扩散在SRT数据中引入了复杂的噪声，超出了单细胞RNA测序（scRNA-seq）中常见的dropout现象。

## Usage
参数选择 https://github.com/illuminate6060/SpotGF/issues/7
如果噪声不高，-p 应该偏大；如果需要精细空间去噪，-b 应该偏小。


```shell
python /data/work/test/SpotGF.py -i /data/work/test/wt10.bin1.gem -o ./wt10 -b 50 -p 0.5 -s 5

cd /data/work/spotgf
python /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/test/SpotGF.py \
-i /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/0526/bins/bin50/spotgf/gem/plp14.gem -o ./plp15 -b 10 -p 0.8 -s 5


python /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/test/SpotGF.py \
-i /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/0526/bins/bin50/spotgf/gem/wt10.gem -o ./wt10 -b 10 -p 0.8 -s 5

cd /data/work/spotgf
for input_path in /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/0526/bins/bin50/spotgf/gem/wt*.gem; do
    prefix=$(basename "$input_path" .gem)
    echo "$prefix"
    python /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/test/SpotGF.py \
    -i $input_path \
    -o ${prefix} -b 10 -p 0.8 -s 5
done

cd /data/work/spotgf
for input_path in /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/0526/bins/bin50/spotgf/gem/plp*.gem; do
    prefix=$(basename "$input_path" .gem)
    echo "$prefix"
    python /data/users/yangdong/yangdong_85a6a0468616448d95c32fbb272a0e62/online/test/SpotGF.py \
    -i $input_path \
    -o ${prefix} -b 10 -p 0.8 -s 5
done
```

```shell
positional arguments:

  -i                Input SRT data files path (.gem). 

  -o                Outpath for saving results (dir). 
    
  -b                Denoising resolution binsize, must int type, default=10.

  -lower            Lower limit for tissue structures capturing optimization, default=0.

  -upper            Upper limit for tissue structures capturing optimization', default=sys.float_info.max.

  -max_iterations   maximum number of iterations when capturing tissue structures', default=10000.

  -p                Proportion of gene numbers, must float type [0,1], default=0.5.

  -auto_threshold   if True, return a denoised gem file using automatic threshold, default=True.

  -v                Visualize SpotGF-denoised data, default=True.

  -s                Spot size used for spatial expression figure, default=5.

  -a                Alpha for tissue boundary detection, default use auto optimizealpha, default=0.
```

```shell
mamba create -n st python=3.13 pandas pot numpy scipy matplotlib descartes scanpy seaborn alphashape shapely -y
git clone https://github.com/illuminate6060/SpotGF.git
cd SpotGF
pip install .
```

- 2024|Cell systems SpotGF *使用基于最优传输的基因过滤算法去噪空间分辨转录组数据* https://mp.weixin.qq.com/s/0jDVNFJht5XTnsVDVg43YQ https://github.com/illuminate6060/SpotGF