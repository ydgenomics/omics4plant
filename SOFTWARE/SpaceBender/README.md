# SpaceBender

SpaceBender 是将 CellBender 的深度生成能力成功搬到了“空间维度”。对于处理复杂的植物、动物或肿瘤组织切片中普遍存在的 RNA 漂移和扩散污染，它提供了一个极具鲁棒性的前沿解决方案。**SpaceBender（生成模型派）**： 它的目的是“保点、保基因、抠噪声”。通过 VAE 深度学习去估算每个 Spot、每个基因里有多少分子是漂移过来的，然后把这部分分子数“扣掉”，返回的是一个一模一样维度但数值更干净的 Count 矩阵。**SpotGF（信息过滤派）**： 它的核心是利用最优输运（OT）计算每个基因的空间分布是“弥散的”还是“聚集的”。它给基因打分，核心操作是通过过滤掉那些没有空间结构的“弥散噪声基因”来达到去噪目的。

```shell
# install packages
!pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu130
!pip install pyro-ppl scanpy tables seaborn matplotlib
!git clone https://github.com/danielgchen/SpaceBender.git
```

```python
# import SpaceBender package, attach to path via "sys" package
import sys
sys.path.insert(0, '/path/to/<PACKAGE>')
from <PACKAGE> import spacebender
# import local SpaceBender package
from SpaceBender import spacebender
```

- 2026|biorxiv SpaceBender: *Denoising Spatial Transcriptomics Data to Enhance Biological Signals* https://github.com/danielgchen/SpaceBender