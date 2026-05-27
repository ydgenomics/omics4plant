# SpaceBender

SpaceBender 是将 CellBender 的深度生成能力成功搬到了“空间维度”。对于处理复杂的植物、动物或肿瘤组织切片中普遍存在的 RNA 漂移和扩散污染，它提供了一个极具鲁棒性的前沿解决方案。

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