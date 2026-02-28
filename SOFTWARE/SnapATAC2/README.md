

snapATAC2似乎只用fragments(chromosom/start/end/barcode/count)就可以构建出来
```shell
# zless -S GSM5276944_3533EL_ATAC_fragments.tsv.gz
chr1    10157   10184   GCCAGACGTGTCCCAG-1      1
chr1    10345   10519   TGCTATTCAGGTCCTG-1      1
chr1    10594   11036   AGCTATGAGTGTCCCG-1      1
```

```python
# custom reference genome and annotation support
import snapatac2 as snap
myCustomGenome = snap.genome.Genome(fasta='/path/to/your/fasta.fasta', annotation='path/to/your/annotation.gtf')
```

GLUE
scJoint https://github.com/SydneyBioX/scJoint

```shell
pip install anndata==0.9.2
pip install scFountain
pip install scib
```

## Ref
- snapATAC2 分析教程 https://mp.weixin.qq.com/s/uNbjxMCPiv9n5BtcIzSM4g
- scATAC整合
  > - snapATAC2 多样本整合及去批次效应 https://mp.weixin.qq.com/s/rlIumT7xP0Yy0gCx6ERJdg [citation](https://pubmed.ncbi.nlm.nih.gov/40236024/)
  > - Nat Mach Intell | 南开陈盛泉团队提出scATAC-seq数据批次效应校正方法**Fountain** [github](https://github.com/BioX-NKU/Fountain) [tutorial](https://scglue.readthedocs.io/en/latest/)
- 跨模态RNA和ATAC整合
  > - QB期刊 | 单细胞多组学“双剑”合璧，scRNA-seq + scATAC-seq整合方法大盘点 https://mp.weixin.qq.com/s/CpGLdBnI_c58Rz_qqoLSFg
  > - Nat Methods | scMultiBench: 单细胞多组学整合方法的多任务基准测试 https://mp.weixin.qq.com/s/VoEASFyG2CyLxBIyNFTkSw
  > - 未配对和配对单细胞RNA-seq和ATAC-seq数据联合整合的基准算法 https://mp.weixin.qq.com/s/gMZnxy6pNi1dSfBkO3tAkA