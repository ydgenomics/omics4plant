根据前面的对话，我帮你梳理了 OrthoFinder 结果的核心解读框架，可以当作速查手册。

---

## 一、先看核心统计（日志里最关键的两句话）

你运行完 OrthoFinder 后，日志里会印出类似这样的信息：

> OrthoFinder assigned **73691 genes (78.9% of total)** to **17077 orthogroups**. Fifty percent of genes were in orthogroups with **4 or more genes (G50 was 4)** and were contained in the largest **6405 orthogroups (O50 was 6405)**. There were **10923 orthogroups with all species present** and **2956 of these consisted entirely of single-copy genes**.

| 指标 | 含义 | 你的数据 | 评价 |
|------|------|----------|------|
| **基因分配率** | 多少基因被划入了某个直系同源群 | 78.9% | 不错，说明物种和注释质量较好 |
| **G50** | 把基因按所在群的大小降序排列，累计到一半基因时，该群至少包含的基因数 | 4 | 偏小，说明有大量小群（可能是远缘物种） |
| **O50** | 承载上述一半基因所需的最大群数量 | 6405 | 印证了小群很多 |
| **所有物种共有群** | 每个物种都有至少一个基因的群 | 10923 | 数量丰富 |
| **完全单拷贝群** | 所有物种都恰好只有一个基因 | 2956 | **这是最优质的进化分析素材** |

---

## 二、结果目录按用途速查

结果路径：`/data/work/orthofinder/OrthoFinder/Results_Jul15/`

| 你想要的... | 去哪里找 |
|-------------|----------|
| **总体统计** | `Comparative_Genomics_Statistics/Statistics_Overall.tsv` |
| **每个物种的统计** | `Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv` |
| **所有直系同源群的基因列表** | `Orthogroups/Orthogroups.tsv`（一行一个群，每个物种列里是该群的基因ID，**会有多对多、一对多、多对一**） |
| **未分配基因** | `Orthogroups/Orthogroups_UnassignedGenes.tsv` |
| **完全一对一的表** | 需手动筛选（见下方命令） |
| **有根物种树** | `Species_Tree/SpeciesTree_rooted.txt` |
| **单拷贝同源群序列（做树用）** | `Single_Copy_Orthologue_Sequences/` |
| **基因树** | `Gene_Trees/` 和 `Resolved_Gene_Trees/` |
| **基因复制事件** | `Gene_Duplication_Events/` |
| **分层直系同源群** | `Phylogenetic_Hierarchical_Orthogroups/`（`N0.tsv` 是最高层级） |
| **水平基因转移候选** | `Putative_Xenologs/` |

---

## 三、`Orthogroups.tsv` 的理解（你之前疑惑的核心）

**这个表里的关系是“群”而不是“对”。**

- **直系同源群（Orthogroup）**：从一个共同祖先基因继承下来的所有后代基因的集合。
- **为什么会有多对一、一对多？**
  - 物种形成 → 产生一对一（直系同源）
  - 基因复制 → 产生一对多（旁系同源也被包含在内）
- **你的观察“拟南芥多对一景天，拟南芥一对多景天”完全正常**，反映了基因在不同支系的特异性扩张或丢失。

---

## 四、如何获得严格一对一的表

`Orthogroups_SingleCopyOrthologues.txt` 里的编号是 `N0.HOG*` 格式，不能直接用于筛选 `OG*` 格式的 `Orthogroups.tsv`。

**正确做法**：用以下命令，从总表中筛出每个物种恰好只有一个基因的行：

```bash
cd /data/work/orthofinder/OrthoFinder/Results_Jul15/Orthogroups/
awk -F'\t' 'NR==1 {print; next} {
    single=1
    for(i=2;i<=NF;i++) {
        if($i=="" || split($i,arr,",")!=1) {single=0; break}
    }
    if(single) print
}' Orthogroups.tsv > Orthogroups_One2One.tsv
```

验证：
```bash
wc -l Orthogroups_One2One.tsv   # 应该约 2957 行（含表头）
head -n 3 Orthogroups_One2One.tsv | cut -f1-3
```

---

## 五、后续分析快速指南

- **想建更好的物种树** → 用 `Single_Copy_Orthologue_Sequences/` 里的序列
- **想分析某个基因家族的扩张/收缩** → 从 `Orthogroups.tsv` 里找对应群，配合 `Gene_Duplication_Events/`
- **想研究基因水平转移** → 查看 `Putative_Xenologs/`
- **想理解古老基因家族的演化层级** → 查看 `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`

如果你有具体想做的分析（比如用那2956个单拷贝基因建物种树），我可以帮你梳理下一步的详细操作流程。