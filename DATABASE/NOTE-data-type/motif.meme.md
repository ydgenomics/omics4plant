# motif.meme
这是 **MEME Suite**（Motif-based sequence analysis tools）使用的标准motif文件格式，广泛应用于转录因子结合位点分析。

---

## 一、格式概述

**motif.meme** 是 MEME Suite 的 **Minimal MEME Motif Format** 文件，用于存储序列motif（如转录因子结合位点）的位置特异性概率矩阵（Position-Specific Probability Matrix, PSPM）。

该格式是**纯文本格式**，可被MEME Suite中的多个工具直接读取，包括MAST、Tomtom、FIMO、GOMO等。

---

## 二、文件结构详解

一个标准的 `.meme` 文件包含以下**5个主要部分**：

### 1. 版本声明（必需）
```
MEME version 4
```
- 必须位于文件开头
- 用于确认文件格式版本，确保兼容性

### 2. 字母表声明（推荐）
```
ALPHABET= ACGT          # DNA序列
ALPHABET= ACDEFGHIKLMNPQRSTVWY   # 蛋白质序列
```
- 定义motif的序列类型
- 支持DNA（ACGT）、RNA（ACGU）和蛋白质序列

### 3. 链信息（可选，仅DNA）
```
strands: + -            # 双链（正负链）
strands: +              # 仅正链
```
- 指示motif是否来源于双链DNA序列

### 4. 背景频率（推荐）
```
Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25
```
- 各碱基/氨基酸在源序列中的背景频率
- 用于概率矩阵与对数几率矩阵的相互转换
- 若未提供，默认使用均匀分布

### 5. Motif数据（必需）
每个motif包含以下子部分：

#### **Motif名称行**
```
MOTIF MA0002.1 RUNX1
```
- `MA0002.1`：motif标识符（需唯一）
- `RUNX1`：可选的替代名称

#### **字母概率矩阵**（推荐，大多数工具使用）
```
letter-probability matrix: alength= 4 w= 18 nsites= 18 E= 1.1e-006
0.611111 0.000000 0.055556 0.333333
0.555556 0.000000 0.111111 0.333333
...
```
- `alength`：字母表长度（DNA=4，蛋白质=20）
- `w`：motif长度（位置数）
- `nsites`：用于构建motif的位点数量（默认20）
- `E`：motif的E-value（默认0）
- **每行代表一个位置**，列按字母表顺序排列（A,C,G,T）

#### **对数几率矩阵**（可选，主要用于MAST）
```
log-odds matrix: alength= 4 w= 18 E= 1.1e-006
101  -1081  -182    13
87   -1081   -82    13
...
```
- 计算公式：`score = round((log(p/b) / log(2)) * 100)`
- `p`：motif中该位置的概率，`b`：背景频率

#### **URL行**（可选）
```
URL http://example.com/motif_info
```
- 指向该motif的详细信息网页

---

## 三、完整示例文件

### DNA Motif 示例
```
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF RUNX1 MA0002.1
letter-probability matrix: alength= 4 w= 6 nsites= 100 E= 1.5e-10
0.10 0.80 0.05 0.05
0.90 0.02 0.03 0.05
0.05 0.10 0.80 0.05
0.05 0.05 0.10 0.80
0.70 0.10 0.10 0.10
0.20 0.30 0.30 0.20
URL http://jaspar.genereg.net/matrix/MA0002.1
```

### 蛋白质 Motif 示例
```
MEME version 4

ALPHABET= ACDEFGHIKLMNPQRSTVWY

Background letter frequencies
A 0.071 C 0.029 D 0.069 E 0.077 F 0.043 G 0.057 H 0.026 I 0.048
K 0.085 L 0.087 M 0.018 N 0.053 P 0.032 Q 0.029 R 0.031 S 0.058
T 0.048 V 0.069 W 0.017 Y 0.050

MOTIF ELM_LIG_3 3
letter-probability matrix: alength= 20 w= 5 nsites= 50 E= 0.001
0.10 0.05 0.60 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05
...
```

---

## 四、关键特性与注意事项

| 特性 | 说明 |
|------|------|
| **无间隙** | 单个motif不包含间隙（gaps），可变长度间隙会被拆分为多个motif |
| **概率约束** | 概率矩阵每行必须和为1 |
| **格式兼容** | 支持MEME、STREME、DREME等工具的输出直接作为输入 |
| **格式转换** | MEME Suite提供从TRANSFAC、JASPAR等格式转换的脚本 |

---

## 五、应用场景

1. **转录因子结合位点预测**：使用FIMO扫描基因组序列
2. **Motif比较**：使用Tomtom与已知数据库（如JASPAR、HOCOMOCO）比对
3. **序列分析**：使用MAST进行motif搜索和评分
4. **富集分析**：使用AME或CentriMo进行motif富集分析

如需进一步了解如何解析或生成此类文件，可以参考Bioconductor的`universalmotif`包或Python的`Bio.motifs`模块。