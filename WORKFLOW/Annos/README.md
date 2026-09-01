## tools
- SAMap 抽样提速
- SATURN

## Db
- https://www.tobaccodb.org/pcmdb/download


> [!NOTE]
> - [P-scAnno](../../PROJECT-open/P-scAnno-20260310/README.md)
> - [WORKFLOW/Annos](../../WORKFLOW/Annos/)
> - [Annos#protocol-plant-single-cell-annotation](https://github.com/ydgenomics/Annos#protocol-plant-single-cell-annotation)


# scAnno
**scAnno**: **s**ingle **c**ell **Anno**tation pipeline. sxRNA-seq data are increasing like flooding for plants and animals, and good annotation for clustered and labeled data is basic foundament of downstream analysis. Now I supposed a rubust, avaiable, and easy piepeline to solve annotation problem of plants' sc-cells, like a diverse toolkie.

- 为什么要做好细胞注释？这既是数据的优势又是下游分析的关键。开发一个python版本的单细胞注释工具的必要性，大数据集合在python生态分析的优势。
- 项目很多，做的很多的也是细胞类型的注释，而这个也是最基础最根本的部分，要做好、做的快不容易。

## workflow
> cell types, states, and other biologically relevant patterns with the goal of creating annotated cell maps.
- 无参注释，接入大语言模型。
- 有参注释，harmony和scvi方法整合后，基于KNN标签转移。singler?
- 基于marker，自动与手动。sctype （也可接入大语言模型）
- Mapping Cell Names to the Cell Ontology/Taxonomy。植物能不能map到plantscrna，细胞相似性，metaneighbor，cellwalker2
- MetaTiME [omicverse](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_metatime.html#celltype-auto-annotation-with-metatime)
- TOSICA [omicverse](https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_tosica.html#celltype-annotation-migration-mapping-with-tosica)
- 大语言模型训练数据模型用于注释 scMulan
- Consensus annotation with CellVote


## References
- tips
  - [2026|掌握这条基础规则后，细胞注释再也不会成为单细胞分析的限速步骤](https://mp.weixin.qq.com/s/EySwQGrFWu4m5z42gqYV1w)
- [OMG Browser (Orthologous Marker Groups Browser)](https://www.omicsempower.com/blog/plant-single-cell-cell-type-annotation-omg-browser/)
- https://omicverse.readthedocs.io/en/latest/tutorials/index_single_annotation.html#annotation
- [Annos](https://github.com/ydgenomics/Annos)
- 2026 https://github.com/illuminate6060/CellBLASTer
- Online website for annotation
  - [scPlantAnnotate](https://scplantannotate.missouri.edu/)
- Markers
  - 2025|Molecular Plant | PlantscRNAdb 4.0上线：首次建立植物细胞“通用语言”，破解跨物种识别难题 [wechat](https://mp.weixin.qq.com/s/rC7RDZZpFDyaDa8BED7YNw) [PCmaster_anno](https://github.com/bioinplant/PCmaster) [HCMarker](https://github.com/daidai905/HCMarker)
- Reference & Query
  - 细胞“名片”：5大工具+全流程代码，解锁细胞注释最优解！ [wechat](https://mp.weixin.qq.com/s/fz0txK_mYAP0jxZ40I0gkw)
- iMetaOmics | 北京科技大学杜宏武组-细胞类型映射注释 [wechat]
