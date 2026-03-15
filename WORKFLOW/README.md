# WORKFLOW

```mermaid
flowchart TB
1[Remove ambient] --- 1.1[Dataget-SoupX] 
1 ==> 2[Remove doublet and Merge] --- 2.1[Dataget-scrublet] 
2 --- 2.2[Dataget-scDblfinder]
2 ==> 3[scAnno]
3 --- |Marker| 3.1[scAnno-sctype] --- 3.1.1[Enrich] --- 3.1.2[AddModuleScore]
3 --- |Reference| 3.2[scAnno-singler]
3 --- 3.3[scAnno-?]
3 ==> 4[MetaNeighbor]
4 ==> 5[Integration]
5 ==> 6[DEA]
5 ==> 7[Trajectory]
5 ==> 8[GRN]
```

---
## [ydFormat](./ydFormat/)
> rds与h5ad文件的转换

## [Dataget](./Dataget/) :arrow_up: 


## [Cluster](./Cluster/) :arrow_up: 

---

## [MetaNeighbor](./MetaNeighbor/) :arrow_up: 

---

## [Annotation](./Annotation/) (optional)

---

## [Enrich](./Enrich/) (optional)

---

## [Integration](./Integration/)

---

## [DEA](./DEA/) :arrow_down: 

## [GRN](./GRN/) :arrow_down: 

## [CEA](./CEA/) :arrow_down: 

## [Trajectory](./Trajectory/) :arrow_down: 
