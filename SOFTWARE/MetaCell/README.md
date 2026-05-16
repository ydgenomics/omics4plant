# metacells
> —defined as disjoint groups of highly similar cells which are aggregated together—have been proposed as a way to simultaneously reduce the size and improve the signal-to-noise ratio in single-cell genomics data

main aspects of the usage of metacells: 
- (i) enhancing signal in sparse scRNA-seq data 
- and (ii) lowering computational burden due to the large size of single-cell genomics data.

Process: The aggregation is usually done by either summing raw counts (Baran et al, 2019; Iacono et al, 2018; Persad et al, 2023) or by averaging normalized counts.

Overall, these observations suggest that the choice of graining levels depends on both the complexity and size of the data. For large and low-complexity data, a relatively high graining level may be used. For datasets with higher complexity or lower size, it is recommended to use lower graining levels to preserve the underlying heterogeneity and ensure that distinct cell populations remain distinguishable.