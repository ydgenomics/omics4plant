[orthofinder](https://github.com/OrthoFinder/OrthoFinder)


[usage](https://github.com/OrthoFinder/OrthoFinder#simple-usage)


拿到一对一的基因
```shell
awk -F'\t' 'NR==1 {print; next} {single=1; for(i=2;i<=NF;i++) {if($i=="") {single=0; break} if(split($i,arr,",")!=1) {single=0; break}} if(single) print}' Orthogroups.tsv > Orthogroups_One2One.tsv
```