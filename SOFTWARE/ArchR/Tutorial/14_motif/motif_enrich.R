# 260405

library(ArchR)
set.seed(1)
library(BSgenome.rice.test)

archr_project='/data/work/rice/ArchR/work/Save-ZHH-0d'
pwm_list_rdata="/data/work/rice/ref/motif/pwm_list.rdata"


proj <- loadArchRProject(archr_project)
load(pwm_list_rdata)
proj <- addMotifAnnotations(ArchRProj = proj, motifPWMs=pwm_list, name = "Motif")