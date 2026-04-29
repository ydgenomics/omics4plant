# editor: yangdong
# image: ArchR_Macs2_ChromVARmotifs
# 260427

library(ArchR)
# library(BSgenome.rice.test) # required
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4){stop('
### Usage:
Rscript footprinting.R <archr_project> <atac_key> <threads> <bsgenome_path>
### Example:
archr_project="EFH-0d"
atac_key="Clusters"
threads=8
bsgenome_path="/data/work/archr/region_annotation/BSgenome.species_1.0.0.tar.gz"
Rscript ../footprinting.R \
$archr_project $atac_key $threads $bsgenome_path
')}

archr_project <- args[1]
atac_key <- args[2]
threads <- as.integer(args[3])
bsgenome_path <- args[4]

system(paste0("R CMD INSTALL ", bsgenome_path))
bsgenome_name <- sub("_1.0.0.tar.gz$", "", basename(bsgenome_path))
do.call("library", list(bsgenome_name))

# markerMotifs=NULL
prefix <- basename(archr_project)
addArchRThreads(threads = threads)


proj <- loadArchRProject(archr_project)
proj

## motif footprinting
motifPositions <- getPositions(proj, name = "Motif")
motifPositions

head(names(motifPositions))
# [1] "LOC-Os01g03720" "LOC-Os01g07120" "LOC-Os01g07480" "LOC-Os01g08160"
# [5] "LOC-Os01g09550" "LOC-Os01g09640"

if(is.null(proj@projectMetadata$GroupCoverages[[atac_key]])){
  proj <- addGroupCoverages(ArchRProj = proj, groupBy = atac_key, force = TRUE)
}

# ?markerMotifs(gene names): only footprint interested motifs, e.g. TFs of interest
if (!exists("markerMotifs") || is.null(markerMotifs) || length(markerMotifs) == 0) {
  markerMotifs <- names(motifPositions)
}
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = atac_key
)

## different methos to normalize the Tn5 bias, e.g. "Subtract", "Divide", "None"
# --- Subtracting the Tn5 Bias ---
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = paste0(prefix, "_Footprints-Subtract-Bias"),
  addDOC = FALSE,
  smoothWindow = 5
)

# --- Dividing by the Tn5 Bias ---
plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "Divide",
    plotName = paste0(prefix, "_Footprints-Divide-Bias"),
    addDOC = FALSE,
    smoothWindow = 5
)

# --- No Normalization ---
plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "None",
    plotName = paste0(prefix, "_Footprints-No-Normalization"),
    addDOC = FALSE,
    smoothWindow = 5
)

## ? Feature Footprinting
# 可对感兴趣的区段看是否有凹下去的情况，例如新发现的调控binding
# seTSS <- getFootprints(
#   ArchRProj = projHeme5, 
#   positions = GRangesList(TSS = getTSS(projHeme5)), 
#   groupBy = "Clusters2",
#   flank = 2000
# )

saveArchRProject(proj)