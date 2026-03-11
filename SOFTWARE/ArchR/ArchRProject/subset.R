projSubset <- subsetArchRProject(
  ArchRProj = proj,
  cells = proj$cellNames[idxSample],
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  force = TRUE
  )