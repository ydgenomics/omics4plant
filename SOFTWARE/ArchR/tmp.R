getMetadata <- function(proj){
    print(proj)
    print(head(proj@cellColData))
    df <- proj@cellColData
    rownames(df) <- gsub('#', '_', rownames(df))
    write.csv(df, 'archr_metadata.csv')
    # df <- read.csv('archr_metadata.csv')
}