getGeneAnno.old <- function(annoDb='hgu95av2.db', geneID, type){
  # example useage: getGeneAnno(geneID='ENSG00000181555', type='ensembl')
  kk <- unlist(geneID)
  require(annoDb, character.only = TRUE)
  annoDb <- eval(parse(text=annoDb)) 
  if (type == "Entrez Gene ID") {
    kt <- "ENTREZID"
  } else if (grepl('Ensem',type, ignore.case=TRUE)) {
    kt <- "ENSEMBL"
  } else {
    message("geneID type is not supported...\tPlease report it to developer...\n")
    return(NA)
  }

  i <- which(!is.na(kk))
  kk <- gsub("\\.\\d+$", "", kk)

  ann <- tryCatch(
    suppressWarnings(AnnotationDbi::select(annoDb,
                                           keys=unique(kk[i]),
                                           keytype=kt,
                                           columns=c("ENSEMBL", "SYMBOL"))),
    error = function(e) NULL)

  if (is.null(ann)) {
    warning("ID type not matched, gene annotation will not be added...")
    return(NA)
  }


  return(ann)
}

getGeneAnno <- function(geneID,annoDb='EnsDb.Hsapiens.v75'){
    type='ensembl'
    # example useage: getGeneAnno(geneID='ENSG00000181555', type='ensembl')
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))

    if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
    } else if (grepl('Ensem',type, ignore.case=TRUE)) {
        #kt <- "ENSEMBL"
        kt <- "GENEID"
    } else {
        message("geneID type is not supported...\tPlease report it to developer...\n")
        return(NA)
    }

    i <- which(!is.na(kk))
    kk <- gsub("\\.\\d+$", "", kk)

    ann <- tryCatch(
        suppressWarnings(AnnotationDbi::select(annoDb,
                                               keys=unique(kk[i]),
                                               keytype=kt,
                                               columns=c("GENEID", "SYMBOL"))),
        error = function(e) NULL)

    if (is.null(ann)) {
        warning("ID type not matched, gene annotation will not be added...")
        return(NA)
    }
    colnames(ann) = c('gene_id', 'gene_name')
    return(ann)
}


getGeneAnnoEn <- function(geneID, type='ensembl',annoDb='EnsDb.Hsapiens.v75'){
    # example useage: getGeneAnno(geneID='ENSG00000181555', type='ensembl')
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text=annoDb))

    if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
    } else if (grepl('Ensem',type, ignore.case=TRUE)) {
        #kt <- "ENSEMBL"
        kt <- "GENEID"
    } else if (type=='tx') {
      kt <- 'TXID'
    } else {
        message("geneID type is not supported...\tPlease report it to developer...\n")
        return(NA)
    }

    i <- which(!is.na(kk))
    kk <- gsub("\\.\\d+$", "", kk)

    output_col = c(kt, "GENEID", "SYMBOL", "GENEBIOTYPE", "ENTREZID")
    output_col = output_col[!duplicated(output_col)]
    ann <- tryCatch(
        suppressWarnings(AnnotationDbi::select(annoDb,
                                               keys=unique(kk[i]),
                                               keytype=kt,
                                               columns=output_col)),
        error = function(e) NULL)

    if (is.null(ann)) {
        warning("ID type not matched, gene annotation will not be added...")
        return(NA)
    }
    if (type=='tx'){
      colnames(ann) = c('tx_id','gene_id', 'gene_name', 'biotype', 'entrez_id')
    } else{
      colnames(ann) = c('gene_id', 'gene_name', 'biotype', 'entrez_id')
    }

    return(ann)
}
