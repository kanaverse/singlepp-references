library(celldex)
library(matrixStats)

#########################################
#########################################

formatter <- function(x, genes, prefix, dir = "v3.0.0") {
    dir.create(dir, showWarnings=FALSE)

    mat <- assay(x, "logcounts")
    rmat <- colRanks(mat, ties.method="min")
    mhandle <- gzfile(file.path(dir, paste0(prefix, "_matrix.csv.gz")))
    write.table(file=mhandle, rmat, sep=",", col.names=FALSE, row.names=FALSE)

    for (g in colnames(genes)) {
        ghandle <- gzfile(file.path(dir, paste0(prefix, "_genes_", g, ".csv.gz")))
        column <- genes[[g]]
        if (!is.character(column)) {
            column <- paste(column, collapse="\t")
        }
        write(column, file=ghandle)
    }

    for (lab in colnames(colData(x))) {
        suff <- sub("label\\.", "", lab)
        curlabs <- colData(x)[[lab]]

        ulabs <- sort(unique(curlabs))
        uhandle <- gzfile(file.path(dir, paste0(prefix, "_label_names_", suff, ".csv.gz")))
        write(ulabs, file=uhandle, ncolumns=1)

        lhandle <- gzfile(file.path(dir, paste0(prefix, "_labels_", suff, ".csv.gz")))
        write(match(curlabs, ulabs) - 1L, file=lhandle, ncolumns=1)

        Y <- SingleR::getClassicMarkers(x, curlabs)
        curhandle <- gzfile(file.path(dir, paste0(prefix, "_markers_", suff, ".gmt.gz")), open="wb")
        for (i in names(Y)) {
            for (j in names(Y[[i]])) {
                if (i == j) {
                    next
                }
                thing <- c(match(c(i, j), ulabs), match(Y[[i]][[j]], rownames(x))) - 1L
                write(paste(thing, collapse="\t"), file=curhandle, append=TRUE)
            }
        }
        close(curhandle)
    }
}

#########################################
#########################################

library(AnnotationHub)
ahub <- AnnotationHub(cache = "dump", ask = FALSE)

species <- list(
    # ahub$ah_id[grepl("Ensembl.*EnsDb.*Mus musculus", ahub$title)]
    `10090` = c(
       "AH64944",
       "AH67971",
       "AH69210",
       "AH73905",
       "AH75036",
       "AH78811",
       "AH79718",
       "AH83247",
       "AH89211",
       "AH89457",
       "AH95775", 
       "AH98078",
       "AH100674",
       "AH104895",
       "AH109367"
    ),

    # ahub$ah_id[grepl("Ensembl.*EnsDb.*Homo sapiens", ahub$title)]
    `9606` = c(
        "AH64923",
        "AH67950",
        "AH69187",
        "AH73881",
        "AH73986",
        "AH75011",
        "AH78783", 
        "AH79689",
        "AH83216",
        "AH89180",
        "AH89426",
        "AH95744",
        "AH98047",
        "AH100643",
        "AH104864",
        "AH109336"
    )
)

#########################################
#########################################

stripNAs <- function(df) {
    drop <- Reduce("|", lapply(as.list(df), is.na))
    df[!drop,,drop=FALSE]
}

for (spec in names(species)) {
    if (spec == "9606") {
        datasets <- c("BlueprintEncode", "HumanPrimaryCellAtlas", "MonacoImmune", "DatabaseImmuneCellExpression", "NovershternHematopoietic", "MonacoImmune")
    } else {
        datasets <- c("ImmGen", "MouseRNAseq")
    }

    for (ref in datasets) {
        X <- get(paste0(ref, "Data"))()

        # Accumulating all of the gene IDs.
        if (!all(grepl("^ENS", rownames(X)))) {
            to_ensembl <- to_entrez <- list()
            for (x in species[[spec]]) {
                ens <- ahub[[x]]
                to_ensembl[[x]] <- stripNAs(select(ens, keys = rownames(X), keytype = "SYMBOL", column = "GENEID"))
                to_entrez[[x]] <- stripNAs(select(ens, keys = rownames(X), keytype = "SYMBOL", column = "ENTREZID"))
            }

            to_ensembl <- S4Vectors::unique(DataFrame(do.call(rbind, to_ensembl)))
            to_entrez <- unique(DataFrame(do.call(rbind, to_entrez)))
            genes <- DataFrame(
                symbol = rownames(X), 
                ensembl = splitAsList(to_ensembl$GENEID, factor(to_ensembl$SYMBOL, rownames(X))),
                entrez = splitAsList(to_entrez$ENTREZID, factor(to_entrez$SYMBOL, rownames(X)))
            )

        } else {
            to_symbol <- to_entrez <- list()
            for (x in species[[spec]]) {
                ens <- ahub[[x]]
                to_entrez[[x]] <- stripNAs(select(ens, keys = rownames(X), keytype = "GENID", column = "ENTREZID"))
                to_symbol[[x]] <- stripNAs(select(ens, keys = rownames(X), keytype = "GENEID", column = "SYMBOL"))
            }

            to_symbol <- unique(DataFrame(do.call(rbind, to_symbol)))
            to_entrez <- unique(DataFrame(do.call(rbind, to_entrez)))
            genes <- DataFrame(
                symbol = splitAsList(to_symbol$SYMBOL, factor(to_symbol$SYMBOL, rownames(X))),
                ensembl = rownames(X), 
                entrez = splitAsList(as.character(to_entrez$ENTREZID), factor(to_entrez$SYMBOL, rownames(X)))
            )
        }

        formatter(X, genes, ref)
    }
}
