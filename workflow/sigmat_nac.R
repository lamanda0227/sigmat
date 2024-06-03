## Performing DE analysis using Seurat
##
## scdata: The single cell data matrix
## cell_types: A vector of the cell type annotations
##
## return: List with the differentially expressed genes for each cell type

DEAnalysis <- function(scdata, cell_types) {
    list.DE.group <- as.list(rep(0, length(cell_types)))
    for (c in cell_types) {
        de_group <- Seurat::FindMarkers(object = scdata,
                                        ident.1 = c,
                                        ident.2 = NULL,
                                        only.pos = TRUE,
                                        test.use = "bimod")

        index <- which(cell_types == c)
        list.DE.group[[index]] <- de_group
    }
    
    return(list.DE.group)
}

## Building the signature matrix using Seurat
##
## scdata: The single cell data matrix
## cell_types: A vector of the cell type annotations
## diff.cutoff: The FC cutoff
## pval.cutoff: The p-value cutoff
##
## return: Computed signature matrix

buildSignatureMatrixUsingSeurat <- function(scdata,
                                            cell_types,
                                            diff.cutoff = 0.5,
                                            pval.cutoff = 0.01){
    # perform differential expression analysis
    print("performing differential expression analysis")
    list.DE.group <- DEAnalysis(scdata, cell_types)
    saveRDS(list.DE.group, "list.DE.group.rds")

    num_genes <- c()
    for (c in cell_types) {
        index <- which(cell_types == c)
        de_group <- list.DE.group[[index]]

        DEGenes <- rownames(de_group)[intersect(
            which(de_group$p_val_adj < pval.cutoff),
            which(de_group$avg_log2FC > diff.cutoff)
          )]

        ## why removing those genes ?
        nonMir <- grep("MIR|Mir", DEGenes, invert = TRUE)
        assign(
          paste("cluster_lrTest.table.", c, sep = ""),
          de_group[which(rownames(de_group) %in% DEGenes[nonMir]), ]
        )
        num_genes <- c(num_genes, length(DEGenes[nonMir]))
    }
    
    print("calculating optimal number of genes")
    conditionNumbers <- c()
    for (g in 50:200) {
        Genes <- c()
        j <- 1
        for (i in cell_types) {
          if (num_genes[j] > 0) {
            temp <- paste("cluster_lrTest.table.", i, sep = "")
            temp <- get(temp)
            temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
            Genes <-
              c(Genes, (rownames(temp)[1:min(g, num_genes[j])]))
          }
          j <- j + 1
        }

        ## make signature matrix
        ExprSubset <- scdata_nac_mat[Genes, , drop = FALSE]
        Sig <- NULL
        for (i in as.list(cell_types)) {
          Sig <-
            cbind(Sig, (apply(ExprSubset, 1, function(y) {
              mean(y[names(y) == i])
            })))
        }
        colnames(Sig) <- cell_types
        conditionNumbers <- c(conditionNumbers, kappa(Sig))
    }
    
    # g is optimal gene number
    g <- which.min(conditionNumbers) + min(49, num_genes - 1)
    Genes <- c()
    j <- 1
    for (i in cell_types) {
        if (num_genes[j] > 0) {
          temp <- paste("cluster_lrTest.table.", i, sep = "")
          temp <- get(temp)
          temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
          Genes <-
            c(Genes, (rownames(temp)[1:min(g, num_genes[j])]))
        }
        j <- j + 1
    }

    Genes <- unique(Genes)
    ExprSubset <- scdata_nac_mat[Genes, , drop = FALSE]
    Sig <- NULL
    for (i in as.list(cell_types)) {
      Sig <-
        cbind(Sig, (apply(ExprSubset, 1, function(y) {
          mean(y[names(y) == i])
        })))
    }

    colnames(Sig) <- cell_types
    rownames(Sig) <- Genes
    
    return(Sig)
}

setwd("~/project/git/sigmat/")

library(SingleCellExperiment)
library(dplyr)
library(Seurat)

load("./data/SCE_NAc-n8_tran-etal.rda")

# convert SCE to Seurat object
scdata_nac <- as.Seurat(sce.nac.tran)

# set cell types
Idents(scdata_nac) <- scdata_nac@meta.data$cellType
cell_types <- scdata_nac@meta.data$cellType %>% unique() %>% as.character()

# retrieve expression matrix
scdata_nac_mat <- GetAssayData(scdata_nac, assay = "originalexp", slot = "data")
colnames(scdata_nac_mat) <- scdata_nac@meta.data$cellType

sigmat <- buildSignatureMatrixUsingSeurat(scdata_nac, cell_types)
saveRDS(sigmat, "signature_matrix_nac.rds")
