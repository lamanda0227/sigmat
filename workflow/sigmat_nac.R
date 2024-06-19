## Performing DE analysis using Seurat
##
## scdata: The single cell data matrix
## cell_types: A vector of the cell type annotations
##
## return: List with the differentially expressed genes for each cell type

DEAnalysis <- function(scdata, cell_types) {
    print(sprintf("%s cell types available for DE analysis", length(cell_types)))
    print(cell_types)
    
    scdata_cell_types <- scdata@meta.data$cellType %>% unique() %>% as.character()
    print(sprintf("%s cell types available from Seurat object", length(scdata_cell_types)))
    print(scdata_cell_types)
    
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
## out_dir: Output directry for saving intermediary outputs
## diff.cutoff: The FC cutoff
## pval.cutoff: The p-value cutoff
##
## return: Computed signature matrix

buildSignatureMatrixUsingSeurat <- function(scdata,
                                            cell_types,
                                            out_dir = ".",
                                            diff.cutoff = 0.5,
                                            pval.cutoff = 0.01){
    # perform differential expression analysis
    print("==============================")
    print("performing differential expression analysis")
    print(format(Sys.time(), "%H:%M:%S"))
    list.DE.group <- DEAnalysis(scdata, cell_types)
    saveRDS(list.DE.group, sprintf("%s/list.DE.group.rds", out_dir))
    
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
    
    print("==============================")
    print("calculating optimal number of genes")
    print(format(Sys.time(), "%H:%M:%S"))

    scdata_nac_mat <- GetAssayData(scdata, assay = "originalexp", layer = "data")
    colnames(scdata_nac_mat) <- scdata@meta.data$cellType
    
    conditionNumbers <- c()
    for (g in 50:200) {
        Genes <- c()
        j <- 1
        for (i in cell_types) {
          if (num_genes[j] > 0) {
            temp <- paste("cluster_lrTest.table.", i, sep = "")
            temp <- get(temp)
            temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
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
    
    print("==============================")
    print("generating signature matrix")
    print(format(Sys.time(), "%H:%M:%S"))
    
    # g is optimal gene number
    g <- which.min(conditionNumbers) + min(49, num_genes - 1)
    Genes <- c()
    j <- 1
    for (i in cell_types) {
        if (num_genes[j] > 0) {
          temp <- paste("cluster_lrTest.table.", i, sep = "")
          temp <- get(temp)
          temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
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

# remove unwanted cell types
unwanted_cell_types <- c("drop.lowNTx", "drop.doublet_A", "drop.doublet_B", "drop.doublet_D", "drop.doublet_C")
scdata_nac_filtered <- subset(scdata_nac, idents = setdiff(Idents(scdata_nac), unwanted_cell_types))

# set cell types for filtered data
Idents(scdata_nac_filtered) <- scdata_nac_filtered@meta.data$cellType
cell_types_filtered <- scdata_nac_filtered@meta.data$cellType %>% unique() %>% as.character()

# calculating signature matrix
out_dir <- "~/project/git/sigmat/data"
cell_types <- scdata_nac_filtered@meta.data$cellType %>% unique() %>% as.character()
sigmat <- buildSignatureMatrixUsingSeurat(scdata_nac_filtered, cell_types, out_dir)
saveRDS(sigmat, sprintf("%s/signature_matrix_nac.rds", out_dir))
