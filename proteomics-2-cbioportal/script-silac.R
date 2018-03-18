setwd("~/opl/_out/R44-breast-Tyanova-silac")

tp.write <- function(d, file) {
    write.table(d,
                file = file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}


tab <- read.delim("opl-protein-report.txt", header = TRUE)

ratio <- grep("^Ratio.H.L.normalized.", colnames(tab))

d <- tab[, ratio]

colnames(d) <- paste0("x.", substring(colnames(d), 22, 1000))

d_final  <- cbind(Uniprot = as.character(tab[, "Protein.IDs"]),
                  HGNC = gsub("[\"=]", "", as.character(tab[, "HGNCplus.gene.names"])),
                  d)

tp.write(d_final, "quant.txt")

tp.write(cbind(Sample = colnames(d_final)[3:ncol(d_final)],
               Sample_type = rep("Tissue", ncol(d_final)-2),
               Oncotree_code = rep("BRCA", ncol(d_final)-2),
               Cancer_type = rep("Breast Cancer", ncol(d_final)-2),
               Cancer_type_detailed = rep("Breast Cancer", ncol(d_final)-2),
               Sub_type = rep("[Not Available]", ncol(d_final)-2)),
         file = "quant-meta-tmp.txt")

