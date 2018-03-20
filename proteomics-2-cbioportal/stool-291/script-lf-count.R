
tp.write <- function(d, file) {
    write.table(d,
                file = file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}


tab <- read.delim("M17-1068_Supplement2-tab2-normalized-count.txt", header = TRUE)

msms <- 2:ncol(tab)

d <- tab[, msms]

d_final  <- cbind(Uniprot = gsub("[\"=]", "", as.character(tab[, 1])),
                  HGNC = gsub("[\"=]", "", as.character(tab[, 1])),
                  d)

tp.write(d_final, "quant.txt")

cancer_types <- colnames(d)
cancer_types[grep("^Control.Sample.", colnames(d))] <- "Control"
cancer_types[grep("^Adenoma", colnames(d))] <- "Adenoma"
cancer_types[grep("^Advanced", colnames(d))] <- "Advanced Adenoma"
cancer_types[grep("^Colo", colnames(d))] <- "Colorectal Cancer"

tp.write(cbind(Sample = colnames(d_final)[3:ncol(d_final)],
               Sample_type = rep("Tissue", ncol(d_final)-2),
               Oncotree_code = rep("COAD", ncol(d_final)-2),
               Cancer_type = rep("Colorectal Cancer", ncol(d_final)-2),
               Cancer_type_detailed = rep("Colorectal Cancer", ncol(d_final)-2),
               Sub_type = cancer_types),
         file = "quant-meta.txt")

