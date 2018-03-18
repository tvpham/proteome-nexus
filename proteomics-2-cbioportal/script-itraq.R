setwd("~/opl/_out/R44-breast-Bouchal-itraq_raw")

# check Fig. 1, 3 runs of itraq 8plex.


tp.write <- function(d, file) {
    write.table(d,
                file = file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}

log2_na <- function(d) {
    d_log2 <- log2(d)
    d_log2[d_log2 == -Inf]  <- NA
    return(d_log2)
}


tab <- read.delim("opl-protein-report.txt", header = TRUE)

reporter <- grep("^Reporter.intensity.corrected...", colnames(tab))

colnames(tab)[reporter]

# all samples
{
    selected <- colnames(tab)[reporter]

    n <- unique(substring(selected, 32, 100))


    d_final  <- cbind(Uniprot = as.character(tab[, "Protein.IDs"]), HGNC = gsub("[\"=]", "", as.character(tab[, "HGNCplus.gene.names"])))

    for (ms_run in c("E1", "E2", "E3")) {

        ref_columns <- paste0("Reporter.intensity.corrected.", 0:7, ".", ms_run)
        d_ref <- log2_na(tab[, ref_columns])
        d_ref_mean <- rowMeans(d_ref, na.rm = TRUE)
        
        for (s in 0:7) {
            column <- paste0("Reporter.intensity.corrected.", s, ".", ms_run)

            vec <- log2_na(tab[, column])
            
            for (i in 1:length(vec)) {
                if (!is.na(vec[i])) {
                    if (!is.na(d_ref_mean[i])) {
                        vec[i] <- vec[i] - d_ref_mean[i]
                    } else {
                        stop("why")
                    }
                }
            }

            tmp <- matrix(vec, ncol=1)

            tmp_name <- paste0("x.", s, ".", substring(ms_run, 1, 20))

            while (tmp_name %in% colnames(d_final)) {
                cat("duplicate ", tmp_name, "\n")
                tmp_name <- paste0(tmp_name, "_")
            }
        
            
            colnames(tmp) <- tmp_name
        
            d_final <- cbind(d_final, tmp) 
        }
        cat(ms_run,"\n")
    }

}



tp.write(d_final, "quant.txt")

tp.write(cbind(Sample = colnames(d_final)[3:ncol(d_final)],
               Sample_type = rep("Tissue", ncol(d_final)-2),
               Oncotree_code = rep("BRCA", ncol(d_final)-2),
               Cancer_type = rep("Breast Cancer", ncol(d_final)-2),
               Cancer_type_detailed = rep("Breast Cancer", ncol(d_final)-2),
               Sub_type = rep("Breast cancer", ncol(d_final)-2)),
         file = "quant-meta-tmp.txt")

