
# create
tp.write <- function(d, file) {
    write.table(d,
                file = file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
}

opl.norm <- function(d) {

    total <- apply(d, 2, sum)
    m <- mean(total)
    factor <- total / m;

    tmp <-  d / (matrix(1, nrow(d), 1) %*% factor)

    return(tmp)

}

# collapse to first gene symbol
semicolon_split <- function(comma_separated_text) {
    return(gsub(" ", "", unlist(strsplit(paste0(" ", comma_separated_text), ";"))))
}

collapse_to_first_gene <- function(quant_mode = "lf",
                                   gene.names,             # list of genes
                                   d.final,                # data table
                                   tc = colSums(d.final)   # total count used for normalization
                                   ) {

    gn <- gene.names

    for (i in 1:length(gn)) {
        gn[i] <- semicolon_split(gn[i])[1]
    }

    gn.1 <- unique(gn)

    dd.1 <- NULL

    if (quant_mode == "lf") {
        for (i in 1:length(gn.1)) {
            ind <- which(gn == gn.1[i])
            dd.1 <- rbind(dd.1, colSums(d.final[ind,]))
        }
    } else {
        for (i in 1:length(gn.1)) {
            ind <- which(gn == gn.1[i])
            dd.1 <- rbind(dd.1, apply(d.final[ind,], 2, FUN = median, na.rm = TRUE))
        }
    }

    row_ok <- gn.1 != ""
    
    rownames(dd.1)[row_ok] <- gn.1[row_ok]

    #dd.1 <- dd.1[rowSums(dd.1) > 0,]


    if (quant_mode == "lf") {

        m <- mean(tc)
        factor <- tc / m

        tmp <-  dd.1 / (matrix(1, nrow(dd.1), 1) %*% factor)
    } else {
        tmp <- dd.1
    }

    return(list(dd = dd.1[row_ok,], dd.norm = tmp[row_ok,]))

}


# cd ~/cbioportal/cbioportal/core/src/main/scripts/importer
# ./validateData.py -s ~/opl/R/crc1 -n -html ~/opl/R/test-import.html
# ./metaImport.py -s ~/opl/R/crc1 -n -o

# to remove
# ./cbioportalImporter.py -c remove-study -meta ~/cbioportal/R/crc-r44/meta_study.txt

to_cbioportal_2 <- function(txt_file,
                            clinical_data,
                            quant_mode,
                            type_of_cancer,
                            study_id,
                            study_name,
                            study_name_short,
                            study_desc,
                            out_dir = study_id) {

    tab <- read.delim(txt_file, header = TRUE)

    gn <- gsub("[\"=]", "", as.character(tab[, "HGNC"]))

    quant_col <- colnames(tab)[3:ncol(tab)]

    tmp <- collapse_to_first_gene(quant_mode, gn, tab[, quant_col]) 
        
    d <- tmp$dd.norm

    colnames(d) <- paste0(study_id, "-", colnames(d))

    
    

    gn2 <- paste0(rownames(d), "|", rownames(d))

    
    dir.create(out_dir)

    # ----- write meta study
    
    cat("type_of_cancer:", type_of_cancer,
"\ncancer_study_identifier: ", study_id,
"\nname: ", study_name,
"\nshort_name: ", study_name_short,
"\ndescription: ", study_desc,
"\ngroups: PUBLIC\n",
sep = "",
file = paste0(out_dir, "/", "meta_study.txt"))

  cat("Meta study file written.\n")
    
    #----- write normalized count data

    tp.write(cbind(Composite.Element.REF = gn2, d),
             paste0(out_dir, "/", "data_quant.txt"))

    cat("cancer_study_identifier: ", study_id, "\n",
        "genetic_alteration_type: PROTEIN_LEVEL
datatype: CONTINUOUS
stable_id: protein_quantification
profile_description: Protein levels (mass spectrometry)
show_profile_in_analysis_tab: true
profile_name: Protein levels (mass spectrometry)
data_filename: ", "data_quant.txt", "\n",
sep = "",
file = paste0(out_dir, "/", "meta_quant.txt"))


    cat("Data files written.\n")


    #----- write z-score data

    tp.write(cbind(Composite.Element.REF = gn2, t(scale(t(d)))),
             paste0(out_dir, "/", "data_quant_Zscores.txt"))

    cat("cancer_study_identifier: ", study_id, "\n",
        "genetic_alteration_type: PROTEIN_LEVEL
datatype: Z-SCORE
stable_id: protein_quantification_zscores
profile_description: Protein levels (mass spectrometry) Z-score
show_profile_in_analysis_tab: true
profile_name: Protein levels (mass spectrometry) Z-score
data_filename: ", "data_quant_Zscores.txt", "\n",
sep = "",
file = paste0(out_dir, "/", "meta_quant_Zscores.txt"))


    cat("Z-score data files written.\n")

    #----- write clinical data

    cat("cancer_study_identifier: ", study_id,
        "\ngenetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_samples.txt\n",
sep = "",
file = paste0(out_dir, "/", "meta_samples.txt"))


    # data
    
    #Row 1: The attribute Display Names: The display name for each clinical attribute
    #Row 2: The attribute Descriptions: Long(er) description of each clinical attribute
    #Row 3: The attribute Datatype: The datatype of each clinical attribute (must be one of: STRING, NUMBER, BOOLEAN)
    #Row 4: The attribute Priority: A number which indicates the importance of each attribute. In the future, higher priority attributes will appear in more prominent places than lower priority ones on relevant pages (such as the Study View). A lower number indicates a higher priority.

#    cat("#Patient Identifier\tSample Identifier\tSample Type\tOncotree Code\tCancer Type\tCancer Type Detailed
#Patient Identifier\tSample Identifier\tThe type of sample (i.e., normal, primary, met, recurrence)\tOncotree code\tCancer type\tCancer type detailed
#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING
#1\t1\t1\t1\t1\t1
#PATIENT_ID\tSAMPLE_ID\tSAMPLE_TYPE\tONCOTREE_CODE\tCANCER_TYPE\tCANCER_TYPE_DETAILED\n",
#sep = "",
#file = paste0(out_dir, "/", "data_samples.txt"))

    tab2  <- read.delim(clinical_data, header = TRUE)
    rownames(tab2) <- paste0(study_id, "-", as.character(tab2[,1]))

    # row 1
    cat("#Patient Identifier\tSample Identifier",
        sep = "",
        file = paste0(out_dir, "/", "data_samples.txt"))
    for (i in 2:ncol(tab2)) {
        cat("\t", colnames(tab2)[i],
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
    }

    # row 2
    cat("\n#Patient Identifier\tSample Identifier",
        sep = "",
        file = paste0(out_dir, "/", "data_samples.txt"),
        append = TRUE)
    for (i in 2:ncol(tab2)) {
        cat("\t", colnames(tab2)[i],
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
    }

    # row 3
    cat("\n#STRING\tSTRING",
        sep = "",
        file = paste0(out_dir, "/", "data_samples.txt"),
        append = TRUE)
    for (i in 2:ncol(tab2)) {
        cat("\t", "STRING",
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
    }

    # row 4
    cat("\n#1\t1",
        sep = "",
        file = paste0(out_dir, "/", "data_samples.txt"),
        append = TRUE)
    for (i in 2:ncol(tab2)) {
        cat("\t", "1",
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
    }

    # row 5
    cat("\nPATIENT_ID\tSAMPLE_ID",
        sep = "",
        file = paste0(out_dir, "/", "data_samples.txt"),
        append = TRUE)
    for (i in 2:ncol(tab2)) {
        cat("\t", toupper(colnames(tab2)[i]),
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
    }

    # real data
    for (i in 1:ncol(d)) {
        cat("\n", colnames(d)[i], "\t", colnames(d)[i],
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)

        for (j in 2:ncol(tab2)) {
            cat("\t", as.character(tab2[colnames(d)[i], j]),
            sep = "",
            file = paste0(out_dir, "/", "data_samples.txt"),
            append = TRUE)
        }
    }
    
    cat("Clinical info written.\n")
    

    #----- case list
    dir.create(paste0(out_dir, "/", "case_lists"))

    cat("cancer_study_identifier: ", study_id,
        "\nstable_id: ", study_id, "_all",
        "\ncase_list_name: All tumors
case_list_description: All tumor samples (", ncol(d)," samples)
case_list_category: all_cases_in_study
case_list_ids: ", paste0(colnames(d), collapse = "\t"),
            sep = "",
            file = paste0(out_dir, "/", "case_lists/cases_all.txt"))
}



to_cbioportal_2(txt_file = "./quant.txt",
                    clinical_data = "./quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "coad",
                    study_id = "stool 291",
                    study_name = "Colorectal cancer stool proteomics study",
                    study_name_short = "Stool proteomics",
                    study_desc = "Stool proteomics 291 samples")



