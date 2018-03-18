setwd("/home/t.pham/opl/R")

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

#03-08-17
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

to_cbioportal_2(txt_file = "../_out/R44-breast-tcga-itraq_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-tcga-itraq_raw/quant-meta.txt",
                    quant_mode = "itraq",
                    type_of_cancer = "brca",
                    study_id = "breast-tcga-mq",
                    study_name = "TCGA breast cancer proteomics study",
                    study_name_short = "TCGA breast cancer proteomics",
                    study_desc = "TCGA breast proteomics processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-crc-tcga_raw/quant.txt",
                    clinical_data = "../_out/R44-crc-tcga_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "coad",
                    study_id = "crc-tcga-mq",
                    study_name = "TCGA colorectal cancer proteomics study",
                    study_name_short = "TCGA colorectal cancer proteomics",
                    study_desc = "TCGA colorectal proteomics processed with MaxQuant")


# 04-08-17
to_cbioportal_2(txt_file = "../_out/R44-ovarian-tcga-itraq_raw/quant.txt",
                    clinical_data = "../_out/R44-ovarian-tcga-itraq_raw/quant-meta.txt",
                    quant_mode = "itraq",
                    type_of_cancer = "ovt",
                    study_id = "ovarian-tcga-mq",
                    study_name = "TCGA ovarian cancer proteomics study",
                    study_name_short = "TCGA ovarian cancer proteomics",
                    study_desc = "TCGA ovarian proteomics processed with MaxQuant")

to_cbioportal_2(txt_file = "../_out/R44-bladder-Latosinska-itraq-lf_itraq/quant.txt",
                    clinical_data = "../_out/R44-bladder-Latosinska-itraq-lf_itraq/quant-meta.txt",
                    quant_mode = "itraq",
                    type_of_cancer = "bladder",
                    study_id = "bladder-Latosinska-itraq",
                    study_name = "Bladder Latosinska itraq",
                    study_name_short = "Bladder Latosinska itraq",
                    study_desc = "Bladder Latosinska itraq processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-bladder-Latosinska-itraq-lf_lf-cid-hcd/quant.txt",
                    clinical_data = "../_out/R44-bladder-Latosinska-itraq-lf_lf-cid-hcd/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "bladder",
                    study_id = "bladder-Latosinska-lf",
                    study_name = "Bladder Latosinska label-free",
                    study_name_short = "Bladder Latosinska label-free",
                    study_desc = "Bladder Latosinska label-free processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-brain-Saratsis_raw/quant.txt",
                    clinical_data = "../_out/R44-brain-Saratsis_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "brain",
                    study_id = "brain-Saratsis",
                    study_name = "Brain Saratsis label-free",
                    study_name_short = "Brain Saratsis label-free",
                    study_desc = "Brain Saratsis label-free processed with MaxQuant")

to_cbioportal_2(txt_file = "../_out/R44-breast-Bouchal-itraq_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-Bouchal-itraq_raw/quant-meta.txt",
                    quant_mode = "itraq",
                    type_of_cancer = "brca",
                    study_id = "breast-Bouchal",
                    study_name = "Breast cancer Bouchal study",
                    study_name_short = "Breast cancer Bouchal",
                    study_desc = "Breast cancer Bouchal study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-breast-DeMarchi_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-DeMarchi_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "brca",
                    study_id = "breast-DeMarchi",
                    study_name = "Breast cancer DeMarchi study",
                    study_name_short = "Breast cancer DeMarchi",
                    study_desc = "Breast cancer DeMarchi study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-breast-Groessl_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-Groessl_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "brca",
                    study_id = "breast-Groessl",
                    study_name = "Breast cancer Groessl study",
                    study_name_short = "Breast cancer Groessl",
                    study_desc = "Breast cancer Groessl study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-breast-Lawrence_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-Lawrence_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "brca",
                    study_id = "breast-Lawrence",
                    study_name = "Breast cancer Lawrence study",
                    study_name_short = "Breast cancer Lawrence",
                    study_desc = "Breast cancer Lawrence study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-breast-Liu_raw/quant.txt",
                    clinical_data = "../_out/R44-breast-Liu_raw/quant-meta.txt",
                    quant_mode = "lf",
                    type_of_cancer = "brca",
                    study_id = "breast-Liu",
                    study_name = "Breast cancer Liu study",
                    study_name_short = "Breast cancer Liu",
                    study_desc = "Breast cancer Liu study processed with MaxQuant")


#07-08-17
to_cbioportal_2(txt_file = "../_out/R44-liver-Xing-itraq8_raw/quant.txt",
                    clinical_data = "../_out/R44-liver-Xing-itraq8_raw/quant-meta.txt",
                    quant_mode = "itraq",
                    type_of_cancer = "liver",
                    study_id = "liver-Xing",
                    study_name = "Liver cancer Xing study",
                    study_name_short = "Liver cancer Xing",
                    study_desc = "Liver cancer Xing study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-crc-Wisniewski_raw/quant.txt",
                    clinical_data = "../_out/R44-crc-Wisniewski_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "coad",
                study_id = "crc-Wisniewski",
                    study_name = "Colorectal cancer Wisniewski study",
                    study_name_short = "Colorectal cancer Wisniewski",
                    study_desc = "Colorectal cancer Wisniewski study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-head-neck-Sepiashvili_raw/quant.txt",
                clinical_data = "../_out/R44-head-neck-Sepiashvili_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "head_neck",
                study_id = "head-neck-Sepiashvili",
                    study_name = "Head and Neck Sepiashvili study",
                    study_name_short = "Head and Neck Sepiashvili",
                    study_desc = "Head and Neck Sepiashvili study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-kidney-Neely_raw/quant.txt",
                clinical_data = "../_out/R44-kidney-Neely_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "kidney",
                study_id = "kidney-Neely",
                study_name = "Kidney Cancer Neely study",
                study_name_short = "Kidney Cancer Neely",
                study_desc = "Kidney Cancer Neely study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-leukemia-Tong_raw/quant.txt",
                clinical_data = "../_out/R44-leukemia-Tong_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "leuk",
                study_id = "leukemia-Tong",
                study_name = "Leukemia Tong study",
                study_name_short = "Leukemia Tong",
                study_desc = "Leukemia Tong study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-lung-Ahn_raw/quant.txt",
                clinical_data = "../_out/R44-lung-Ahn_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "lung",
                study_id = "lung-Ahn",
                study_name = "Lung Ahn study",
                study_name_short = "Lung Ahn",
                study_desc = "Lung Ahn study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-lung-Kikuchi-12_raw/quant.txt",
                clinical_data = "../_out/R44-lung-Kikuchi-12_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "nsclc",
                study_id = "lung-Kikuchi",
                study_name = "Lung Kikuchi study",
                study_name_short = "Lung Kikuchi",
                study_desc = "Lung Kikuchi study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-lung-Kim-PXD002523_raw/quant.txt",
                clinical_data = "../_out/R44-lung-Kim-PXD002523_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "lung",
                study_id = "lung-Kim",
                study_name = "Lung Kim study",
                study_name_short = "Lung Kim",
                study_desc = "Lung Kim study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-lung-Hsu-itraq_raw/quant.txt",
                clinical_data = "../_out/R44-lung-Hsu-itraq_raw/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "lung",
                study_id = "lung-Hsu",
                study_name = "Lung Hsu itraq study",
                study_name_short = "Lung Hsu",
                study_desc = "Lung Hsu itraq study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-lung-Steward-lf-tmt-dia_dda-single-shot/quant.txt",
                clinical_data = "../_out/R44-lung-Steward-lf-tmt-dia_dda-single-shot/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "lusc",
                study_id = "lung-Steward-ss",
                study_name = "Lung Steward single shot proteomics study",
                study_name_short = "Lung Steward single shot",
                study_desc = "Lung Steward single shot proteomics study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-lung-Steward-lf-tmt-dia_lf/quant.txt",
                clinical_data = "../_out/R44-lung-Steward-lf-tmt-dia_lf/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "lusc",
                study_id = "lung-Steward-lf",
                study_name = "Lung Steward 12-fraction label-free proteomics study",
                study_name_short = "Lung Steward 12-fraction label-free",
                study_desc = "Lung Steward 12-fraction label-free proteomics study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-lung-Steward-lf-tmt-dia_tmt/quant.txt",
                clinical_data = "../_out/R44-lung-Steward-lf-tmt-dia_tmt/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "lusc",
                study_id = "lung-Steward-tmt",
                study_name = "Lung Steward TMT proteomics study",
                study_name_short = "Lung Steward TMT",
                study_desc = "Lung Steward TMT proteomics study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-ovarian-Hughes-tmt_tmt10/quant.txt",
                clinical_data = "../_out/R44-ovarian-Hughes-tmt_tmt10/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "ovary",
                study_id = "ovary-Hughes-tmt10-FFPE",
                study_name = "Ovary Hughes FFPE proteomics study",
                study_name_short = "Ovary Hughes FFPE",
                study_desc = "Ovary Hughes FFPE TMT10 study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-ovarian-Hughes-tmt_tmt8/quant.txt",
                clinical_data = "../_out/R44-ovarian-Hughes-tmt_tmt8/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "ovary",
                study_id = "ovary-Hughes-tmt8",
                study_name = "Ovary Hughes fronzen sample study",
                study_name_short = "Ovary Hughes fronzen sample",
                study_desc = "Ovary Hughes fronzen sample TMT8 study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-ovarian-Hughes-tmt_tmt10_2/quant.txt",
                clinical_data = "../_out/R44-ovarian-Hughes-tmt_tmt10_2/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "ovary",
                study_id = "ovary-Hughes-tmt10_2",
                study_name = "Ovary Hughes 3-subtype study",
                study_name_short = "Ovary Hughes 3-subtype",
                study_desc = "Ovary Hughes 3-subtype TMT10 study processed with MaxQuant")




to_cbioportal_2(txt_file = "../_out/R44-cll-Eagle-itraq_raw-wiff/quant.txt",
                clinical_data = "../_out/R44-cll-Eagle-itraq_raw-wiff/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "cll",
                study_id = "cll-Eagle",
                study_name = "Chronic Lymphocytic Leukemia Eagle study",
                study_name_short = "Chronic Lymphocytic Leukemia Eagle",
                study_desc = "Chronic Lymphocytic Leukemia Eagle study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-retinoblastoma-Danda-itraq_raw/quant.txt",
                clinical_data = "../_out/R44-retinoblastoma-Danda-itraq_raw/quant-meta.txt",
                quant_mode = "itraq",
                type_of_cancer = "rbl",
                study_id = "rbl-Danda",
                study_name = "Retinoblastoma Danda study",
                study_name_short = "Retinoblastoma Danda",
                study_desc = "Retinoblastoma Danda study processed with MaxQuant")



to_cbioportal_2(txt_file = "../_out/R44-pancreas-Padden_raw/quant.txt",
                clinical_data = "../_out/R44-pancreas-Padden_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "pancreas",
                study_id = "pancreas-Padden",
                study_name = "Pancreas Padden study",
                study_name_short = "Pancreas Padden",
                study_desc = "Pancreas Padden study processed with MaxQuant")


to_cbioportal_2(txt_file = "../_out/R44-skin-Welinder_raw/quant.txt",
                clinical_data = "../_out/R44-skin-Welinder_raw/quant-meta.txt",
                quant_mode = "lf",
                type_of_cancer = "skin",
                study_id = "skin-Welinder",
                study_name = "Skin Welinder study",
                study_name_short = "Skin Welinder",
                study_desc = "Skin Welinder study processed with MaxQuant")
