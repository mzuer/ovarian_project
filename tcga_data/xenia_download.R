

# download 05.09.2021
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_GTEX_category.txt
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz
# from
# https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=http%3A%2F%2F127.0.0.1%3A7222


annot_dt <- read.delim("../tcga_data/TCGA_GTEX_category.txt", header=TRUE, stringsAsFactors = FALSE)
annot_dt <- annot_dt[grepl("Ovar", annot_dt$TCGA_GTEX_main_category),]



suppressMessages(library(UCSCXenaTools))
suppressMessages(library(dplyr))
tcga_ov_cohort = XenaData %>%
  filter(XenaHostNames == "toil") %>% # select TCGA Hub
  XenaScan("TCGA Ovarian Serous Cystadenocarcinoma")   # select LUAD cohort
tcga_ov_cohort


XenaGenerate(subset = XenaHostNames=="toilHub") %>% 
  
  XenaFilter(filterDatasets = "TcgaTargetGtex_rsem_gene") 
  
  XenaScan("Ovary") %>% 
  XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo

XenaData %>%
  filter(XenaHostNames == "Toil") %>% # select TCGA Hub
  XenaScan("TCGA Ovarian")


XenaData  %>% XenaScan("Toil")  

  XenaFilter(filterDatasets = "Ovary")

%>% 
  XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo
