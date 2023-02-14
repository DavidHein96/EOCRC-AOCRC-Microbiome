
# ------------ Libraries & Setup ------------ #

library(phyloseq)
library(ANCOMBC)
library(vegan)
library(plyr)
library(microbiomeMarker)
library(ggpubr)
library(tidyverse)
library(Maaslin2)
library(ggtext)
set.seed(2023)



# ------------ Function to create phyloseq object ------------ #

build_phylo_obj <- function(selected_patients, counts_all) {
  # selected patients is semi-prefiltered clinical data with only 
    # select variables and filtered patients, needs SampleID column, counts_all is full ASV table with patients as columns
  
  obj_list <- list()
  
  # Get only patients that we want to include
  counts_all <- counts_all %>%dplyr::select("Kingdom","Phylum","Class","Order","Family","Genus","Species",any_of(selected_patients$SampleID))
  
  # Filter out Archea and ASVs that have all zeros
  counts_all <- counts_all %>%filter(Kingdom !="Archaea")%>%filter(rowSums(across(where(is.numeric)))!=0)
  
  # add in signifier of taxa level
  counts_all<-counts_all%>%mutate(Kingdom=paste0("d__",Kingdom),Phylum=paste0("p__",Phylum),Class=paste0("c__",Class),Order=paste0("o__",Order),
                                  Family=paste0("f__",Family),Genus=paste0("g__",Genus),Species=paste0("s__",Species))
  
  # add up the identical ASVs
  # first make a taxonomy variable in a single column and remove the multiple columns
  counts_all_combinedASV <- counts_all %>% mutate(taxonomy = base::paste(Kingdom,Phylum,Class,Order,Family,Genus,Species, sep=";")) %>% relocate(taxonomy)
  counts_all_combinedASV2 <- counts_all_combinedASV %>% dplyr::select(-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species)
  
  # use ddply to sum up the identical taxa
  counts_all_combinedASV3<-counts_all_combinedASV2 %>%ddply("taxonomy",numcolwise(sum))
  
  # rejoin the taxa names and reorder
  taxa_names<- counts_all_combinedASV%>%dplyr::select(taxonomy,Kingdom,Phylum,Class,Order,Family,Genus,Species)%>%distinct()
  taxa_order <- counts_all_combinedASV3$taxonomy
  taxa_names_ordered <- taxa_names[match(taxa_order,taxa_names$taxonomy),]
  
  # Build Phyloseq object
  asvs <- as.matrix(counts_all_combinedASV3[,-1])
  row.names(asvs)<-counts_all_combinedASV3$taxonomy
  otu<-otu_table(asvs,taxa_are_rows =TRUE)
  
  clin <- selected_patients%>%dplyr::select(-SampleID)
  row.names(clin)<-selected_patients$SampleID
  samples<-sample_data(clin)
  
  tax <- as.matrix(taxa_names_ordered[,-1])
  row.names(tax) <- taxa_names_ordered$taxonomy
  tax<-tax_table(tax)
  
  obj<-phyloseq(otu,samples,tax)
  
  # 10% prev filter
  obj_filt <- filter_taxa(obj,function(x) sum(x>0)>nrow(clin)*0.1,TRUE)
  
  # return both filtered and un filt phylo objs
  obj_list[["unfilt"]]<-obj
  obj_list[["filt"]]<-obj_filt
  
  return(obj_list)
  
}



# ------------ Function to run ancombc2 on all taxa levels ------------ #
run_ancom_all_levels <- function(phylo_obj,variable){
  
  ancom_results<-data.frame()
  tax_levels<-c("Phylum","Class","Order","Family","Genus","Species")
  
  # runs ancombc at all levels and binds results
  for (tax_level in tax_levels){
    
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    
    a4<-ancombc(obj_sub,formula=paste0(variable),prv_cut = 0)
    res<-data.frame(a4$res)
    ancom_results<-rbind(ancom_results,res)
  }
  
  # adjust p value for multiple testing
  ancom_bhcor_p<-p.adjust(ancom_results[,12],method = "BH")
  ancom_results$ancom_bhcor_p<-ancom_bhcor_p
  return(ancom_results)
}



# ------------ Function to run MaAsLin2 on all taxa levels ------------ #
run_mas_all_levels <- function(phylo_obj,variable){
  
  mas_results<-data.frame()
  tax_levels<-c("Phylum","Class","Order","Family","Genus","Species")
  
  for (tax_level in tax_levels){
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run maaslin2 and bind results
    fit_mas <- Maaslin2(input_data = df_input_data,input_metadata = df_meta,output = "test",fixed_effects = c(paste0(variable)), min_prevalence = 0.0,
                        plot_heatmap = FALSE,plot_scatter = FALSE)
    res<-data.frame(fit_mas$results)
    mas_results<-rbind(mas_results,res)
  }
  
  # adjust p values for multiple testing
  mas_bhcor_p<-p.adjust(mas_results[,6],method = "BH")
  mas_results$mas_bhcor_p<-mas_bhcor_p
  return(mas_results)
}



# ------------ function to combine mas and ancom result into a single df ------------ #
combine_res <- function(mas_res,ancom_res){
  
  # select and rename relevant columns so can bind them
  m <- mas_res %>% select(feature, value, coef, stderr, pval, mas_bhcor_p) %>% dplyr::rename(taxa = feature, bhcorp = mas_bhcor_p)%>%mutate(test="mas")
  colnames(ancom_res)[1] = "taxa"
  colnames(ancom_res)[3] = "coef"
  colnames(ancom_res)[6] = "stderr"
  colnames(ancom_res)[12] = "pval"
  colnames(ancom_res)[19] = "bhcorp"
  
  ancom_res<-ancom_res%>%mutate(test="ancom", value = colnames(ancom_res)[9])
  a <- ancom_res%>% select(taxa, value, coef, stderr, pval, bhcorp, test)
  
  combined <- rbind(a,m)
  combined <- combined %>% mutate(taxa = str_replace_all(taxa, "\\."," "),
                     taxa = str_replace_all(taxa, "-"," "),
                     taxa = str_replace_all(taxa, "\\]"," "),
                     taxa = str_replace_all(taxa, "\\["," "))
  
  # pivot wider
  overall_results_wide <- combined%>%pivot_wider(names_from = test,values_from = c(coef,stderr,pval,bhcorp,value))
  
  # Use this table for the supplemental tables
  significant_results<-overall_results_wide%>%#filter(bhcorp_ancom<0.5 & bhcorp_mas<0.5)%>%
                                      select(taxa,coef_ancom,stderr_ancom,bhcorp_ancom,coef_mas,stderr_mas,bhcorp_mas)%>%
                                      mutate(taxa = ifelse( str_detect(taxa,"^d__"), str_extract(taxa,"g__.+"),taxa))%>%
                                      arrange( (bhcorp_ancom+bhcorp_mas)/2)
                  
  #all_results<-overall_results_wide%>%select(taxa,coef_ancom,stderr_ancom,bhcorp_ancom,coef_mas,stderr_mas,bhcorp_mas)
  
  #combined_res = list()
  #combined_res[["all_results"]]<-all_results
  #combined_res[["significant_results"]]<-significant_results
  return(significant_results)
}


