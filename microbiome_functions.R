#Microbiome Tools 

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(scales)
library(phyloseq)
library(cowplot)

# Function: Create phyloseq object ----------------------------------------
physeq.create <- function(asvdf, #Dataframe; Rownames = sequences (e.g CCTGTT...) . Columns = sequencing counts.
                          taxdf, #Dataframe; Rownames = sequences (e.g CCTGTT...). Columns = Sequence taxonomic identification (Kingdom, Phylum, Class, Order, Family Genus)
                          mapdf, #Dataframe; No specific rownames required. Rows = Samples, Columns = Sample metadata (e.g depression scores, sex, age, etc.)
                          ID.col, #name of the column containing study IDs 
                          export.asvsp = TRUE, #export CSV file that maps phyloseq assigned species ID (sp####) to ASV sequences (CCTGGT...)
                          export.wd  = get_wd() #directory exports will be saved
                          )
  {
  
  #=== Wrangling === 
  #Standardize ID column 
  mapdf <- mapdf %>% dplyr::rename(StudyID = paste(ID.col))
  
  #=== asvdf === 
  #Remove asvdf rownames
  seqs = rownames(asvdf)
  rownames(asvdf) = NULL
  
  # === Taxdf  === 
  #Unit test: taxdf is aligned with asvdf
  qc.physeq1 <- all(rownames(taxdf) == seqs) #pass = TRUE
  if(qc.physeq1 == FALSE){
    print("Error: taxdf and asvdf rows are not aligned")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  
  
  #Remove taxdf rownames for physeq formatting
  taxdf.asv <- taxdf #save ASV key
  rownames(taxdf) = NULL #remove rows
  
  # Unit test: 
  #Taxa nrow = ASV seq nrow 
  qc.physeq2 <- nrow(taxdf) == nrow(asvdf) #pass = TRUE
  if(qc.physeq2 == FALSE){
    print("Error: taxdf and asvdf have different number of rows")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  
  
  
  
  #=== mapdf === 
  #Unit test: 
  #nrow of mapdf$StudyID is equal to ncol of asvdf
  #Checks that there are same number of StudyIDs in mapdf and asvdf
  qc.physeq3 <- nrow(mapdf) == ncol(asvdf)  #pass = TRUE
  if(qc.physeq3 == FALSE){
    print("Error: mismatched number of samples between asvdf and mapdf" )
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  
  #Unit test: 
  #Set of StudyIDs in asvdf (cols) are also in mapdf (rows)
  qc.physeq4 <- all(colnames(asvdf) %in% mapdf$StudyID) %>% all #pass = TRUE
  qc.physeq5 <- all(mapdf$StudyID %in% colnames(asvdf)) %>% all #pass = TRUE
  if(qc.physeq4 == FALSE | qc.physeq5 == FALSE){
    print("Error: mismatched sample IDs between asvdf and mapdf" )
    print("Samples in asvdf but not mapdf:")
    qc.df <- data.frame(colnames(asvdf)) %>% filter(!(colnames.asvdf. %in% mapdf$StudyID)) #show all samples in asvdf except those shared between asvdf and mapdf
    print(qc.df)
    
    print("Samples in mapdf but not asvdf:")
    qc.df <- data.frame(mapdf$StudyID) %>% filter(!(mapdf.StudyID %in% colnames(asvdf))) #show all samples in mapdf except those shared between asvdf and mapdf
    
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  #StudyID becomes rownames
  rownames(mapdf) = mapdf$StudyID
  
  
  #Unit test: 
  #Set of StudyIDs in asvdf (cols) are also in mapdf (rownames)
  qc.physeq6 <- all(rownames(mapdf) %in% colnames(asvdf)) %>% all #pass = TRUE
  qc.physeq7 <- all(colnames(asvdf) %in% rownames(mapdf)) %>% all #pass = TRUE
  if(qc.physeq6 == FALSE | qc.physeq7 == FALSE){
    print("Error: set of asvdf colnames is different from set of mapdf rownames" )
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  #=== Create phyloseq object === 
  # Convert the ASV and taxonomy data frames to matrices. 
  tax_mat = as.matrix(taxdf)
  asv_mat = as.matrix(asvdf)
  
  dat = phyloseq(otu_table(asv_mat, taxa_are_rows = TRUE), tax_table(tax_mat), sample_data(mapdf))
  str(dat)
  
  
  

  # === Export ASV -> phyloseq sp ID# conversion key === 
  taxdf.asv #taxa df without rows (ASV sequences) removed 
  
  if(export.asvsp == TRUE){
    #Unit test: 
    #Physeq and taxdf taxonomic assignments are in the same order
    qc.asv_sp <- (dat %>% tax_table() %>% data.frame() %>% replace(is.na(.), "UNKNOWN") == taxdf %>% replace(is.na(.), "UNKNOWN")) %>% 
      as.logical %>% as.data.frame() %>% rename(qc = ".") 
    if(unique(qc.asv_sp) %>% length != 1 | unique(qc.asv_sp)[1,1] != TRUE){
      print("Error: taxdf.asv and phyloseq object are not aligned.")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
      
  }
    
    #Create asv_sp key
    print(paste("Phyloseq taxa table (dat) and ASV sequence taxa df (taxdf.asv) aligned with original taxa df (taxdf) :", qc))
    print("Create ASV -> phyloseq sp ID# conversion key")
    sp_asv <- (taxdf.asv %>% rownames) %>% data.frame() %>% rename(asv_seq = ".") %>% mutate(sp = dat %>% tax_table %>% rownames())
    
    #Export
    write.csv(sp_asv, "asv_spCode.csv")
  }
  
  #=== Remove host mitochondrial sequences ===
  dat2 = subset_taxa(dat,
                     !(Kingdom == 'Eukaryota' |
                         is.na(Phylum) |
                         (!is.na(Family) &
                            Family == 'Mitochondria')))
  
  
  # === Return phyloseq object === 
  return(dat2)
}



# Function: Add taxa names to sp ID# ------------------------------------------------
#Get sp ID and taxonomic identification from phyloseq object (df)

sp_taxa_ID <- function(physeq, #phyloseq object
                       taxa.col = NA #Character vector; taxonomic levels to be included in phyloseq taxa ID (e.g c("Family", "Genus")). If no input given, defaults to all levels in physeq object
                       )
{
  asvtab <- psmelt(physeq) %>% data.frame

  #Set taxa names if not given
  if(is.na(taxa.col)){
    taxa.col <- rank_names(physeq)
  }

  #List of phyloseq sp IDs and their taxonomic assignment
  asvtaxa <- asvtab[,c("OTU", taxa.col)] %>% unique() %>% 
    mutate(tname = paste(OTU, Kingdom, Phylum, Class, Order, Family, Genus, sep = "_"), #Create new col tname with OTU and full taxonomic details 
           OTU_order = parse_number(OTU) %>% as.numeric()) %>%
    arrange(OTU_order) %>% dplyr::select(-OTU_order)
  
  #Unit test: 
  #asvtaxa ASV (colname = OTU) aligns with phyloseq object ASV rownames 
  qc.taxa_names <- (asvtaxa$OTU == tax_table(physeq) %>% row.names()) %>% data.frame() %>% filter(FALSE %in% .) %>% nrow()
  if(qc.taxa_names != 0){
    print("QC passed, full ASV taxonomic names will be added to phyloseq object")
  }
  
  #Add new taxa names to phyloseq object
  taxa_names(physeq) <- asvtaxa$tname
  print(head(asvtaxa$tname)) #print example of new taxa names for user
  
  
  return(physeq)
  }



# Functions: relative abundance   -------------------------------------------
#phyloseq counts -> phyloseq relative abundance and wide ASV dataframe
physeq.convert.rel <- function(physeq, warning.in = TRUE)
  {
  dat_rel <- transform_sample_counts(physeq,function(x) x/sum(x))
  
  rel <- physeqrel.dataframe(dat_rel, warning.in)
  return(list(dat_rel, rel))
}



#phyloseq relative abundance -> wide ASV dataframe
physeqrel.dataframe <- function(physeq, #relative abundance phyloseq object
                                warning = TRUE #if relative abundance is not 1, issue a warning instead of killing script. TRUE = issue warning, FALSE = kill script.
                                )
  #input relative abundance phyloseq object
  #ouput relative abundance ASV table (rel)
  {
  
  #Unit test: is it relative abundance data? 
  qc <- physeq %>% otu_table %>% data.frame %>% pull(1) %>% sum 
  if(qc > 1){
    print("ERROR: first sample looks like count data. Input must be relative abundance data")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  
  
  #Psmelt relative abundanc physeq 
  rel <- psmelt(physeq)
  
  
  #Unit test: 
  #Total relative abundance sums to 1
  #Returns warning. Script will continue even if test fails, unless warning = FALSE 
  qc <- rel %>% group_by(Sample) %>% dplyr::summarize(total_abundance = sum(Abundance))
  if(min(qc$total_abundance) != 1 | max(qc$total_abundance) != 1){
    if(warning == TRUE){
    warning("Relative abundance does not sum to 1 for all samples.")
    print("Total relative abundance summary:")
    print(summary(qc))
    }
    
    if(warning == FALSE){
      print("Relative abundance does not sum to 1 for all samples.")
      print("Total relative abundance summary:")
      print(summary(qc))
      
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
    
  }
  
  rel <- rel %>% dplyr::select(c(OTU,Sample, Abundance, StudyID)) %>%
    arrange(Sample, desc(Abundance))%>%
    group_by(Sample) %>%
    tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
    replace(is.na(.), 0) %>% ungroup() %>% dplyr::select(-Sample)
  
  return(rel)
}

#Get mean relative abundance for each taxa
rel.mean.extract <- function(rel){
  relmean <- rel %>% select(-StudyID) %>% summarise_all(mean) %>% mutate(stat = "mean") %>% 
    pivot_longer(!stat, names_to = "OTU", values_to = "mean") %>% dplyr::select(-stat)
  
  return(relmean)
}

#Get taxa prevalence 
rel.prev.extract <- function(rel){
  relprev <- rel %>% select(-StudyID) %>% 
    mutate_all(function(x) (ifelse(x>0, 1, 0))) %>% #converts taxa presence as binary (1 = yes, 0 = no)
    summarise_all(mean) %>%
    mutate(stat = "prev") %>% pivot_longer(!stat, names_to = "OTU", values_to = "prev") %>% dplyr::select(-stat)
  
  return(relprev)
}

#Combine mean relative abundance and prevalence into one dataframe
rel.criteria.extract <- function(relmean, relprev){
  rel.criteria <- relmean %>% left_join(relprev, by = "OTU") %>% dplyr::rename(mean_relabund = mean)
  summary(rel.criteria %>% select(-OTU))
  return(rel.criteria)
}


#Mean relative abundance and prevalence data exploration
mra.prev.vis <- function(rel.criteria, #output from rel.criteria.extract
                         export_prefix  #file prefix
)
  #Exports density and scatter plots of prevalence and mean relative abundance to files in the current wd
  {
  density(log10(rel.criteria$mean_relabund)) %>% plot(main = "Log 10 Relative Abundance") #density
  p1 <- recordPlot()
  ggdraw(p1)
  
  density(rel.criteria$prev) %>% plot(main = "Prevalence")
  p2 <- recordPlot()
  ggdraw(p2)
  
  
  
  #Scatter
  p <- ggplot(rel.criteria, aes(x = prev, y = mean_relabund)) + geom_point() + 
    theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
    annotation_logticks(sides = "l")  
  
  #Export
  tiff(paste0(export_prefix, "_mra_prev.tiff"), height = 500, width = 1800, res = 90)
  print(cowplot::plot_grid(p1, p2, p, ncol = 3))
  dev.off()
}





physeq.mra.prevalence <- function(physeq,
                                  type, #one of: "counts" or "relabund" 
                                  vis.export = TRUE, #export mra/prev density function and scatter plots
                                  export_prefix.in #prefix for export file names 
                                  )
  #Summary function that uses other microbiome relative abundance functions 
  #Input phyloseq object of type: "counts" or "relabund"
  #Outputs taxa mean relative abundance and taxa prevalence (rel.criteria [[1]]), rel abund ASV table [[2]], and abundance phyloseq object (dat_rel [[3]])
  #Option to export mean relative abundance and prevalence visualization plots (mra/prev probability density functions/distribution, scatter plots) 
  {
  #Unit test: type has valid input
  if(type != "counts" & type != "relabund"){
    print("Invalid type input.")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
    
  }
  
  #quality control object for determining if relative abundance or count data was input
  qc <- physeq %>% otu_table %>% data.frame %>% pull(1) %>% sum 

  
  if(type == "counts"){
    #Unit test: is it count data? 
    if(qc < 1){
      print("ERROR: first sample looks like relative abundance data. Current input is count data")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
      
    }
    
    #Convert phyloseq counts to relative abundance phyloseq and ASV table
    rel.results <- physeq.convert.rel(physeq = physeq)
    
    #Extract mean relative abundance and prevalence
    rel <- rel.results[[2]]
    relmean <- rel.mean.extract(rel = rel)
    relprev <- rel.prev.extract(rel = rel)
    rel.criteria <- rel.criteria.extract(relmean, relprev)
    
    #dat_rel
    dat_rel <- rel.results[[1]]
  }
  if(type == "relabund"){
    
    #Unit test: is it relative abundance data? 
    if(qc > 1){
      print("ERROR: first sample looks like count data. Current input is relative abundance data")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
    
    #Extract mean relative abundance and prevalence
    rel <- physeqrel.dataframe(physeq)
    relmean <- rel.mean.extract(rel = rel)
    relprev <- rel.prev.extract(rel = rel)
    rel.criteria <- rel.criteria.extract(relmean, relprev)
    
    #dat_rel
    dat_rel <- physeq
  }
  
  
  #Visualize and export current mra and prevalence
  if(vis.export == TRUE){
    mra.prev.vis(rel.criteria = rel.criteria, export_prefix = export_prefix.in)
  }
  
  return(list(rel.criteria, #1 = mean relative abundance and prevalence of each taxa
              rel, #2 = relative abundance ASV table 
              dat_rel #3 = relative abundance phyloseq object
              ))
}

# Function: Glom algorithm -----------------------------------------------
retain.resolve_genus <- function(physeq, #count phyloseq object
                                 dir = getwd(),
                                 relabund_prev.ASV = c(0.001, 0.1), #threshold for mean relative abundance [1] AND prevalence [2]. Taxa must pass both criteria for retention. ASV level filter. 
                                 prev.ASV = 0.5, #threshold for prevalence-only. ASV level filter. 
                                 relabund_prev.genus = c(0.0001, 0.1), #threshold for mean relative abundance [1] AND prevalence [2]. Taxa must pass both criteria for retention. genus level filter. ,
                                 prev.genus = 0.5 #threshold for prevalence-only. Genus level filter. 
                                 )
  
  #If taxa pass relabund_prev or prev filters, they will be retained.
  {
  
  # ====== Unit tests: User input ====== 
  #Unit test: did user input count data? 
  qc <- physeq %>% otu_table %>% data.frame %>% pull(1) %>% sum 
  if(qc < 1){
    print("ERROR: first sample looks like relative abundance data. Only count data is accepted")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
  
  #Unit test: is relative abundance 1 for input physeq? 
    #If total relative abundance is not 1, will kill script
  physeq.convert.rel(physeq = physeq,
                     warning.in = FALSE) 
  
  
  
  
  # ====== Set new directory and create subdirectories ======
  
  #dir = grandparent directory
  #retain_resolve.path-resolve = parent directory
  
  #Create retain-resolve (parent)
  setwd(dir)
  newdir <- "retain-resolve"
  dir.create(newdir)
  retain_resolve.path <- paste0(dir, "/", newdir)
  setwd(retain_resolve.path)
  
  print(retain_resolve.path)
  
  
  #Abundance 
  newdir <- "totalabundance"
  dir.create(newdir)
  abund_diff.path <- paste0(retain_resolve.path, "/", newdir)
  
  #Vis
  newdir <- "plotresults"
  dir.create(newdir)
  vis.path <- paste0(retain_resolve.path, "/", newdir)
  
  #Taxa lists
  newdir <- "taxa_lists"
  dir.create(newdir)
  taxa_lists.path <- paste0(retain_resolve.path, "/", newdir)
  
  
  #mra_prev
  newdir <- "mra_prevalence"
  dir.create(newdir)
  mra_prev.path <- paste0(retain_resolve.path, "/", newdir)
  
  
  # ====== Unfiltered dataset - mra and prevalence ======
  print("Step 1: Create dat_empty. Phyloseq object with absent taxa (mra = 0) removed.")
  dat <- physeq
  
  setwd(vis.path)
  mra.prev.results <- physeq.mra.prevalence(dat, type = "counts", vis.export = TRUE, export_prefix.in = "unfiltered")
  setwd(retain_resolve.path)
  #mra.prev.results[[1]] = mean relative abundance and prevalence of each taxa (rel.criteria)
  #mra.prev.results[[2]] = relative abundance ASV table (rel)
  #mra.prev.results[[3]] = relative abundance phyloseq object (dat_rel)
  
  rel.criteria <- mra.prev.results[[1]]
  rel <-  mra.prev.results[[2]]
  dat_rel_unfiltered <- mra.prev.results[[3]]
  
  
  #Export and save unfiltered data  
  setwd(mra_prev.path)
  write.csv(file = "mraprev_unfiltered.csv", x = rel.criteria, row.names = FALSE)
  setwd(retain_resolve.path)
  
  # === Identify and remove empty taxa === 
  #List of empty taxa (prev = 0)
  emptytaxa <- rel.criteria %>% filter(prev == 0) %>% pull(OTU)
  emptytaxa.df <- rel %>% select(all_of(emptytaxa)) %>% summarise_all(mean) %>% mutate(dummy = "dummy") %>% pivot_longer(!dummy, names_to = "OTU", values_to = "mean_abund") %>% select(-dummy)
  summary(emptytaxa.df %>% select(-OTU))
  
  
  
  #Pruned phyloseq object - empty taxa removed 
  emptytaxa.prune <- rel.criteria %>% filter(!(OTU %in% emptytaxa))
  dat_empty <- prune_taxa(emptytaxa.prune$OTU, dat)
  
  #Log results
  log.taxa <- data.frame("level" = "empty",
                         "pass" = rel.criteria %>% filter(!(OTU %in% emptytaxa)) %>% nrow(), 
                         "fail" = length(emptytaxa))
  
  
  
  #Extract mean relative abundance and prevalence - empty taxa removed 
  setwd(vis.path)
  mra.prev.results <- physeq.mra.prevalence(dat_empty, 
                                            type = "counts", vis.export = TRUE, export_prefix.in = "rm-empty")
  setwd(retain_resolve.path)
  rel.criteria <- mra.prev.results[[1]]
  
  #Save dat_empty_rel
  dat_empty_rel <- mra.prev.results[[3]]
  
  
  
  
  
  
  # ====== Identify and separate ASVs that pass/fail criteria ======
  print("Step 2: Identify ASVs that pass or fail criteria.")
  
  pass <- rel.criteria %>% filter(mean_relabund > relabund_prev.ASV[1] & prev > relabund_prev.ASV[2] | prev > prev.ASV) %>% select(OTU) #Taxa that do not need to be glommed (meet criteria)
  fail <- rel.criteria %>% filter(!OTU %in% pass$OTU & !OTU %in% emptytaxa.df$OTU) %>% select(OTU) #Taxa that require glom (below abund/prevalence thresholds) & with empty taxa removed 
  
  #Log results
  log.taxa.asv <- data.frame("level" = "ASV",
                             "pass" = nrow(pass), 
                             "fail" = nrow(fail))
  log.taxa <- rbind(log.taxa, log.taxa.asv)
  
  #Pass criteria phyloseq
  dat_g1_pass <- prune_taxa(pass$OTU, dat_empty)
  
  #Fail criteria phyloseq
  dat_g1_fail <- prune_taxa(fail$OTU, dat_empty)
  
  #Export pass/failed taxa
  fail.tax.preglom <- as(access(dat_g1_fail, "tax_table"), "matrix")[, 6] %>% as.data.frame() %>% rename(genus = ".")
  pass.tax <- as(access(dat_g1_pass, "tax_table"), "matrix")[, 6] %>% as.data.frame() %>% rename(genus = ".")
  
  setwd(taxa_lists.path)
  write.csv(x = fail.tax.preglom, file = "failASV.csv")
  write.csv(x = pass.tax, file = "passASV.csv")
  setwd(retain_resolve.path)
  
  # ====== Resolve to genus ====== 
  print("Step 3: Agglomerate failed ASVs to genus.")
  
  #Glom
  dat_g1_fail_glom <- tax_glom(dat_g1_fail, NArm = FALSE, taxrank = "Genus")
  
  #Postglom taxa list 
  tax.postglom <- dat_g1_fail_glom %>% tax_table %>% as.data.frame() %>% mutate(taxa_full = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "_")) %>% rownames_to_column("sp")
  
  #Export
  setwd(taxa_lists.path)
  write.csv(tax.postglom, "glomgenus_unfiltered.csv", row.names = FALSE)
  setwd(retain_resolve.path)
  
  #Merge retained ASVs and resolved genus'
  dat_retres_unfiltered <- merge_phyloseq(dat_g1_pass, dat_g1_fail_glom) #relative abundance should sum to 1 for each sample
  
  # ====== Unit tests: Total abundance post-glom ====== 
  #Abundance difference between dat_empty and dat_resolved_unfiltered
  abund_diff <- cbind(dat_empty %>% otu_table %>% data.frame() %>% colSums() %>% data.frame(), 
                      dat_retres_unfiltered %>% otu_table %>% data.frame() %>% colSums() %>% data.frame())
  colnames(abund_diff) <- c("original_abund", "glom_abund")
  abund_diff <- abund_diff %>% mutate(Other = original_abund - glom_abund) %>% mutate(percent_diff = round(Other/original_abund * 100, 2))
  
  print("Abundance difference between original dataframe (dat_empty) and unfiltered retain-resolved: ")
  print(abund_diff)
  print("Note: Other = count difference between total sums")
  
  setwd(abund_diff.path)
  write.csv(x = abund_diff, file = "abund_diff_glomunfiltered.csv")
  setwd(retain_resolve.path)

  #Unit test: is abund_diff 0 for all samples? I.e. no counts went missing during the glom. 
  qc.abund_diff <- abund_diff %>% filter(Other != 0)
  if(nrow(qc.abund_diff) != 0){
    print("Error: at least 1 sample has a different total abundance after agglomerating")
    print("Samples of concern:")
    print(qc.abund.diff)
    
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
    
  }
  
  
  # ====== Filter genus that fail criteria ======
  print("Step 4: Filter genus-level taxa.")
  #Merged ASV and genus relative abundance mra and prevalence 
  setwd(vis.path)
  mra.prev.results <- physeq.mra.prevalence(physeq = dat_retres_unfiltered,
                                            type = "counts", 
                                            vis.export = TRUE, 
                                            export_prefix.in = "glom_unfiltered")
  setwd(retain_resolve.path)
  
  rel.criteria <- mra.prev.results[[1]] %>% #relative abundance is calculated with both ASV and genus-level taxa 
    filter(OTU %in% (dat_g1_fail_glom %>% tax_table %>% data.frame %>% rownames)) #but only include taxa that are genus-level for physeq pruning
  
  #Export rel.criteria
  setwd(mra_prev.path)
  write.csv(x = rel.criteria, file = "mraprev_glomgenus.csv", row.names = FALSE)
  setwd(retain_resolve.path)
  
  #Identify genus-level taxa that pass/fail criteria
  pass <- rel.criteria %>% filter(mean_relabund > relabund_prev.genus[1] & prev > relabund_prev.genus[2] | prev > prev.genus) 
  fail <- rel.criteria %>% filter(!OTU %in% pass$OTU)
  
  
  log.taxa.genus <- data.frame("level" = "genus",
                             "pass" = nrow(pass), 
                             "fail" = nrow(fail))
  log.taxa <- rbind(log.taxa, log.taxa.genus)
  
  
  #Filter genus-agglomerated taxa; pass criteria only
  dat_g2_pass <- prune_taxa(pass$OTU, dat_g1_fail_glom)
  
  #Export log.taxa
  write.csv(x = log.taxa, file = "ntaxa_glom_log.csv", row.names = FALSE)
  
  # ====== Save phyloseq species ID of ASV-level and genus-level taxa ====== 
  print("Step 5: Save whether phyloseq species ID are ASV-level or genus-level taxa")
  genus.taxa <- data.frame("taxa" = dat_g2_pass %>% tax_table %>% rownames, 
                           "level" = rep("genus", nrow(dat_g2_pass %>% tax_table)))
  asv.taxa <- data.frame("taxa" = dat_g1_pass %>% tax_table %>% rownames, 
                         "level" = rep("ASV", nrow(dat_g1_pass %>% tax_table)))
  glom.spID.levels <- rbind(asv.taxa, genus.taxa)
  write.csv(glom.spID.levels, "glom_spID_levels.csv", row.names = FALSE)
  
  
  
  
  # ====== Filtered and merged retain-resolve ASV and genus-level phyloseq object ======
  print("Step 6: Final retained resolve phyloseq object (dat_glom).") 
  # === Merge ASV and genus-level phyloseq objects === 
  dat_glom <- merge_phyloseq(dat_g1_pass, #ASV-level
                             dat_g2_pass #genus-level; filtered
                             )

  # === Identify abundance differences === 
    #Total relative abundance will no longer sum to 1 because some counts are lost from genus filtering
  abund_diff <- cbind(dat_empty %>% otu_table %>% data.frame() %>% colSums() %>% data.frame(), 
                      dat_glom %>% otu_table %>% data.frame() %>% colSums() %>% data.frame())
  colnames(abund_diff) <- c("original_abund", "glom_abund")
  abund_diff <- abund_diff %>% mutate(Other = original_abund - glom_abund) %>% mutate(percent_diff = round(Other/original_abund * 100, 2))
  
  print("Abundance difference between original dataframe (dat_empty) and filtered retain-resolved: ")
  print(abund_diff)
  print("Note: Other = count difference between total sums")
  
  #Export abundance differences
  setwd(abund_diff.path)
  write.csv(x = abund_diff, file = "abund_diff_glomfiltered.csv")
  setwd(retain_resolve.path)
  
  # ====== Other taxa ======
  print("Step 7: Add other taxa for analyses that require accurate total counts - CLR, relabund - (dat_glom_aldex).")
  #Add "Other" taxa label that accounts for filtered counts
    #Required for relative abundance and CLR/aldex
  #mapdf
  mapdf.other <- dat %>% sample_data %>% data.frame
  
  #ASV df
  other_taxa <- abund_diff %>% dplyr::select(Other) %>% t() #other = other taxa
  asvdf.other <- other_taxa %>% as.matrix()
  row.names(asvdf.other) <- "other"
  
  #Taxa df
  taxdf.other <- data.frame(c("Other", "Other", "Other", "Other", "Other", "Other")) %>% t() %>% data.frame()
  row.names(taxdf.other) <- "other"
  colnames(taxdf.other) <-  colnames(physeq %>% tax_table %>% data.frame)
  taxdf.other <- taxdf.other %>% as.matrix()
  
  #Create phyloseq object 
  dat_other <- phyloseq(otu_table(asvdf.other, taxa_are_rows = TRUE), tax_table(taxdf.other), sample_data(mapdf.other))

  
  #Merge phyloseq other
  dat_glom_aldex <- merge_phyloseq(dat_glom, #Filtered ASVs and Genus
                                   dat_other) #Lost taxa
  
  #Unit test: dat_glom_aldex total counts = dat_empty total counts
  abund_diff <- cbind(dat_empty %>% otu_table %>% data.frame() %>% colSums() %>% data.frame(), 
                      dat_glom_aldex %>% otu_table %>% data.frame() %>% colSums() %>% data.frame())
  colnames(abund_diff) <- c("original_abund", "glom_abund")
  abund_diff <- abund_diff %>% mutate(Other = original_abund - glom_abund) %>% mutate(percent_diff = round(Other/original_abund * 100, 2))
  
  print("Abundance difference between original dataframe (dat_empty) and filtered retain-resolved with Other taxa: ")
  print(abund_diff)
  print("Note: Other = count difference between total sums")
  
  #Export abundance differences
  setwd(abund_diff.path)
  write.csv(x = abund_diff, file = "abund_diff_glomfiltered_other.csv")
  setwd(retain_resolve.path)
  
  qc.abund_diff <- abund_diff %>% filter(Other != 0)
  if(nrow(qc.abund_diff) != 0){
    print("Error: at least 1 sample has a different total abundance. Step: dat_glom_aldex")
    print("Samples of concern:")
    print(qc.abund.diff)
    
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
    
  }
  
  
  
  # ====== Relative abundance phyloseq objects ====== 
  print("Step 8: Relative abundance conversions - (dat_empty_rel, dat_rel_glom).")
  
  
  # === dat_rel_glom === 
  setwd(vis.path)
  
  mra.prev.results <- physeq.mra.prevalence(physeq = dat_glom_aldex,
                                            type = "counts",
                                            vis.export = TRUE, 
                                            export_prefix.in = "glom-filtered-other")
  setwd(retain_resolve.path)
  
  dat_rel_glom <- mra.prev.results[[3]]
  
  #Remove other taxa
  rm.other <- dat_rel_glom %>% tax_table %>% data.frame %>% rownames_to_column(var = "OTU") %>% filter(OTU != "other") %>% pull(OTU)
  dat_rel_glom <- prune_taxa(rm.other, dat_rel_glom)
  

  # === dat_empty_rel === 
  convert.rel.results <- physeq.convert.rel(physeq = dat_empty)
  dat_empty_rel <- convert.rel.results[[1]]
  
  
  
  # ====== Export phyloseq objects to file ======
  print("Step 9: Export physeq objects")
  #Set wd
  setwd(dir)
  newdir <- "Physeq_glommed_objects"
  dir.create(newdir)
  setwd(paste0(dir, "/", newdir))
  
  print(paste0("Export to:",  paste0(dir, "/", newdir)))
  #Save
  save(dat_glom, dat_rel_glom, dat_glom_aldex, 
       dat_empty, dat_empty_rel, dat, file = "physeq_ResolveGenus.RData")
  
  
  
}


