#Microbiome Tools 

# Load Packages -----------------------------------------------------------



# Function: Create phyloseq object ----------------------------------------
physeq.create <- function(asvdf, #Dataframe; Rownames = sequences (e.g CCTGTT...) . Columns = sequencing counts.
                          taxdf, #Dataframe; Rownames = sequences (e.g CCTGTT...). Columns = Sequence taxonomic identification (Kingdom, Phylum, Class, Order, Family Genus)
                          mapdf, #Dataframe; No specific rownames required. Rows = Samples, Columns = Sample metadata (e.g depression scores, sex, age, etc.)
                          ID.col, #name of the column containing study IDs 
                          export.asvsp = TRUE,
                          export.wd  = get_wd()
                          )
  {
  
  #=== Wrangling === 
  #Standardize ID column 
  mapdf <- mapdf %>% rename(StudyID = paste(ID.col))
  
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



# Function: relative abundance   -------------------------------------------

#phyloseq counts -> phyloseq relative abundance and wide ASV dataframe
physeq.convert.rel <- function(physeq
)
  {
  dat_rel <- transform_sample_counts(physeq,function(x) x/sum(x))
  
  rel <- psmelt(dat_rel)
  
  #Unit test: 
  #Total relative abundance sums to 1
  #Returns warning. Script will continue even if test fails. 
  qc <- rel %>% group_by(Sample) %>% dplyr::summarize(total_abundance = sum(Abundance))
  if(min(qc$total_abundance) != 1 | max(qc$total_abundance) != 1){
    warning("Relative abundance does not sum to 1 for all samples.")
    print("Total relative abundance summary:")
    print(summary(qc))
  }
  

  
  #Convert to wide ASV table
  rel <- rel %>% dplyr::select(c(OTU,Sample, Abundance, StudyID)) %>%
    arrange(Sample, desc(Abundance))%>%
    group_by(Sample) %>%
    tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
    replace(is.na(.), 0) %>% ungroup() %>% select(-Sample)
  

  
  invisible(list(dat_rel, rel))
}


#Get mean relative abundance for each taxa
rel.mean.extract <- function(rel){
  relmean <- rel %>% select(-StudyID) %>% summarise_all(mean) %>% mutate(stat = "mean") %>% 
    pivot_longer(!stat, names_to = "OTU", values_to = "mean") %>% select(-stat)
  
  return(relmean)
}

#Get taxa prevalence 
rel.prop.extract <- function(rel){
  relprop <- rel %>% select(-StudyID) %>% 
    mutate_all(function(x) (ifelse(x>0, 1, 0))) %>% #converts taxa presence as binary (1 = yes, 0 = no)
    summarise_all(mean) %>%
    mutate(stat = "prop") %>% pivot_longer(!stat, names_to = "OTU", values_to = "prop") %>% select(-stat)
  
  return(relprop)
}


#Combine mean relative abundance and prevalence into one dataframe
rel.criteria.extract <- function(relmean, relprop){
  rel.criteria <- relmean %>% left_join(relprop, by = "OTU") %>% dplyr::rename(mean_relabund = mean)
  summary(rel.criteria %>% select(-OTU))
  return(rel.criteria)
}

#Mean relative abundance and prevalence data exploration
mra.prev.vis <- function(rel.criteria, #output from rel.criteria.extract
                         export_prefix  #file prefix
)
  #Exports density and scatter plots of proportion and mean relative abundance to files in the current wd
  {
  density(log10(rel.criteria$mean_relabund)) %>% plot(main = "Log 10 Relative Abundance") #density
  p1 <- recordPlot()
  ggdraw(p1)
  
  density(rel.criteria$prop) %>% plot(main = "Proportion")
  p2 <- recordPlot()
  ggdraw(p2)
  
  
  
  #Scatter
  p <- ggplot(rel.criteria, aes(x = prop, y = mean_relabund)) + geom_point() + 
    theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
    annotation_logticks(sides = "l")  
  
  #Export
  tiff(paste0(pexport_prefix, "mra_prop.tiff", height = 500, width = 1800, res = 90))
  cowplot::plot_grid(p1, p2, p)
  dev.off()
}

# Function: Glom algorithm -----------------------------------------------


retain.resolve_genus <- function(physeq, 
                                 vis.export = TRUE){
  
  # === Identify and remove empty taxa === 
  
  #Convert to relative abundance to identify taxa where mean relative abundance = 0 
  rel.results <- physeq.convert.rel(physeq = physeq)

  
  #Extract mean relative abundance and prevalence
  relmean <- rel.mean.extract(rel.results[[2]])
  relprop <- rel.prop.extract(rel.results[[2]])
  rel.criteria <- rel.criteria.extract(relmean, relprop)
  
  
  #Visualize and export current mra and prevalence
  if(vis.export == TRUE){
  mra.prev.vis(rel.critera, export_prefix = "unfiltered")
  }
  
  
  
  
}


