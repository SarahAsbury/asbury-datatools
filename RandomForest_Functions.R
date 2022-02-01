#Functions to run random forest and process results



# Load Packages -----------------------------------------------------------
library(tidyverse)
library(ggpubr)



# Functions: Split Methods ------------------------------------------------



# Function: Regression RF ----------------------------------------------------



# Function: Classification RF ---------------------------------------------


# Function:  -------------------------------------




# Functions: Variable importance formatting --------------------------------
varimport_import <- function(df)
  #edit varimp df imported from CSV 
  {
  df.out <- df %>% rename(predictors = X)
}

select.varimp <- function(rf.type, metric)
  #Select variable importance metric (gini, mse, mda)
  {
  #Unit test: rf.type and metric inputs compatible.
  try(if((rf.type == "class" & metric == "mse") | (rf.type == "reg" & metric == "mda"))
    stop("rf.type and metric incompatabile. Ensure you are using the correct variable importance metric for your random forest type.")
  )
  
  
  
  if(rf.type == "reg"){
    if(metric == "mse"){
      y <- "X.IncMSE_name"
      ylab <- "Increased Mean Square Error"
    }
    
    if(metric == "gini"){
      y <- "IncNodePurity_name"
      ylab <- "Increased node purity"
    }
  }
  
  if(rf.type == "class"){
    if(metric == "mda"){
      y <- "MeanDecreaseAccuracy"
      ylab <- "Decreased Accuracy"
    }
    
    if(metric == "gini"){
      y <- "MeanDecreaseGini"
      ylab <- "Increased node purity"
    }
    
  }
  return(c(y, ylab))
}


# Function: Variable Importance -------------------------------------------



varimport_plot <- function(varimp, #dataframe of variable importance
                           rf.type, #One of: "class" or "reg"
                           metric = "gini", #One of: "mse", "gini", or "mda)
                           selection_type = "random_top", #One of: random_top, top
                           top = 20, #number of "top" predictor variables to plot
                           xlab = NA, #Labels to use form plotting. NA will use R column names. Accepts a character vector or "extract" instructions.
                           extract.names.df = NA  
)
  #Lollipop plot of variable importance 

  #Varimp df structure: 
    #Predictor variable column must be named ("predictors")
    #Use varimp import 
  
  
  #Selection_type: 
    #top = top 20 (or specified) predictor variable (PV)
    #random_top = #Take top 20 (or specified) predictor variables (PV), bottom 5 (specified*0.25) PV, and 15 (specified*0.75) additional random PV for visualization

{
  #=== X and y variables ===
  select <- select.varimp(rf.type, metric) 
  y <- select[1]
  ylab <- select[2]
  
  print(y) #Check y 
  print(ylab) #Check y label
  
  x <- "predictors" #do not change 

  
  
  
  #=== Variable importance dataframe ===
  if(selection_type == "random_top")
  {
    toppred <- varimp %>% slice_max(order_by = get(y), n = 20)
    minpred <- varimp %>% slice_min(order_by = get(y), n = round(top*0.25, digits = 0))
    randompred <- varimp %>% filter(!(predictors %in% toppred$predictors | predictors %in% minpred$predictors)) %>% 
      slice_sample(n = round(top*0.75, digits = 0))
    
    varimp <- varimp %>% filter(predictors %in% toppred$predictors | predictors %in% minpred$predictors | predictors %in% randompred$predictors)
  }
  
  if(selection_type == "top")
    #Take top 20 taxa (or specified)
  {
    varimp <- varimp %>% slice_max(order_by = get(y), n = top)
  }
  
  varimp <- varimp %>% arrange(desc(y)) #re-order variables
  
  
  
  
  
  
  #=== Predictor variable labels ===
  xlab.edit <- ifelse(class(xlab) == "logical" & xlab[1] %>% is.na, 
                           FALSE, 
                           TRUE) #TRUE = edit predictor labels; FALSE = predictor labels are same as R column names 
  
  if(!(is.na(xlab)) & xlab == "extract")
    #extract method
  {
    #Unit test: extract.names.df provided
    if(class(extract.names.df) == "logical")
    {
      if(extract.names.df[1] %>% is.na){
        print("Error: Provide dataframe to extract variables.")
        print("Running UNDECLARED to kill script.")
        UNDECLARED()
      }
    }
    
    #Align extract.names.df to varimp order
    xlab.extract <- varimp %>% rename(column_names = "predictors") %>% 
      left_join(extract.names.df, by = "column_names")
    
    #Unit test: Top_pred and extract.names.df labels align
    qc <- (xlab.extract$column_names == varimp$predictors) %>% data.frame() %>% rename(qc = ".") %>% count(qc) %>% unique
    print(qc)
    if(qc[1] == FALSE){
      print("Error: Extract names dataframe labels do not align with predictor variable order")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
    
    
    xlab <- xlab.extract 
  }
  
  #Unit test: are there enough xlabels supplied?
  if(class(xlab) != "logical"){
    if(nrow(xlab) != nrow(varimp))
      {
    print("Error: Mismatched number of x labels and predictor variables")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
    }
  }

  
  
  
  
  #=== Plot ===
  ggdotchart(varimp, x = paste(x), y = paste(y),
             sorting = "ascending",                        
             ggtheme = theme_pubr(), 
             add = "segment",
             ylab = ylab, 
             xlab = "") + 
    {if(xlab.edit == TRUE)
      scale_x_discrete("", breaks = xlab$column_names, labels = xlab$xlabel_names)
      } + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
}




# Functions: Density plot loop  --------------------------------------------
#Density plot function:
gdens <- function(df, #name of dataframe
                  x, #region to plot 
                  dvar, #group by (e.g genotype)
                  tsize = 24
                  )
{ 
  ggplot(df, aes_string(x = x, fill = dvar)) + 
    geom_density(alpha=0.5) +
    theme_minimal_hgrid(12) +
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(size=tsize)) 
}





multi.density.plot <- function(df, #dataframe containing predictor variables and response variable by subject
                               response, #response variable (text)
                               tsize = 24,
                               top = 20, #number of "top" predictor variables to plot
                               varimp, #variable importance dataframe imported using 
                               rf.type, #One of: "class" or "reg" 
                               metric, #One of: "mse", "gini", or "mda
                               scale = 0.8, #cowplot scale graphs
                               ncol = 2, #cowplot number of cols
                               xlab = NA, #one of: NA, "extract", "blank", or list of custom x labels listed according to plot order 
                               extract.names.df = NA #provide df if xlab = "extract". 1st col must be desired x label names. 2nd col must be how they appear as R columns.
                               )

#extract.names.df formatting:
  #1st column = xlabel_names 
    #X labels to appear in density plot
  #2nd column = column_names
    #How the variable appears in R columns

  {
  #Set variable importance metric 
  select <- select.varimp(rf.type, metric) 
  metric.name <- select[1]

  top_pred <- varimp %>% slice_max(order_by = get(metric.name), n = top) %>% pull(predictors)
  
  #Set x labels
  if(xlab == "blank" & !(is.na(xlab))){
    xlab <- rep("", length(top_pred))
  }
  
  if(xlab == "extract" & !(is.na(xlab)))
    #extract method
    {
    #Unit test: extract.names.df provided
    if(class(extract.names.df) == "logical")
      {
      if(extract.names.df[1] %>% is.na){
        print("Error: Provide dataframe to extract variables.")
        print("Running UNDECLARED to kill script.")
        UNDECLARED()
      }
    }
    
    #Align extract.names.df to top_pred order
    xlab.extract <- top_pred %>% data.frame %>% rename(column_names = ".") %>% 
      left_join(extract.names.df, by = "column_names")

    #Unit test: Top_pred and extract.names.df labels align
    qc <- (xlab.extract$column_names == top_pred) %>% data.frame() %>% rename(qc = ".") %>% count(qc) %>% unique
    if(qc[1] == FALSE){
      print("Error: Extract names dataframe labels do not align with predictor variable order")
      print("Running UNDECLARED to kill script.")
      UNDECLARED()
    }
      
    
    #Extract labels
    xlab <- xlab.extract %>% pull(xlabel_names)
    
  }
  
  #Unit test: are there enough xlabels supplied?
  if(class(xlab) != "logical"){
    if(length(xlab) != length(top_pred)){
    print("Error: Mismatched number of x labels and predictor variables")
    print("Running UNDECLARED to kill script.")
    UNDECLARED()
  }
}
  
  #Multiple plotting:
  dens.df <- data.frame()
  plotlist <- list()
  for (i in 1:length(top_pred))
    {
    j <- top_pred[i]
    print(i)
    print(j)

    p <- gdens(df = df, 
               x = j, 
               dvar = response, 
               tsize = tsize) + theme(legend.position = "none") + 
      {if(class(xlab) != "logical") 
                 xlab(xlab[i])} 
    
    
    assign(paste("dens", i, sep ="."), p)
    plotname <- paste("dens", i, sep =".")
    dens.df <- rbind(dens.df, plotname)
    plotlist[[i]] <- p
  }
  
  #Extract legend
  plotlist[[length(plotlist) + 1]] <- gdens(df = df, x = top_pred[1], dvar = response, tsize = tsize) %>% get_legend()
  
  p_grid <- cowplot::plot_grid(plotlist = plotlist, scale = scale, ncol = ncol)
  return(p_grid)
}

