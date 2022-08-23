################################################################################
# PARENCHYMA
################################################################################

extract_parenchyma <- function(parenchyma){
  
  # extract data
  max_size <- max(parenchyma$area)
  min_size <- min(parenchyma$area)
  mean_size <- mean(parenchyma$area)
  sd_size <-  sd(parenchyma$area)
  cell_diameter <- sqrt((parenchyma$area)/pi)/1000  # mean radius = mean(sqrt(A)/pi) and not the formula with mean area
  mean_diameter <- mean(cell_diameter) 
  
  
  max_size_df <- data.frame(name = "stele", type = "max_size", value = max_size)
  min_size_df <- data.frame(name = "stele", type = "min_size", value = min_size)
  mean_size_df <- data.frame(name = "stele", type = "mean_size", value = mean_size)
  mean_diameter_df <- data.frame(name = "stele", type = "cell_diameter", value = mean_diameter)
  sd_size_df <- data.frame(name = "stele", type = "SD", value = sd_size)
  
  # store in a dataframe
  data_parenchyma <- rbind(mean_diameter_df, 
                           sd_size_df,
                           mean_size_df)
  
  # return the df
  return(data_parenchyma)
}


################################################################################
# XYLEM
################################################################################

extract_xylem <- function(xylem){
  
  # extract data
  n_files <- length(xylem$xylem)
  max_size <- max(xylem$area)
  min_size <- min(xylem$area)
  mean_size <- mean(xylem$area)
  order <- 1.5
  
  n_files_df <- data.frame(name = "xylem", type = "n_files", value = n_files)
  max_size_df <- data.frame(name = "xylem", type = "max_size", value = max_size)
  min_size_df <- data.frame(name = "xylem", type = "min_size", value = min_size)
  mean_size_df <- data.frame(name = "xylem", type = "mean_size", value = mean_size)
  order_df <- data.frame(name = "xylem", type = "order", value = order)
  
  # store in a dataframe
  data_xylem <- rbind(n_files_df, 
                      max_size_df, 
                      mean_size_df,
                      order_df)
  
  # return the df
  return(data_xylem)
}


################################################################################
# PHLOEM
################################################################################


extract_phloem <- function(phloem){
  
  # extract data
  
  # Size as area
  #max_size <- max(phloem$area)
  #min_size <- min(phloem$area)
  #mean_size <- mean(phloem$area)
  
  # Size as diameter
  max_size <- max(phloem$length)/1000 # [mm]
  
  max_size_df <- data.frame(name = "phloem", type = "max_size", value = max_size)
  # min_size_df <- data.frame(name = "phloem", type = "min_size", value = min_size)
  # mean_size_df <- data.frame(name = "phloem", type = "mean_size", value = mean_size)
  
  
  # store in a dataframe
  data_phloem <- rbind(max_size_df)
  
  # return the df
  return(data_phloem)
}

################################################################################
# STELE
################################################################################

extract_stele <- function(stele){
  
  # extract data
  totarea <- stele$area
  order <- 1.0
  
  totarea_df <- data.frame(name = "stele",
                     type = "totarea",
                     value = totarea)
  
  order_df <- data.frame(name = "stele", type = "order", value = order) 
  
  # store in a dataframe
  data_stele <- rbind(totarea_df, 
                      order_df)
  
  # return the df
  return(data_stele)
}


################################################################################
# cortex
################################################################################

extract_cortex <- function(cortex){
  
  cell_diameter <- mean(cortex$length)/1000
  order <- 4.0                         # by default
  
  data_cortex = data.frame(name = c("cortex", "cortex"),
                              type = c("cell_diameter", "order"),
                              value = c(cell_diameter, order))
  
  return(data_cortex)
}

################################################################################
# PERICYCLE
################################################################################

extract_pericycle <- function(pericycle){
  
  cell_diameter <- mean(pericycle$length)/1000
  order <- 2.0
  n_layers <- 1.0
  
  cell_diameter_df <- data.frame(name = "pericycle", type = "cell_diameter", value = cell_diameter)
  order_df <- data.frame(name = "pericycle", type = "order", value = order)
  n_layers_df <- data.frame(name = "pericycle", type = "n_layers", value = n_layers)
  
  data_pericycle = rbind(cell_diameter_df,
                         order_df,
                         n_layers_df)
  
  return(data_pericycle)
}

################################################################################
# ENDODERMIS
################################################################################

extract_endodermis <- function(endodermis){
  
  cell_diameter <- mean(endodermis$length)/1000
  n_layers <- 1.0                    # by default
  order <- 3.0                         # by default
  
  data_endodermis = data.frame(name = c("endodermis", "endodermis", "endodermis"),
                              type = c("cell_diameter", "n_layers", "order"),
                              value = c(cell_diameter, n_layers, order))
  
  return(data_endodermis)
}

################################################################################
# EPIDERMIS
################################################################################

extract_epidermis <- function(epidermis){
  
  cell_diameter <- mean(epidermis$length)/1000
  n_layers <- 1.0                    # by default
  order <- 6.0                         # by default
  
  data_epidermis = data.frame(name = c("epidermis", "epidermis", "epidermis"),
                              type = c("cell_diameter", "n_layers", "order"),
                              value = c(cell_diameter, n_layers, order))
  
  
  return(data_epidermis)
}


################################################################################
# EXODERMIS
################################################################################

extract_exodermis <- function(exodermis){
  
  cell_diameter <- mean(exodermis$length)/1000
  n_layers <- 1.0                    # by default
  order <- 5.0                         # by default
  
  data_exodermis = data.frame(name = c("exodermis", "exodermis", "exodermis"),
                              type = c("cell_diameter", "n_layers", "order"),
                              value = c(cell_diameter, n_layers, order))
  
  return(data_exodermis)
}


################################################################################
# Extract BARRIERS
################################################################################

extract_barriers <- function(barriers){
  
  b_endo <- data_frame(name = "endodermis", type = "barrier", value = barriers$stage[barriers$type == "endodermis"])
  b_exo <- data_frame(name = "exodermis", type = "barrier", value = barriers$stage[barriers$type == "exodermis"])
  
  data_barriers <- rbind(b_endo, b_exo)
  
  return(data_barriers)
  
}


################################################################################
# Extract LAYERS
################################################################################

extract_layers <- function(layers){
  
  l_stele <- layers$layers[layers$name == "stele"]
  d_stele <- 2*layers$radius[layers$name == "stele"]/1000 # µm -> mm
  l_cortex <- layers$layers[layers$name == "cortex"]
  r_cortex <- layers$radius[layers$name == "cortex"]  
  phloem_layer <- layers$layers[layers$name == "phloem"]
  phloem_radius <- layers$radius[layers$name == "phloem"]/1000 # µm -> mm
  
  # proportion of phloem in the stele
  phloem_prop <- phloem_radius/(d_stele/2)
  
  l_stele_df <- data.frame(name = "stele", type = "n_layers", value = l_stele)
  d_stele_df <- data.frame(name = "stele", type = "layer_diameter", value = d_stele)
  l_cortex_df <- data.frame(name = "cortex", type = "n_layers", value = l_cortex)
  #r_cortex_df <- data.frame(name = "cortex", type = "radius", value = r_cortex)
  phloem_layer_df <- data.frame(name = "phloem", type = "n_files", value = phloem_layer)
  phloem_prop_df <- data.frame(name = "phloem", type = "proportion", value = phloem_prop)
  
  data_layers <- rbind(l_stele_df, 
                       d_stele_df, 
                       l_cortex_df,
                       phloem_layer_df,
                       phloem_prop_df)
  
  return(data_layers)
}


################################################################################
# Extract all data
################################################################################

extract_all <- function(path){
  # input : a path of xlsx file with quantification (ex : A01_S1_C01_Z04_01.xlsx)
  # output : a dataframe with the parameters of GRANAR
  
  # Basic parameters : 
  b_params <- data.frame(name = c("secondarygrowth", "randomness", "planttype"),
                         type = c("param", "param", "param"),
                         value = c(1, 3, 2))
  
  # aerenchyma
  data_aerenchyma <- data.frame(name = c("aerenchyma", "aerenchyma"),
                                type = c("proportion", "n_files"),
                                value = c(0.0, 10.0))
  
  
  # quantified tissues
  xylem <- read_xlsx(path, sheet = "Xylem")
  parenchyma <- read_xlsx(path, sheet = "Parenchyma")
  stele <- read_xlsx(path, sheet = "Stele")
  cortex <- read_xlsx(path, sheet = "Cortex")
  pericycle <- read_xlsx(path, sheet = "Pericycle")
  endodermis <- read_xlsx(path, sheet = "Endodermis")
  epidermis <- read_xlsx(path, sheet = "Epidermis")
  phloem <- read_xlsx(path, sheet = "Phloem")
  barriers <- read_xlsx(path, sheet = "Barriers")
  layers <- read_xlsx(path, sheet = "Layers")
  exodermis <- read_xlsx(path, sheet = "Exodermis")
  
  # Extract information from subvariables
  data_xylem <- extract_xylem(xylem)
  data_parenchyma <- extract_parenchyma(parenchyma)
  data_phloem <- extract_phloem(phloem)
  data_stele <- extract_stele(stele)
  data_cortex <- extract_cortex(cortex)
  data_endodermis <- extract_endodermis(endodermis)
  data_epidermis <- extract_epidermis(epidermis)
  data_exodermis <- extract_exodermis(exodermis)
  data_pericycle <- extract_pericycle(pericycle)
  data_barriers <- extract_barriers(barriers)
  data_layers <- extract_layers(layers)

  # Extract each variables
  data_all <- rbind(b_params, 
                    data_aerenchyma,
                    data_xylem,
                    data_parenchyma,
                    data_phloem,
                    data_stele,
                    data_cortex,
                    data_endodermis,
                    data_exodermis,
                    data_epidermis,
                    data_pericycle,
                    data_barriers,
                    data_layers)  

  return(data_all)  
}


################################################################################
# Extract identifyiers
################################################################################

extract_ids <- function(list, i){
  #_____________________________________________________________________________
  # @list = list of all xlsx files
  # @i    = iteration
  # this function will extract metadata from the name of the xlsx file 
  # and generate a dataframe row :
  # id | plant_id | segment | cut | zoom | repetition
  #_____________________________________________________________________________
  
  # take the i string
  temp_string <- list[i]
  
  # split to list
  temp_string2 <- strsplit(temp_string, "_")
  
  # define params
  temp_id <- i
  temp_plant_id   <- temp_string2[[1]][1]
  temp_segment    <- temp_string2[[1]][2]
  temp_cut        <- temp_string2[[1]][3]
  temp_zoom       <- temp_string2[[1]][4]
  temp_repetition <- strsplit(temp_string2[[1]][5], ".xlsx")[[1]][1]
  
  # create dataframe
  temp_df = data.frame(id = temp_id,
                       plant_id = temp_plant_id,
                       segment = temp_segment,
                       cut = temp_cut,
                       zoom = temp_zoom,
                       repetition = temp_repetition)
  
  return(temp_df)
}


################################################################################
# Nice boxplot
################################################################################

nice_boxplot1 <- function(data, X, Y, Z){
  # @data = a dataframe
  # @X = variable x        (ex : data$plant_id)
  # @Y = variable y        (ex : data$length)
  # @Z = variable to group (ex : data$genotype)
  
  g <- data %>% ggplot(aes(x = X, y = Y)) +
    geom_point() +
    geom_boxplot() +
    
    #facet_wrap(~ Z) +
    
    theme_bw()   
  return(g)
  
}


################################################################################
# Nice boxplot 2
################################################################################

nice_boxplot2 <- function(data = data, X = x, Y = y, Z = z){
  
  
      g <- ggplot(data, aes(x=X, y=Y, fill=Z)) +
        
        geom_boxplot() +
        
        scale_fill_viridis(discrete = TRUE, alpha=0.6) + # Colors
        geom_jitter(color="black", size=0.4, alpha=0.9) + # Add points
        
        # Theme params
        theme(
          legend.position="none",
          plot.title = element_text(size=11)
        ) +
        
        theme_bw() +
        theme(axis.line = element_line(color='black'),
              
              # Remove grid
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              
              # Change font
              text = element_text(family = "A"),
              
              # Remove x ticks
              #axis.text.x=element_blank(),
              #axis.ticks.x=element_blank(),
              
              # Remove legend title
              #legend.title=element_blank()
              
              # Remove legend
              legend.position="none"
              
              # Remove panel borders 
              #panel.border = element_blank()
              
          
        )
        
      return(g)}


