

library(xml2)
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(viridis)
library(cowplot)
library(retistruct)
library(Hmisc)
library(devtools)
library(packcircles)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(sp)
`%!in%` <- compose(`!`, `%in%`)

create_anatomy <- function(path = NULL,  # path to xml file
                           parameters = NULL,
                           verbatim = F,
                           maturity_x = F,
                           paraview = T){
  # Return NULL is no parameters are specified
  if( is.null(path) & is.null(parameters)){
    warning("Please specify a parameter set for the simulation")
    return(NULL)
  }
  if(!is.null(path)){
    params <- read_param_xml(path)
  }
  if(!(is.null(parameters))){
    params <- parameters
    if(nrow(params[params$name == "planttype",]) == 0){
      warning(paste0("Could not find the 'planttype' information in the parameter input"))
      return(NULL)
    }
    # Quality control
    to_find <- c("planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex")
    for(tf in to_find){
      if (nrow(params[params$name == tf,]) == 0){
        warning(paste0("Could not find the '",tf,"' information in the parameter input"))
      }
    }
    cols_to_find <- c("name", "type", "value")
    for(ctf in cols_to_find){
      if (is.null(params[[ctf]])){
        warning(paste0("Could not find the '",ctf,"' column in the parameter input"))
        return(NULL)
      }
    }
    
  }
  
  # set initial time
  t_1 <- Sys.time()
  t1 <- proc.time()
  # set the random factor
  random_fact <- params$value[params$name == "randomness"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  
  data_list <- cell_layer(params)
  layers <- data_list$layers # layers: cell_type, diameter, n_layer, order
  all_layers <- data_list$all_layers # expand layers
  center <- max(all_layers$radius) # center of the cross section
  xylarea<-params$value[params$name == "xylem" & params$type == "mean_size"]*params$value[params$name == "xylem" & params$type == "n_files"]
  insidestele_diameter<-params$value[params$name == "stele" & params$type == "layer_diameter"]-(2*params$value[params$name == "pericycle" & params$type == "cell_diameter"])
  nparenchyma<-(((((insidestele_diameter/2)*1000)^2)*pi)-xylarea)/(params$value[params$name == "stele" & params$type == "mean_size"])
if(nparenchyma>=1000){
  nparenchyma=nparenchyma*0.86}
if(nparenchyma>=300 && nparenchyma<1000){
    nparenchyma=nparenchyma*0.95}  
  schoppach<-  ((((params$value[params$name == "stele" & params$type == "layer_diameter"]/2)*1000)^2)*pi)
  schoppach2<-params$value[params$name == "stele" & params$type == "totarea"]
  
  schoppach3<-(params$value[params$name == "stele" & params$type == "mean_size"]*nparenchyma)+xylarea
  #nparenchyma2<-((params$value[params$name == "stele" & params$type == "totarea"]-xylarea)/params$value[params$name == "stele" & params$type == "mean_size"])*0.81
  #areaparenchyma<-params$value[params$name == "stele" & params$type == "totarea"]
  #params$value[params$name == "stele" & params$type == "layer_diameter"] <- (sqrt(areaparenchyma/pi)/1000)*2

  
  all_cells <- create_cells(all_layers, random_fact)
  summary_cells <- plyr::ddply(all_cells, .(type), summarise, n_cells = length(angle))
  all_cells$type[substr(all_cells$type, 1,6) == "cortex"] <- "cortex"
  
  all_cells$id_group <- 0
  
  all_cells%>%
    ggplot()+
    geom_point(aes(x,y, colour = type))+
    coord_fixed()
  
  
  if(params$value[params$name == "secondarygrowth"] == 2){all_cells <- vascular(all_cells, params, layers, center)}
  if (params$value[params$name == "secondarygrowth"] == 1){
    packing<-pack_xylem (params,nparenchyma, center)
    rm_stele <- all_cells%>%
      filter(type != "stele")
    new_cells <- rbind(packing, rm_stele)
    new_cells$id_cell <- 1:nrow(new_cells)
  }
  
  new_cells%>%
    ggplot()+
    geom_point(aes(x,y, colour = type))+
    coord_fixed()
  
  all_cells <- new_cells 
  
  if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
    all_cells <- make_pith(all_cells, params, center)
  }
  
  # Get the voronio data
  vtess <- deldir(all_cells$x, all_cells$y, digits = 8)
  if(is.null(vtess)){return(NULL)}
  vorono_list <- cell_voro(all_cells, vtess, center)
  all_cells <- vorono_list$all_cells
  rs2 <- vorono_list$rs2
  
  rs1 <- rs2 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan)
  
  rs1$id_point <- paste0(rs1$x,";",rs1$y)
  
  rs1 <- rs1%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    ungroup()
  
  rs1%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
    coord_fixed()+
    theme_classic()
  
  rs1 <- smoothy_cells(rs1)
  
  # message(paste0("a few possible mistake are possible around point", voiz$id_point[voiz$n < 2]))
  
  if(verbatim) message("Merging inter cellular space")
  rs1 <- fuzze_inter(rs1)
  # rs1 <- rs1[rs1$type != "inter_cellular_space", ]
  
  ini_cortex_area <- sum(all_cells$area[all_cells$type %in% c( "cortex" ,"exodermis" , # ,"endodermis"
                                                               "epidermis", "inter_cellular_space")])
  saved_rs1 <- rs1
  rs1 <- saved_rs1
  if(proportion_aerenchyma > 0){
    
    for (i in unique(rs1$id_cell)) {
      tmp <- rs1[rs1$id_cell == i,]
      tmp <- tmp[!is.na(tmp$x), ]
      pol <- Polygon(tmp[, c("x","y")])
      rs1$area[rs1$id_cell == i] <-  pol@area
    }
    # make aerenchyma
    rs1 <- aerenchyma(params, rs1)
    
    # simplify septa
    rs1 <- septa(rs1)
    # id_aerenchyma <- unique(septum$id_cell)
  }else{cortex_area <- ini_cortex_area}
  
  
  # hairy epidermis # add-on 27-02-2020
  #-----------------------------------------
  
  if(length(params$value[params$name == "hair"] )!= 0){
    if(params$value[params$name == "hair" & params$type == "n_files"] > 0 ){
      rs1 <- root_hair(rs1, params, center)
    }
  }
  
  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type, value = params$value)
  
  if(length(which(is.na(rs1$x)))>0){
    print("NA in cell coordinate ... ")
  }
  rs1 <- rs1[!is.na(rs1$x), ]
  for (i in unique(rs1$id_cell)) {
    tmp <- rs1[rs1$id_cell == i,]
    if(nrow(tmp)> 0){
      pol <- Polygon(tmp[, c("x","y")])
      rs1$area[rs1$id_cell == i] <-  pol@area
      if(pol@area == 0){
        print("cell_area = 0")
        print("this id_cell will be removed")
        rs1 <- rs1[rs1$id_cell != i,]
      }
    }
  }
  
  # Reset the ids of the cells to be continuous
  ids <- data.frame(id_cell = unique(rs1$id_cell))
  ids$new <- c(1:nrow(ids))
  rs1 <- merge(rs1, ids, by="id_cell")
  rs1$id_cell <- rs1$new
  
  if(proportion_aerenchyma > 0){
    id_aerenchyma <- unique(rs1$id_cell[rs1$aer == "aer"])
  }else{id_aerenchyma <- NA}
  
  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new
  
  mX <- mean(rs1$area[rs1$type == "xylem"])
  if(params$value[params$name == "planttype"] == 1){
    rs1$type[rs1$type == "xylem" & rs1$area > mX] <- "metaxylem"
  }
  one_cells <- rs1%>%
    filter(!duplicated(id_cell))# , !duplicated(type), !duplicated(id_group), !duplicated(area)
  
  
  
  all_cells <- merge(all_cells, one_cells, by = "id_cell")
  
  
  # adding the outputs by cell layers
  out <- ddply(all_cells, .(type.y), summarise, n_cells=length(type.y),
               layer_area = sum(area.y),
               cell_area = mean(area.y)) %>%
    mutate(name = type.y) %>%
    dplyr::select(-type.y) %>%
    gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
    mutate(io = "output")%>%
    dplyr::select(io, everything())
  output <- rbind(output, out)
  
  
  time <- as.numeric(Sys.time()-t_1)
  
  # finaly we add the outputs for the whole section
  output <-rbind(output, data.frame(io="output", name="all", type="n_cells", value = nrow(all_cells)))
  
  output <-rbind(output, data.frame(io="output", name="stelar", type="layer_area",
                                    value = sum(all_cells$area.y[all_cells$order < 4])))
  TCA <- sum(all_cells$area.y[all_cells$order > 3])
  output <-rbind(output, data.frame(io="output", name="cortex_alive_to_epidermis", type="layer_area",
                                    value = TCA))
  
  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area.y)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="layer_area", value = (ini_cortex_area - cortex_area)))
  # output <-rbind(output, data.frame(io="output", name="aerenchyma", type="proportion", value = (ini_cortex_area - cortex_area)/ini_cortex_area))
  output <-rbind(output, data.frame(io="output", name="simulation", type="time", value = time))
  
  
  
  rs1$sorting <- c(1:nrow(rs1))
  
  nodes <- vertex(rs1)
  nodes <- nodes[!is.na(nodes$x), ]
  
  # In the MECHA python script, to have unmature metaxylem vessels
  # Metaxylem elements are turned into stele cell type
  if(maturity_x){
    tmp_m <- mean(nodes$area[nodes$type == "cortex"])
    nodes$type[nodes$type == "metaxylem" & nodes$area > tmp_m ] <- "stele"
  }
  
  # comment
  if(paraview){
    walls <- pv_ready(rs1)
    wall_length <- walls%>%select(-x, -y, -xx, -yy)%>% # ends_with(as.character(c(0:9)))
      select(starts_with("x"), starts_with("y"))%>%
      colnames()
    print(wall_length)
    wally <- walls[!duplicated(walls[,wall_length]),] %>%
      dplyr::select(wall_length)
    wally$id_wall <- c(1:nrow(wally))
    walls <- merge(walls, wally, by= wall_length)
    walls <- walls %>%
      # filter(!duplicated(id_wall))%>% # 30/09/2020
      arrange(sorting)
    
  }else{
    wally <- nodes[!duplicated(nodes[,c('x1', 'x2', 'y1', 'y2')]),] %>%
      dplyr::select(c(x1, x2, y1, y2))
    
    wally$id_wall <- c(1:nrow(wally))
    walls <- wally
    
    nodes <- merge(nodes, walls, by=c("x1", "x2", "y1", "y2"))
    nodes <- nodes %>%
      arrange(sorting)
  }
  
  id_aerenchyma <- unique(nodes$id_cell[nodes$type %in% c("aerenchyma", "inter_cellular_space")])
  id_aerenchyma <- id_aerenchyma-1
  
  print(Sys.time()-t_1)
  
  return(list(nodes = nodes,
              walls_nodes = walls,
              walls = wally,
              cells=all_cells,
              output = output,
              id_aerenchyma = id_aerenchyma))
  
}

pack_xylem <- function(params,nparenchyma, center){
  
  
  areas <- c(rnorm(round(nparenchyma),params$value[params$name == "stele" & params$type == "mean_size"],params$value[params$name == "stele" & params$type == "SD"]), rbeta(params$value[params$name == "xylem" & params$type == "n_files"],2,4)* params$value[params$name == "xylem" & params$type == "max_size"]) #rnorm (number, mean, SD) rbeta (number, shape1, shape2)* area max
  areas<-as.data.frame(areas)
  areas$type="NA"
  areas$type[1:round(nparenchyma)]="parenchyma"
  areas$type[round(nparenchyma+1):nrow(areas)]="xylem"
  
  #minxyl<-min(areas[round(nparenchyma+1):length(areas)])
  #maxpar<-max(areas[1:round(nparenchyma)])
  r_phlo<-params$value[params$name == "stele" & params$type == "cell_diameter"]* params$value[params$name == "stele" & params$type == "n_layers"]
  r_par<-r_phlo*(1-params$value[params$name == "phloem" & params$type == "proportion"])
  aire_par<-pi*(r_par^2)
  aire_phlo<-(pi*(r_phlo^2))-aire_par
  aire_tot<-pi*(r_phlo^2)
  fact<-aire_phlo/aire_tot
 
  if(nparenchyma>=1000){areas_phlo<-areas[1:round(nrow(areas)*fact),]
  areas_phlo$type="phloem" }
  if(nparenchyma<1000){
    if(fact>=0.4){areas_phlo<-areas[1:round(nrow(areas)*fact),]
    areas_phlo$type="phloem"}
    if(fact<0.4){areas_phlo<-areas[1:round(nrow(areas)*fact),]
    areas_phlo$type="phloem"
    areas_par<-areas[nrow(areas_phlo):round(nrow(areas)*0.4),]
    areas_par$type="parenchyma"
    areas_phlo<-rbind(areas_par, areas_phlo)}}
    
  
  if(nparenchyma>=1000){areas_par_layer<-areas[nrow(areas_phlo):round(nrow(areas_phlo)*1.7),]
  areas_par_layer$type="parenchyma"}
  if(nparenchyma<1000){areas_par_layer<-areas[nrow(areas_phlo):round(nrow(areas)*0.7),]
  areas_par_layer$type="parenchyma"}
  areas_xyl<-tail(areas,n=(nrow(areas)-(nrow(areas_phlo)+ nrow(areas_par_layer))))
  
  
  #areas_xyl<-areas[round(nrow(areas)*0.15)+1:nrow(areas),]
  areas_xyl <- sample_n(areas_xyl, nrow(areas_xyl))
  areas<-rbind(areas_xyl,areas_par_layer,areas_phlo)
  
  areas$areas<-sqrt(areas$areas * areas$areas)
  areas$id<-1:nrow(areas)
  # Generate the layout 
  packing <- circleProgressiveLayout(areas,sizecol ="areas")
  packing<-cbind(packing,areas)
  dat.gg <- circleLayoutVertices(packing, npoints=25, sizetype = "radius")
  dat.gg<-full_join(dat.gg,areas,by="id")
  
  ggplot(data = dat.gg) +
    geom_polygon(aes(x, y, group=id, fill = type), alpha = 0.7) + coord_fixed()
  
  packing <- packing%>%
    mutate(radius = radius/1000,
           x = x/1000 + center,
           mx = x,
           y= y /1000+ center,
           my = y,
           area = pi*radius^2,
           type = type,
           id_cell = id)
  
  xyl <- packing[packing$type == "xylem", ]
  xyl$id_group = 1:nrow(xyl)
  xyl$id_layer = 1.5
  circus <- seq(-0.95,0.95,0.95/4)
  cir <- data.frame(x_cir = rep(circus,2*nrow(xyl)))%>%
    mutate(y_cir = rep(c(sqrt(1-circus^2),-sqrt(1-circus^2)),nrow(xyl)),
           # mx = rep(cor$mx,2*length(circus)),
           # my = rep(cor$my,2*length(circus)),
           id_group = sort(rep(1:nrow(xyl),2*length(circus))))
  
  scaling = 0.95
  xyl_frontier <- merge(cir, xyl[,c("id_group", "radius", "mx", "my", "id_layer")], by = "id_group")%>%
    transmute(radius = radius,
              x = x_cir*radius*scaling+mx, # Check scaling in rondicortex (0.95)
              y = y_cir*radius*scaling+my,
              euc = sqrt((mx-center)^2+(my-center)^2),
              angle = ifelse(my-center > 0,acos((mx - center)/euc),
                             2*pi-acos((mx - center)/euc)) ,
              id_layer = id_layer,
              id_cell = 1,
              type = "xylem",
              order = params$value[params$name == "xylem" & params$type == "order"],
              id_group = id_group
    )%>%
    select(-euc)
  
  xyl_frontier%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group)))+
    coord_fixed()+
    guides(colour = F)
  
  packing%>%
    ggplot()+
    geom_point(aes(x, y , colour = type))+
    coord_fixed()
  head(packing)
  # filter without xyl and check column names
  packing_par <- subset(packing, packing$type =="parenchyma" |packing$type =="phloem"|packing$type =="test") %>%
    transmute(radius = radius,
              x=x,
              y=y,
              euc = sqrt((x-center)^2+(y-center)^2),
              angle = ifelse(y-center > 0,acos((x - center)/euc),
                             2*pi-acos((x - center)/euc)),
              id_layer = 1.5,
              id_cell = 1,
              type = type,
              order = 1,
              id_group = 0
    )%>%
    select(-euc)
  
  packing_par$id_cell<-1:length(packing_par$id_cell)
  xyl_frontier$id_cell<-(length(packing_par$id_cell)+1):(length(packing_par$id_cell)+length(xyl_frontier$id_cell))
  packing <- rbind (packing_par , xyl_frontier)
  return(packing)}


#' Read the parameters for GRANAR
#'
#' Read the parameters for GRANAR from an XML file
#' @param path The path to the XML files ith the parameters
#' @keywords granar
#' @export
#' @import xml2
#' @examples
#' read_param_xml()
#'

read_param_xml <- function(path = NULL){
  
  if( is.null(path) ){
    warning("No path specified")
  }
  
  input <- read_xml(path)
  params <- NULL
  
  # Quality checks. Check if all the needed tags are present in the XML file
  to_find <- c("planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex","secondarygrowth")
  for(tf in to_find){
    if (length(xml_find_all(input, paste0("//",tf))) == 0) warning(paste0("Could not find the '",tf,"' tag in the XML file"))
  }
  
  #Read the file and get the parameters in a table
  for( ch in xml_children(xml_find_all(input, "//*"))){
    att <- xml_attrs(ch)
    for(i in c(1:length(att))){
      params <- rbind(params, data.frame(
        name = xml_name(ch),
        type = names(att)[i],
        value = att[i]
      ))
    }
  }
  row.names(params) <- NULL
  params <- params %>% mutate(value = as.numeric(as.character(value)))
  
  return(params)
}


get_root_section <- function(path){
  
  if(is.null(path) ){warning("No path specified")}
  x <- read_xml(path)
  
  # GET THE CELLS DATA
  mydata <- xml_children(x)[2]
  temp <- xml_find_all(x, ".//cell")
  cells <- NULL
  for(i in c(1:length(temp))){
    wall <- xml_find_all(temp[i], ".//wall")
    n <- length(wall)
    cells <- rbind(cells, data.frame(id_cell = rep(as.numeric(xml_attr(temp[i], "id")), n),
                                     group = rep(as.numeric(xml_attr(temp[i], "group")), n),
                                     wall_id = as.numeric(xml_attr(wall, "id"))))
  }
  
  # GET THE WALL DATA
  mydata <- xml_children(x)[3]
  temp <- xml_find_all(mydata, ".//wall")
  walls <- NULL
  for(i in c(1:length(temp))){
    points <- xml_find_all(temp[i], ".//point")
    n <- length(points)
    walls <- rbind(walls, data.frame(id_cell = rep(as.numeric(xml_attr(temp[i], "id")), n),
                                     group = rep(as.numeric(xml_attr(temp[i], "group")), n),
                                     x = as.numeric(xml_attr(points, "x")),
                                     y = as.numeric(xml_attr(points, "y"))))
  }
  
  
  # GET THE GROUP INFORMATIONS
  mydata <- xml_find_all(x, ".//group")
  groups <- data.frame(id_cell = as.numeric(xml_attr(mydata, "id")),
                       type = xml_attr(mydata, "name")
  )
  
  #MERGE THE DATA
  rs <- merge(cells, walls, by.x = "wall_id", by.y = "id_cell")
  rs <- merge(rs, groups, by.x = "group.x", by.y = "id_cell")
  
  # REORDER THE POINTS WITHIN EACH CELL tO HAVE A NICE POLYGON
  rs1 <- NULL
  for(i in unique(rs$id_cell)){
    temp <- rs[rs$id_cell == i,]
    temp$atan <- atan2(temp$y - mean(temp$y), temp$x - mean(temp$x))
    temp <- temp[order(temp$atan),]
    temp <- temp[!duplicated(temp$atan),]
    rs1 <- rbind(rs1,temp)
  }
  root_section <- rs1[,c("id_cell", "x", "y", "type")]
  return(root_section)
}

aer_in_geom_xml <- function(sim, path = "C:/Users/schoppach/Documents/Tomate/MECHA/Projects/GRANAR/in/Tomato_Geometry.xml"){
  
  require(xml2)
  if (is.null(path)) {
    warning("No path specified")
  }
  if(is.null(sim$id_aerenchyma)){
    warning("No aerenchyma id specified")
  }else{
    id_aerenchyma <- sim$id_aerenchyma
  }
  xml <- read_xml(path)
  
  aer <- xml_children(xml_find_all(xml, "//aerenchyma_range"))
  
  # newbee <- 'aerenchyma id="0"'
  if(length(sim$id_aerenchyma) == 0){
    message("no aerenchyma")
  }else{
    new_siblings <- paste0('aerenchyma id="',id_aerenchyma,'"')
    xml_add_sibling(aer, new_siblings)
  }
  
  xml_remove(aer[1])
  
  path <- paste0(c(unlist(str_split(path, ".xml"))[1]),"_aer.xml")
  
  write_xml(xml, path)
  return(TRUE)
}

cell_layer <- function(params){
  layers <- params %>%
    filter(type %in% c("cell_diameter","n_layers","order")) %>%
    spread(type, value) %>%
    filter(!is.na(n_layers)) %>%
    arrange(order)
  stele_diameter <- params$value[params$name == "stele" & params$type == "layer_diameter"]
  
  # Create and "outside" layer to serve as boundary for the voronoi algorithm.
  layers <- rbind(layers, data.frame(name="outside",
                                     n_layers=2,
                                     cell_diameter=layers$cell_diameter[layers$name == "epidermis"]* 1,
                                     order = max(layers$order)+1))
  
  # Get the number of cell layers for the stele
  layers$n_layers[layers$name == "stele"] <- round((stele_diameter/2) / layers$cell_diameter[layers$name == "stele"]) #
  #layers$size[layers$name == "stele"] <- diam_stele
  
  
  # Get one row per actual cell layer
  all_layers <- NULL
  for(i in c(1:nrow(layers))){
    if(layers[i,"n_layers"] != 0){
      for(j in c(1:layers$n_layers[i])){
        all_layers <- rbind(all_layers, layers[i,])
      }
    }
  }
  
  all_layers <- layer_info(all_layers)
  
  return(list(all_layers = all_layers , layers = layers))
}

layer_info <- function(all_layers){
  
  all_layers$radius <- all_layers$cell_diameter / 2
  all_layers$perim <- all_layers$radius * 2 * pi
  all_layers$n_cell <- 1
  all_layers$angle_inc <- 0
  
  all_layers$radius[1] <- 0
  multi <- 1
  
  for(i in c(2:nrow(all_layers))){
    # Update radius
    all_layers$radius[i] <- all_layers$radius[i-1] + all_layers$cell_diameter[i-1] / 2 + all_layers$cell_diameter[i] / 2
    
    if(all_layers$name[i] == "pericyle" & multi == 1){
      all_layers$radius[i] <- stele_diameter/2 + all_layers$cell_diameter[i] / 2
      multi <- 2 # in case there is more than one pericycle layer
    }
    # Update perimeter
    all_layers$perim[i] <- all_layers$radius[i] * 2 * pi
    
    # Update number of cells in the layers
    all_layers$n_cell[i] <- round(all_layers$perim[i] / all_layers$cell_diameter[i])
    
    # Update the mean angle between cells
    all_layers$angle_inc[i] <- 2 * pi / all_layers$n_cell[i]
  }
  
  return(all_layers)
}

create_cells <- function(all_layers, random_fact){
  
  center <- max(all_layers$radius)
  all_cells <- NULL
  k <- 1
  for(i in c(1:nrow(all_layers))){
    radius <- all_layers$radius[i]
    if(all_layers$angle_inc[i] > 0){
      angles <- seq(from = 0, to = (2*pi), by = all_layers$angle_inc[i])[-1]
    }else{
      angles <- 0
    }
    k1 <- k+all_layers$n_cell[i]-1
    ks <- c(k:k1)
    k <- k1+1
    
    if(all_layers$name[i] == "outside"){
      x <- center + (radius * cos(angles))
      y <- center + (radius * sin(angles))
    }else if(all_layers$name[i] == "stele"){
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)
    }else if(substr(all_layers$name[i], 1,6) == "cortex"){
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact*3, random_fact*3)#* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact*3, random_fact*3)##* runif(all_layers$n_cell[i], 1-(random_fact*2), 1+(random_fact*2))
    }else{
      x <- center + (radius * cos(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-random_fact, 1+random_fact)
      y <- center + (radius * sin(angles)) + runif(all_layers$n_cell[i], -random_fact, random_fact)#* runif(all_layers$n_cell[i], 1-random_fact, 1+random_fact)
    }
    
    all_cells <- rbind(all_cells, data.frame(
      angle = angles,
      radius = radius,
      x = x,
      y = y,
      id_layer = i,
      id_cell = ks,
      type = all_layers$name[i],
      order = all_layers$order[i]
    )
    )
  }
  
  return(all_cells)
}

make_pith <- function(all_cells, params, center){
  if(params$value[params$name == "pith"][1] > 0){
    pith_size <- params$value[params$name == "pith" & params$type == "layer_diameter"]/2
    pcell <- params$value[params$name == "pith" & params$type == "cell_diameter"]
  }else{pith_size <- 0}
  
  if(pith_size > 0){
    
    
    xylem <- all_cells%>%
      filter(type == "xylem")
    
    xylem <- xylem%>%
      dplyr::group_by(id_group)%>%
      dplyr::mutate(mx = mean(x),
                    my = mean(y),
                    euc = sqrt((mx-center)^2+(my - center)^2))
    inner <- unique(xylem$id_group[xylem$euc < pith_size])
    all_cells <- all_cells[all_cells$type != "xylem" | all_cells$id_group %!in% inner,]
    
    
    n_pith_lay <- round(1+(pith_size-pcell/2)/pcell)
    pith_layer <- data.frame(name="stele",
                             n_layers=rep(n_pith_lay, n_pith_lay),
                             cell_diameter=pcell,
                             order = 0.5)
    
    pith_layer <- layer_info(pith_layer)
    new_cells <- create_cells(all_layers = pith_layer, random_fact = 0.001)
    new_cells%>%
      ggplot()+geom_point(aes(x,y))+
      coord_fixed()
    new_center <- mean(new_cells$x[new_cells$angle == 0], new_cells$y[new_cells$angle == 0])
    
    new_cells$x <- new_cells$x-new_center+center
    new_cells$y <- new_cells$y-new_center+center
    new_cells$id_group <- 0
    
    all_cells <- all_cells[sqrt((all_cells$x-center)^2+(all_cells$y-center)^2) > pith_size,]
    all_cells <- rbind(new_cells, all_cells)
    all_cells$id_cell <- 1:nrow(all_cells)
    
    all_cells%>%
      # filter(id_group %!in% inner)%>%
      ggplot()+
      geom_point(aes(x,y, colour = type))+
      # geom_point(aes(x,y), colour = "red", alpha = 0.2, data =     xylem%>%
      #              filter(id_group %in% inner))+
      # geom_point(aes(x,y), colour = "green", data = new_cells)+
      coord_fixed()
  }
  
  return(all_cells)
}

rondy_cortex <- function(params, all_cells, center){
  random_fact <- params$value[params$name == "randomness"] / 10
  cor_d <- params$value[params$name == "cortex" & params$type == "cell_diameter"]
  
  # calibration parameter
  # to_adjust1 <- params$value[params$name == "coefficient" & params$type == "icp_size"]
  # to_adjust2 <- params$value[params$name == "coefficient" & params$type == "icp_ratio"]
  
  all_cortex <- all_cells[all_cells$type %in% c("cortex","endodermis", "exodermis"),]
  
  icp_size <- params$value[params$name == "inter_cellular_space" & params$type == "size"]
  if(length(icp_size)> 0){
    scaling <- 1-(icp_size/cor_d)
    if(scaling >= 0.99){
      scaling = 0.99
    }
  }else{scaling <- 0.95}
  
  
  if(length(all_cells$id_group[all_cells$type == "xylem"]) > 0){
    k_max_xylem <- max(all_cells$id_group[all_cells$type == "xylem"])
  }else{k_max_xylem <- 0}
  
  ctess <- deldir(all_cortex$x, all_cortex$y, digits = 8)
  idc <- unique(all_cortex$id_cell)
  idc <- 1:length(idc)
  rc <- ctess$dirsgs[ctess$dirsgs$ind1 %in% idc |
                       ctess$dirsgs$ind2 %in% idc,]
  rc <- rc%>% arrange(ind1)
  rc2 <- data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind1)
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind1))
  rc2 <- rbind(rc2, data.frame(x = rc$x2, y=rc$y2, id_cell = rc$ind2))
  rc2 <- rbind(rc2, data.frame(x = rc$x1, y=rc$y1, id_cell = rc$ind2))
  
  inner <- min(all_cortex$radius[all_cortex$type %in% c("cortex")]) # [all_cortex$type %in% c("endodermis", "cortex")]
  outer <- max(all_cortex$radius[all_cortex$type %in% c("cortex")]) # no intercellular space between exo cortex and cortex
  
  rc2 <- rc2%>%mutate(euc = sqrt((x-center)^2+(y-center)^2))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(mx = mean(x),
                  my = mean(y),
                  atan = atan2(y-my, x - mx))%>%
    dplyr::arrange(id_cell, atan)
  
  all_cortex$id_cell <- 1:nrow(all_cortex)
  rc1 <- merge(rc2, all_cortex[,c("id_cell", "type", "radius", "id_layer")], by="id_cell")
  
  rcin <- rc2%>%
    filter(euc > inner,
           euc < outer)%>%
    mutate(ID = paste0(x,y))%>%
    filter(!duplicated(ID))
  
  all_inter <- data.frame(angle = ifelse(rcin$y-center >= 0, acos((rcin$x - center)/rcin$euc),
                                         2*pi-acos((rcin$x - center)/rcin$euc)),
                          radius = rcin$euc,
                          x = rcin$x, y = rcin$y,
                          id_layer = all_cortex$id_layer[1]+0.5,
                          id_cell = 1:nrow(rcin),
                          type = "inter_cellular_space",
                          order = params$value[params$name == "cortex" & params$type == "order"]+0.5,
                          id_group = 0 # if too close, they should be merge but not now
  )
  if(length(params$value[params$name =="inter_cellular_space"]) > 0){
    coef_icp <- (10 * params$value[params$name =="inter_cellular_space" & params$type == "ratio"]) # here to modulate icp proportion coeficient
    if(coef_icp > 1){coef_icp = 1}
    inter_cellular_proportion <- coef_icp*nrow(all_inter)
  }else{
    inter_cellular_proportion <- 0.5*nrow(all_inter)
  }
  if(inter_cellular_proportion == 0){
    all_inter <- NULL
  }else{
    to_keep <- sample(1:nrow(all_inter), round(inter_cellular_proportion), replace=F)
    all_inter <- all_inter[all_inter$id_cell %in% to_keep,]
    all_inter$id_point <- paste0(all_inter$x,";",all_inter$y)
    all_inter <- all_inter%>%
      filter(!duplicated(id_point))%>%
      select(-id_point)
    all_inter%>%
      ggplot()+
      geom_point(aes(x,y))+
      coord_fixed()
  }
  
  nodes <- vertex(rc1%>%
                    filter(type == "cortex"))
  nodes <- nodes %>%
    filter(wall_length > 0)%>%
    mutate(m = (y2-y1)/(x2-x1),
           k = y1-m*x1,
           r_dist = abs(k+m*mx-my)/sqrt(1+m^2)) # distance between a point (mx,my) to a segment defined by to point (x1,y1; x2,y2)
  
  nodes%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2))+
    coord_fixed()
  
  cor <- nodes%>%
    filter(wall_length > 0)%>%
    dplyr::group_by(id_cell)%>%
    dplyr::mutate(radius = min(r_dist))%>%
    filter(!duplicated(id_cell))
  
  cor$radius[cor$id_layer %in% c(inner, outer)] <- cor$radius[cor$id_layer %in% c(inner, outer)]*0.45
  
  cor$id_group = 1:nrow(cor)
  
  circus <- seq(-0.95,0.95,0.95/4)
  cir <- data.frame(x_cir = rep(circus,2*nrow(cor)))%>%
    mutate(y_cir = rep(c(sqrt(1-circus^2),-sqrt(1-circus^2)),nrow(cor)),
           # mx = rep(cor$mx,2*length(circus)),
           # my = rep(cor$my,2*length(circus)),
           id_group = sort(rep(1:nrow(cor),2*length(circus))))
  
  cor_frontier <- merge(cir, cor[,c("id_group", "radius", "mx", "my", "id_layer")], by = "id_group")%>%
    transmute(radius = radius,
              x = x_cir*radius*scaling+mx,
              y = y_cir*radius*scaling+my,
              euc = sqrt((mx-center)^2+(my-center)^2),
              angle = ifelse(my-center > 0,acos((mx - center)/euc),
                             2*pi-acos((mx - center)/euc)) ,
              id_layer = id_layer,
              id_cell = 1,
              type = "cortex",
              order = params$value[params$name == "cortex" & params$type == "order"],
              id_group = id_group
    )%>%
    select(-euc)
  
  cor_frontier%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group)))+
    coord_fixed()+
    guides(colour = F)
  
  all_cells <- rbind(all_cells[all_cells$type != "cortex",], cor_frontier)
  all_cells$id_group[all_cells$type == "cortex"] <- all_cells$id_group[all_cells$type == "cortex" & all_cells$id_group != 0] + k_max_xylem
  all_cells <- rbind(all_cells, all_inter)
  
  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))
  return(all_cells)
}

vascular <- function(all_cells, params, layers, center){
  
  
  n_xylem_files <- params$value[params$name == "xylem" & params$type == "n_files"]
  proto_meta_ratio <- params$value[params$name == "xylem" & params$type == "ratio"]
  n_proto_xylem <- round(n_xylem_files*proto_meta_ratio)
  plant_type <- params$value[params$name == "planttype"]
  if(length(all_cells$id_group[all_cells$type == "cortex"])> 0){
    k_max_cortex <- max(all_cells$id_group[all_cells$type == "cortex"])
  }else{k_max_cortex  <- 0}
  
  # length(all_cells$id_cell[all_cells$type == "inter_cellular_space"])
  
  all_cells%>%
    filter(id_group != 0)%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group)))+coord_fixed()+guides(colour = F)
  
  if(plant_type == 2){ # DICOT
    xyl <- data.frame(r=numeric(2), d=numeric(2))
    xyl$r <- c(0, max(all_cells$radius[all_cells$type == "stele"]))
    xyl$d <- c(params$value[params$type == "max_size" & params$name == "xylem"], layers$cell_diameter[layers$name == "stele"])
    
    # Get the cells in between
    fit <- lm(d ~ r, data=xyl)$coefficients
    rnew <- xyl$r[1]
    i <- 1
    rmax <- xyl$r[2]
    dmin <- xyl$d[2]
    keep_going <- T
    while(keep_going){
      xyl <- xyl %>% arrange(r)
      
      rnew <- xyl$r[i] + xyl$d[i] #+ (xyl$r[2]/10)
      dnew <- fit[1] + rnew*fit[2]
      while(rnew+(dnew/2) > rmax-(dmin/2)){
        rnew <- rnew - 0.05
        dnew <- dnew - 0.05
        keep_going = F
      }
      xyl <- rbind(xyl, data.frame(r = rnew,d = dnew))
      i <- i+1
    }
    xyl <- xyl %>% arrange(r)%>%
      filter(d > 0)
    while(xyl$d[nrow(xyl)] >= xyl$d[nrow(xyl)-1]){
      xyl$d[nrow(xyl)] <- xyl$d[nrow(xyl)] - 0.04
      xyl$r[nrow(xyl)] <- xyl$r[nrow(xyl)] - 0.02
    }
    all_xylem <- NULL
    i <- 1
    angle_seq <- seq(from = 0, to = (2*pi), by = (2 * pi) / n_xylem_files)
    x <- center + (xyl$r[1] * cos(angle_seq[1]))
    y <- center + (xyl$r[1] * sin(angle_seq[1]))
    all_xylem <- rbind(all_xylem, data.frame(x = x,
                                             y = y,
                                             d = xyl$d[1],
                                             angle = angle_seq[1],
                                             id_group = i))
    all_cells <- rbind(all_cells, data.frame(
      angle = angle_seq[1],
      radius = xyl$r[1],
      x = x,
      y = y,
      id_layer = 20,
      id_cell = 1,
      type = "xylem",
      order = 1.5,
      id_group = i
    ))
    
    i <- i+1
    for(angle in angle_seq){
      x <- center + (xyl$r[-1] * cos(angle))
      y <- center + (xyl$r[-1] * sin(angle))
      all_xylem <- rbind(all_xylem, data.frame(x = x,
                                               y = y,
                                               d = xyl$d[-1],
                                               angle = angle,
                                               id_group = i))
      all_cells <- rbind(all_cells, data.frame(
        angle = angle,
        radius = xyl$r[-1],
        x = x,
        y = y,
        id_layer = 20,
        id_cell = 1,
        type = "xylem",
        order = 1.5,
        id_group = i
      )
      )
      i <- i+1
    }
    # Phloem vessels are built between xylem ones
    phl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                      d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_ph <- seq(from = ((2 * pi) / n_xylem_files ) /2, to = (2*pi), by = (2 * pi) / n_xylem_files)
    for(angle in angle_seq_ph){
      x1 <- center + (phl$r[1] * cos(angle))
      y1 <- center + (phl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a phloem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "phloem", type))
      
      # Get the compagnion cells
      for (compa in 1:2) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
      }
      
      for (compa in 1:6) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "cambium", type))
      }
      all_cells <- all_cells%>% select(-dist_phl)
      
    }
    
  }else if(plant_type == 1){ # MONOCOT
    if(n_xylem_files == 1){ # One metaxylem in the center of the stele
      r <- 0
      xyl <- data.frame(r = r,
                        d = params$value[params$type == "max_size" & params$name == "xylem"])
    }else{
      
      #modification 05/01
      r= max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])*1.5 -
        (params$value[params$type == "max_size" & params$name == "xylem"])/2
      xyl <- data.frame(r = r,
                        d = params$value[params$type == "max_size" & params$name == "xylem"])
    }
    #  xyl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "max_size" & params$name == "xylem"])/2,
    #                  d = params$value[params$type == "max_size" & params$name == "xylem"])
    all_xylem <- NULL
    angle_seq <- seq(from = 0, to = (2*pi)-(2 * pi) / n_xylem_files, by = (2 * pi) / n_xylem_files)
    i <- 1
    for(angle in angle_seq){
      x <- center + (xyl$r[1] * cos(angle))
      y <- center + (xyl$r[1] * sin(angle))
      all_xylem <- rbind(all_xylem, data.frame(x = x,
                                               y = y,
                                               d = xyl$d[1],
                                               angle = angle,
                                               id_group = i))
      all_cells <- rbind(all_cells, data.frame(
        angle = angle,
        radius = xyl$r[1],
        x = x,
        y = y,
        id_layer = 20,
        id_cell = 1,
        order = 1.5,
        type = "xylem",
        id_group = i))
      i <-i+1
    }
    
    # remove stele cell inside metaxylem vessels
    for(i in c(1:nrow(all_xylem))){
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele"))
    }
    
    # Make circular frontier for xylem
    x_cir <- seq(-0.95,0.95,0.95/4)
    y_p <- sqrt(1-x_cir^2)
    y_m <- -sqrt(1-x_cir^2)
    xyl_frontier <- NULL
    k <- 1
    for (i_xyl in 1:nrow(all_xylem)) {
      tmp <- all_xylem[i_xyl,]
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs(tmp$d*0.8)/2
      cir$x <- cir$x + tmp$x
      cir$y <- cir$y + tmp$y
      xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                     radius = xyl$r[1],
                                                     x = cir$x,
                                                     y = cir$y,
                                                     id_layer = 20,
                                                     id_cell = 1,
                                                     type = "xylem",
                                                     order = 1.5,
                                                     id_group = k))
      
      k <- k + 1
    }
    xyl_frontier%>%
      ggplot()+
      geom_point(aes(x,y, colour = factor(id_group)))+
      geom_point(aes(x,y, colour = factor(id_group)), data = all_xylem)+
      coord_fixed()
    # add xyl frontier
    all_cells <- all_cells[all_cells$type != "xylem",]
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex
    
    # protoxylem vessels are built on the outer stele rim
    protoxyl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                           d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_proto <- seq(from = 0, to = (2*pi)-(2 * pi) / n_proto_xylem, by = (2 * pi) / n_proto_xylem)
    for(angle in angle_seq_proto){
      x1 <- center + (protoxyl$r[1] * cos(angle))
      y1 <- center + (protoxyl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a protoxylem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_protoxyl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_protoxyl = ifelse(type == "stele", dist_protoxyl, 100)) %>%
        mutate(type = ifelse(dist_protoxyl == min(dist_protoxyl), "xylem", type))
    }
    
    # Phloem vessels are built between xylem ones
    phl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                      d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_ph <- seq(from = ((2 * pi) / n_proto_xylem ) /2, to = (2*pi), by = (2 * pi) / n_proto_xylem)
    for(angle in angle_seq_ph){
      x1 <- center + (phl$r[1] * cos(angle))
      y1 <- center + (phl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a phloem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "phloem", type))
      
      # Get the compC"nion cells
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
      
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
      
      
    }
    all_cells <- all_cells%>%select(-dist_phl, -dist_protoxyl)
    
  }
  
  # Change the identity of stele cells to be replaced by xylem cells
  if(plant_type == 2){
    # all_xylem <- all_xylem%>%
    #   filter(d < mean(all_xylem$d))
    # all_cells <- all_cells%>%
    #   filter(type != "xylem" | radius > mean(all_cells$radius[all_cells$type == "xylem"]) )
    #
    # all_xylem%>%
    #   ggplot()+
    #   geom_point(aes(x,y, colour = d))+
    #   coord_fixed()
    
    for(i in c(1:nrow(all_xylem))){
      # print(i)
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
    }
  }
  
  if(plant_type == 2){
    #
    # xyl <- xyl[xyl$r < mean(xyl$r),]
    
    
    all_xylem <- all_cells[all_cells$type == "xylem",]%>%
      mutate(radius = round(radius,2),
             d = c(xyl$d[1], rep(xyl$d[2:length(xyl$d)],
                                 max(unique(all_cells$id_group[all_cells$type == "xylem"]))-1)))# %>%
    # Remove xylem cell that are to close to each other
    # filter(!duplicated(angle), !duplicated(radius))
    
    all_cells <- all_cells[all_cells$type != "xylem",]
    
    # Make circular frontier for xylem
    x_cir <- seq(-0.95,0.95,0.95/4)
    y_p <- sqrt(1-x_cir^2)
    y_m <- -sqrt(1-x_cir^2)
    xyl_frontier <- NULL
    k <- 1
    for (i_xyl in 1:nrow(all_xylem)) {
      tmp <- all_xylem[i_xyl,]
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs(tmp$d*0.8)/2
      cir$x <- cir$x + tmp$x
      cir$y <- cir$y + tmp$y
      xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                     radius = tmp$radius,
                                                     x = cir$x,
                                                     y = cir$y,
                                                     id_layer = 20,
                                                     id_cell = 1,
                                                     type = "xylem",
                                                     order = 1.5,
                                                     id_group = k))
      xyl_frontier%>%
        ggplot()+
        geom_point(aes(x,y, colour = factor(id_group)))+
        coord_fixed()
      k <- k + 1
    }
    all_cells <- rbind(all_cells, xyl_frontier)
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex
    
    
    if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
      all_cells <- make_pith(all_cells, params, center)
    }
    
  }
  all_cells %>%
    filter(id_group != 0)%>%
    ggplot()+
    geom_point(aes(x,y, colour = id_group))+
    coord_fixed()
  
  # reset the cell ids
  all_cells$id_cell <- c(1:nrow(all_cells))
  
  return(all_cells)
}

cell_voro <- function(all_cells, vtess, center){
  
  # Get the size of the cells
  #cell_size <- vtess$summary
  ids <- all_cells$id_cell
  #cell_size$id_cell <- c(1:nrow(cell_size))
  #all_cells <- merge(all_cells, cell_size[,c("id_cell", "dir.area")], by="id_cell")
  all_cells$area <- NA #all_cells$dir.area
  all_cells$dist <- sqrt((all_cells$x - center)^2 + (all_cells$y - center)^2 )
  
  rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                       vtess$dirsgs$ind2 %in% ids,]
  
  # Get the cooridnates for every nodes in the voronoi
  rs <- rs %>% arrange(ind1)
  rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
  rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))
  
  rs2 <- merge(rs2, all_cells[,c("id_cell", "type", "area", "dist", "angle", "radius", "id_layer", "id_group")], by="id_cell")
  rs2 %>%
    filter(type %in% c("cortex", "xylem", "inter_cellular_space"))%>%
    ggplot()+
    geom_point(aes(x,y, colour = factor(id_group), shape = type))+
    coord_fixed()+guides(colour = F)
  rs2 <- rs2%>%filter(type != "outside")
  
  return(list(all_cells = all_cells, rs2 = rs2))
  
}

vertex <- function(rs1){
  nodes <- rs1 %>%
    mutate(id_point= paste0(rs1$x,";",rs1$y))%>%
    group_by(id_cell) %>%
    filter(!duplicated(id_point))%>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))%>%
    select(-id_point)
  
  nodes <- nodes %>%
    ungroup() %>%
    mutate(vertical = ifelse(x == xx, "true", "false")) %>%
    mutate(x1 = ifelse(x > xx, x, xx)) %>%
    mutate(x2 = ifelse(x > xx, xx, x)) %>%
    mutate(y1 = ifelse(x > xx, y, yy)) %>%
    mutate(y2 = ifelse(x > xx, yy, y)) %>%
    # Bug fix when wall is perfectly vertical
    mutate(y1 = ifelse(x == xx,
                       ifelse(y > yy, yy, y), y1)) %>%
    mutate(y2 = ifelse(x == xx,
                       ifelse(y > yy, y, yy), y2)) %>%
    mutate(wall_length = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
    mutate(wall_length2 = sqrt((xx-x)^2 + (yy-y)^2),
           slope = (y2-y1)/(x2-x1),
           intercept = y1 - slope*x1)
  return(nodes)
}

smoothy_cells <- function(rs1){
  # saved_rs <- rs1
  
  # rs1 <- saved_rs
  for(i in c(1:max(rs1$id_group))){
    temp <- rs1[rs1$id_group == i,]
    rs1 <- rs1[rs1$id_group != i,]
    if(nrow(temp)> 0){
      temp <- concavety(temp)
      rs1 <- rbind(rs1,temp)%>%as.tibble()
    }
    
  }
  
  rs1 <- rs1%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    ungroup()
  voiz <- rs1%>%
    dplyr::group_by(id_point)%>%
    dplyr::summarise(n = n())%>%
    ungroup()
  
  rs1 <- rs1[rs1$id_point %!in% voiz$id_point[voiz$n < 3] | rs1$type == "epidermis", ]
  
  return(rs1)
}

concavety <- function(data){
  
  # data %>%
  #   ggplot()+
  #   geom_polygon(aes(x,y, group = id_cell, fill = id_cell), alpha = 0.2,colour = "white")+
  #   geom_point(aes(x,y), data = rock_roll)+
  #   #geom_polygon(aes(x,y),fill = "blue",alpha = 0.5, colour = "white", data = rock_roll)+
  #   coord_fixed()
  #  data <- saved_cells
  
  nodes <- data%>%
    mutate(id_point = paste0(x,";", y))
  id_list <- NULL
  rock_roll <- NULL
  while(length(id_list) != length(unique(data$id_cell))){
    nodes <- nodes%>%
      filter(id_cell %!in% id_list)%>%
      dplyr::group_by(id_cell)%>%
      filter(!duplicated(id_point))
    
    double <- nodes%>%
      dplyr::group_by(id_point)%>%
      dplyr::summarise(n = n())%>%
      filter(n > 1)
    
    
    first <- nodes[nodes$y == max(nodes$y),][1,]# select a point in the outer border
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$x == max(nodes$x),][1,]
    }
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$x == min(nodes$x),][1,]
    }
    if(first$id_point %in% double$id_point){
      first <- nodes[nodes$y == min(nodes$y),][1,]
    }
    if(first$id_point %in% double$id_point){
      stop("huston, we have a problem")
    }
    
    first_point <- first$id_point[1] # full incremental should start somewhere
    first_cell <- first$id_cell[1]
    
    swing <- first_cell
    it <- T
    rock <- first[1,]
    cinp <- nodes[nodes$id_cell == swing,]
    last_pos <- which(cinp$id_point == first_point)-1
    if(length(last_pos) == 1){
      END_point <- ifelse(cinp$id_point[1] == first_point, last(cinp$id_point),
                          cinp$id_point[last_pos])
    }else{message = "multiple end point in the concave shape object"}
    while(it){
      if(is.na(last(rock$id_point))){print("breaks due to NA values in polygons coords")
        stop()}
      # junction point
      if(last(rock$id_point) %in% double$id_point){ # is it a juction point
        # unless they have only one point as a connection between the cells
        battle <- unique(nodes$id_cell[nodes$id_point == last(rock$id_point)])
        battle <- nodes[nodes$id_cell %in% battle,]
        fr <- battle%>%
          dplyr::group_by(id_point)%>%
          dplyr::summarise(n = n())%>%
          filter(n > 1)
        if (nrow(fr) == 1){ # only one point for contact between cells
          # do nothing
        }else{
          # print(swing)
          # switch from the previous id_cell to the next
          swing <- unique(nodes$id_cell[nodes$id_point == last(rock$id_point) &
                                          nodes$id_cell != rock$id_cell[rock$id_point == last(rock$id_point)]])
          
          
          # print(swing)
          if(length(swing) > 1){ # when triple point or quadri point
            # next cell should be the next in a anti-clockwise order
            depa <- nodes[nodes$id_cell %in% swing,]
            last_y <- nodes$y[nodes$id_point == last(rock$id_point)][1]
            last_x <- nodes$x[nodes$id_point == last(rock$id_point)][1]
            
            prev_rock <- rock[rock$id_point != last(rock$id_point),]
            av_la <- prev_rock[prev_rock$id_point == last(prev_rock$id_point),][1,]
            av_la$dist_point <- sqrt((av_la$x-last_x)^2+(av_la$y-last_y)^2)
            crit_angle <- ifelse(av_la$y-last_y >= 0, acos((av_la$x - last_x)/av_la$dist_point),
                                 2*pi-acos((av_la$x - last_x)/av_la$dist_point))
            
            depa <- depa%>%
              mutate(dist_point = sqrt((mx-last_x)^2+(my-last_y)^2),
                     crit = ifelse(my-last_y >= 0, acos((mx - last_x)/dist_point),
                                   2*pi-acos((mx - last_x)/dist_point)))
            swing <- ifelse(crit_angle < min(depa$crit), depa$id_cell[depa$crit == min(depa$crit)], #
                            ifelse(crit_angle < max(depa$crit), depa$id_cell[depa$crit == max(depa$crit)],
                                   depa$id_cell[depa$crit == min(depa$crit)] )
            )
          }
          # print(swing) # only one value
          yep <- nodes[nodes$id_cell == swing & nodes$id_point == last(rock$id_point),][1,] # take the next point
          rock <- rbind(rock, yep)
        }
      }
      
      cinp <- nodes[nodes$id_cell == swing,] # in this cell selection
      roll <- rbind(cinp[cinp$atan >= last(rock$atan),],cinp[cinp$atan < last(rock$atan),]) # sort points
      # print(roll)
      if(roll$id_point[1] == END_point){
        # print("will finish convave hull shortly")
        it <- F
        break()
      }
      pn <- roll[roll$id_point %!in% rock$id_point,][1,] # take the next point
      #print(pn)
      if(is.na(pn$x)){
        print(roll)
        pl <- nodes%>%
          ggplot()+
          geom_polygon(aes(x,y, group = id_cell, fill = id_cell),colour = "white", alpha = 0.5)+
          geom_point(aes(x,y), size = 2, alpha = 0.5, data = nodes[nodes$id_point %in% rock$id_point,])+
          
          coord_fixed()
        
        print(pl)
      }
      rock <- rbind(rock, pn)
      tata <- pn$atan
      
      
    }
    
    rock$id_cell <- first_cell
    rock$atan <- seq(-pi,pi, 2*pi/nrow(rock))[1:nrow(rock)]
    tmp_id_list <- unique(nodes$id_cell[which(point.in.polygon(nodes$x, nodes$y, rock$x,rock$y) > 0)])
    id_list <- c(id_list, tmp_id_list)
    rock_roll <- rbind(rock_roll, rock)
    
    for (i in unique(rock_roll$id_cell)) {
      tmp_check <- rock_roll[rock_roll$id_cell == i,]
      tmp_check <- tmp_check[!is.na(tmp_check$x), ]
      pol <- Polygon(tmp_check[, c("x","y")])
      rock_roll$area[rock_roll$id_cell == i] <-  pol@area
    }
    if(0 %in% rock_roll$area){
      message(paste0("id ",unique(data$id_cell), "have null area, try to fix bug"))
    }
    #
    # pl <- rock_roll%>%
    #   ggplot(aes(x,y))+
    #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)), fill = "grey",alpha = 0.3,colour = "white", data = data)+
    #   geom_polygon(aes(x,y, group = id_cell, fill = factor(id_cell)),alpha = 0.6, colour = "white")+
    #   coord_fixed()+
    #   guides(fill = F)
    # print(pl)
    
  }
  
  
  return(rock_roll%>%select(colnames(data))%>%ungroup())
  
  
}

fuzze_inter <- function(rs1){
  # saved_rs <- rs1
  space <- rs1[rs1$type == "inter_cellular_space",]
  if(nrow(space)> 0){
    
    double <- space%>%
      dplyr::group_by(id_point)%>%
      dplyr::summarise(n = n())%>%
      filter(n > 1)
    to_correct <- unique(space$id_cell[space$id_point %in% double$id_point])
    
    if(nrow(double) > 0){
      done <- NULL
      itm <- 0
      for (btw in to_correct) {
        comu1 <- space$id_point[space$id_cell == btw & space$id_point %in% double$id_point]
        nei <- unique(space$id_cell[space$id_cell != btw & space$id_point %in% comu1])
        bou <- c(btw, nei)
        ke <- length(bou)
        te <- 0
        while(ke > te){
          te <- length(bou)
          comu1 <- space$id_point[space$id_cell %in% bou & space$id_point %in% double$id_point]
          nei <- unique(space$id_cell[space$id_cell %!in% bou & space$id_point %in% comu1])
          bou <- unique(c(bou, nei))
          ke <- length(bou)
        }
        itm <- bou
        
        if(itm[1] %in% done){next()}
        # print(bou)
        tmp_cell <- space[space$id_cell %in% itm, ]
        for (i in unique(tmp_cell$id_cell)) {
          tmp <- tmp_cell[tmp_cell$id_cell == i,]
          tmp <- tmp[!is.na(tmp$x), ]
          pol <- Polygon(tmp[, c("x","y")])
          tmp_cell$area[tmp_cell$id_cell == i] <-  pol@area
        }
        
        tmp_cell <- tmp_cell[tmp_cell$area > 0, ]
        if(nrow(tmp_cell) > 0){
          tmp_cell <- tmp_cell%>%
            dplyr::group_by(id_cell)%>%
            dplyr::filter(!duplicated(id_point))
          
          rs1 <- rs1[rs1$id_cell %!in% itm,]
          tmp_cell <- concavety(tmp_cell)
          rs1 <- rbind(rs1,tmp_cell)
        }else{
          rs1 <- rs1[rs1$id_cell %!in% itm,]
        }
        done <- c(done, itm)
      }
      
    }
  }
  
  return(rs1)
}

aerenchyma <- function(params, rs1){
  
  # saved_rs <- rs1
  # rs1 <- saved_rs
  
  aer_type <- params$value[params$name == "aerenchyma" & params$type == "type"]
  if(length(aer_type)== 0){
    aer_type <- params$value[params$name == "planttype" & params$type == "param"]
  }
  n_aerenchyma_files <- params$value[params$name == "aerenchyma" & params$type == "n_files"]
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]
  
  angle_inc <- (2 * pi) / n_aerenchyma_files
  safe_cortex_layer <- c(min(rs1$id_layer[rs1$type == "cortex"]))
  last_cortex_layer <- c(max(rs1$id_layer[rs1$type == "cortex"]))
  area_all <- rs1%>%
    filter(type %in% c("cortex", "inter_cellular_space", "exodermis", "epidermis"))%>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))
  ini_cortex_area <- sum(area_all$area)
  surface_to_kill <- ini_cortex_area*proportion_aerenchyma
  stk_zone = surface_to_kill/n_aerenchyma_files
  # (2 * pi * proportion_aerenchyma / 30) / n_aerenchyma_files
  if (aer_type == 1){
    small_r <- mean(rs1$dist[rs1$id_layer == safe_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == safe_cortex_layer])
    big_R <- mean(rs1$dist[rs1$id_layer == last_cortex_layer])+0.5*mean(rs1$radius[rs1$id_layer == last_cortex_layer])
    angle_range_inc <- stk_zone/(big_R^2-small_r^2)
    
  }
  if(aer_type == 2){
    angle_range_inc <- (2 * pi * proportion_aerenchyma / 100) / n_aerenchyma_files
  }
  
  angle_range_inc_ini <- angle_range_inc
  cortex_area <- ini_cortex_area
  angle <- runif(1, 0.6, 1) * pi/n_aerenchyma_files
  angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
  id_group_max <- max(rs1$id_group)
  
  
  for(j in c(1:n_aerenchyma_files)){
    # saved_rs1 <- rs1
    possi <- rs1[rs1$id_layer %!in% safe_cortex_layer & rs1$type %in% c("cortex", "inter_cellular_space"),]
    gotogo <- T
    n_try <- 0
    while(gotogo){
      angle_range <- c(angle - angle_range_inc, angle + angle_range_inc)
      ma_cell <- rs1 %>%
        filter(id_layer %!in% safe_cortex_layer ,
               type == "cortex" ,
               angle > angle_range[1] ,
               angle < angle_range[2])
      are <- sum_area(ma_cell)
      if(are > stk_zone){
        angle_range_inc <- angle_range_inc*0.9
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
        next()
      }
      if(nrow(ma_cell) > 0){
        while(are <= stk_zone ){
          are <- sum_area(ma_cell)
          pointy <- unique(ma_cell$id_point)
          done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
          nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
          icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"]) # neigh inter_cell
          pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)]) # merge potential and inter_cell
          
          done <- c(done, nei) # full potential
          ma_cell <- possi[possi$id_cell %in% done,]
          are <- sum_area(ma_cell)
          
          # pl <- ggplot()+
          #   geom_polygon(aes(x,y,group = id_cell, fill = id_cell), colour = "white", alpha = 0.8,data = ma_cell)+
          #   geom_polygon(aes(x,y,group = id_cell, fill = id_cell), colour = "white", alpha = 0.1,data = rs1)+
          #   geom_polygon(aes(x,y,group = id_cell), fill = "red",colour = "white", alpha = 0.1,data = rs1%>%
          #                  filter(id_layer == safe_cortex_layer))+
          #   coord_fixed()
          # print(pl)
          
          n_try <- n_try + 1
          if(n_try > 100){
            are = stk_zone
          }
        }
        # once it overflow the threshold limit
        ma_cell <- possi[possi$id_cell %in% pre,] # this default potential cell and inter_cell are selected
        are_b <- sum_area(ma_cell)
        # adjust to include the right amount of cells as aerenchyma
        if(are_b > stk_zone){to_be_added <- NULL}else{
          
          next_area <- possi[possi$id_cell %in% done & possi$id_cell %!in% pre,]
          alm <- next_area%>%
            dplyr::group_by(id_cell)%>%
            dplyr::summarise(area = mean(area))%>%
            mutate(sarea = cumsum(area),
                   near = abs(sarea+are_b-stk_zone))
          jk <- which(alm$near == min(alm$near))
          
          to_be_added <- alm$id_cell[1:jk]
          ma_cell <- possi[possi$id_cell %in% c(pre,to_be_added),]
          pointy <- unique(ma_cell$id_point)
          done <- unique(ma_cell$id_cell) # id from potential cell to be named aer
          nei <- unique(possi$id_cell[possi$id_cell %!in% done & possi$id_point %in% pointy]) # neigh cells and inter_cell
          icp <- unique(possi$id_cell[possi$id_cell %in% nei & possi$type == "inter_cellular_space"])
          pre <- unique(possi$id_cell[possi$id_cell %in% c(done,icp)])
          ma_cell <- possi[possi$id_cell %in% pre,]
        }
        
        rs1 <- rs1[rs1$id_cell %!in% pre,]
        gotogo <- F # we are done for this part
        #ma_cell = saved_ma_cell <- ma_cell
        aer <- concavety(ma_cell)
        aer$type <- as.character("aerenchyma")
        aer$id_group <- j+id_group_max
        
        
        # print(unique(rs1$id_cell[rs1$id_cell %in% pre]))
        # print(str(aer))
        # print(str(rs1))
        rs1 <- rbind(rs1, aer)
        
        # print(aer%>%
        #   ggplot(aes(x,y))+
        #   geom_polygon(aes(x,y,group = id_cell, fill = type), colour = "white", alpha = 0.8,data = rs1)+
        #   geom_polygon(colour = "white", fill = "red")+
        #   coord_fixed())
        #
      }else{
        angle_range_inc <- angle_range_inc*1.1
        n_try <- n_try + 1
        if(n_try > 100){
          gotogo <- F
        }
      }
      
    }
    angle <- angle + angle_inc
    message(paste0(j,"/",n_aerenchyma_files))
    # j <- j +1
  }
  rs <- rs1
  
  return(rs)
  
}

sum_area <- function(cells){
  area_tmp <- cells %>%
    dplyr::group_by(id_cell)%>%
    dplyr::summarise(area = mean(area))%>%
    arrange(area, decreasing = T)
  are <- sum(area_tmp$area)
  return(are)
  
}

septa <- function(rs1){
  data <- rs1
  data$id_point <- paste0(data$x,";",data$y)
  data$neib1 = data$neib2 <- NA
  data$neib1_type = data$neib2_type <- NA
  
  for(i in which(data$type == "aerenchyma")){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- data$type[data$id_cell == ne[1]][1]
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- data$type[data$id_cell == ne[2]][1]
    if(!is.na(ne[3])){message("quadri point")
      print(ne)}
  }
  # all triple points and all points next to cortex cells
  must <- data[!is.na(data$neib2),]
  must <- rbind(must, data[data$neib1_type != "aerenchyma",])
  
  # add some noise
  noise <- data[!(data$neib1_type != "aerenchyma" | !is.na(data$neib2)),]
  noise_point <- unique(noise$id_point)
  n_noise <- round(length(noise_point)*0.05)
  noise_keep <- sample(noise_point, n_noise, replace=F)
  noise <- noise[noise$id_point %in% noise_keep,]
  must <- rbind(must, noise)
  
  data <- data[data$type != "aerenchyma",]
  data <- rbind(data,must)%>%
    arrange(id_cell,atan)
  
  
  data %>%
    filter(type == "aerenchyma",
           !is.na(neib1))%>%
    ggplot()+
    geom_polygon(aes(x,y, group = id_cell, fill = type),alpha = 0.6, colour = "white", data = rs1)+
    geom_polygon(aes(x,y, group = id_cell, fill = type), colour = "white")+
    geom_point(aes(x,y), colour = "red", data = data[!is.na(data$neib2),])+
    geom_point(aes(x,y), colour = "blue", data = noise)+
    coord_fixed()
  
  
  return(data%>%select(colnames(rs1)))
}

pv_ready <- function(rs1){
  data <- rs1
  
  data$neib1 = data$neib2 = data$neib3 <- NA
  data$neib1_type = data$neib2_type = data$neib3_type <- NA
  
  for(i in 1:nrow(data)){
    ne <- unique(data$id_cell[data$id_point == data$id_point[i] & data$id_cell != data$id_cell[i]])
    data$neib1[i] <- ne[1]
    data$neib1_type[i] <- as.character(data$type[data$id_cell == ne[1]][1])
    data$neib2[i] <- ne[2]
    data$neib2_type[i] <- as.character(data$type[data$id_cell == ne[2]][1])
    data$neib3[i] <- ne[3]
    data$neib3_type[i] <- as.character(data$type[data$id_cell == ne[3]][1])
  }
  
  junc_point <- data[!is.na(data$neib2) | (data$type == "epidermis" & !is.na(data$neib1)),]
  
  data %>%
    ggplot()+
    geom_polygon(aes(x,y,fill = type, group = id_cell), colour = "white")+
    geom_point(aes(x,y), colour = "red", data = junc_point)+
    geom_point(aes(x,y), colour = "blue", data = data[data$id_point %!in% junc_point$id_point,])+
    coord_fixed()+
    theme_classic()
  
  nodes <- data
  
  nodes<- nodes%>%
    dplyr::group_by(id_cell)%>%
    dplyr::filter(!duplicated(id_point))%>%
    dplyr::mutate(x1 = x,
                  y1 = y,
                  x2 = c(x[-1],x[1]),
                  y2 = c(y[-1],y[1]),
                  id_point2 = paste0(x2,";",y2),
                  w_length = sqrt((x1-x2)^2+(y1-y2)^2))
  
  on_go <- T
  while(on_go){
    
    last_point <- nodes%>%
      select(id_cell, last(starts_with("id_point")),
             last(starts_with("x")),
             last(starts_with("y"))) # take last colmn
    n <- colnames( last_point)
    n <- unique(parse_number(n[-1])) # get the incremental value
    colnames(last_point) <- c("id_cell","id_point", "x","y") # generic col names
    
    last_point <- last_point%>%
      dplyr::group_by(id_cell)%>%
      dplyr::mutate(xx = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(x[-1],x[1])),
                    yy = ifelse(id_point %in% junc_point$id_point | id_point == "NA;NA",NA,c(y[-1],y[1])),
                    id_pointxy = paste0(xx,";",yy))%>%
      ungroup()%>%
      select(xx,yy,id_pointxy)
    
    on_go <- length(which(last_point$id_pointxy != "NA;NA"))>0
    c_name <- paste0(t(c("x", "y", "id_point")),n+1)
    colnames(last_point) <- c(c_name)
    
    nodes <- cbind(nodes%>%
                     ungroup(),last_point)
  }
  
  nodes <- nodes%>%
    dplyr::filter(id_point %in% junc_point$id_point,
                  w_length > 0)# only one wall per junction point
  
  nodes%>%
    ggplot()+
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), alpha = 0.5)+
    geom_segment(aes(x = x2, xend = x3, y = y2, yend = y3), alpha = 0.5)+
    geom_segment(aes(x = x3, xend = x4, y = y3, yend = y4), alpha = 0.5)+
    coord_fixed()
  
  # at this points all walls are double
  
  nodus <- NULL
  more <- nodes
  k <- 2
  while(nrow(more) > 0){
    print(k)
    last_x <- paste0("x",k) # the last point has coord_x
    last_y <- paste0("y",k) # the last point has coord_y
    tag_x <- paste0("x",k+1)  # the next possible points
    
    tmp <- more[is.na(more[,tag_x]),] # if next point is "nan" then
    tmp <- tmp%>%
      # dplyr::group_by(id_cell) %>%
      dplyr::mutate(xx = tmp[,last_x])%>%
      dplyr::mutate(yy = tmp[,last_y])
    
    deto <- tmp%>%select(-x1,-y1,-c(last_x),-c(last_y)) # remove points
    h <- 1
    # dealing with first and last point
    tmp_h <- tmp%>%
      mutate(xh = ifelse(x > xx, x, xx)) %>%
      mutate(yh = ifelse(x > xx, y, yy)) %>%
      mutate(yh = ifelse(x == xx, ifelse (y > yy, yy, y), yh))%>%
      mutate(xlh = ifelse(x > xx, xx, x)) %>%
      mutate(ylh = ifelse(x > xx, yy, y)) %>%
      mutate(ylh = ifelse(x == xx, ifelse (y > yy, y, yy), ylh))%>%
      select(xh,yh, xlh, ylh)
    colnames(tmp_h) <- c(paste0("x",h),paste0("y",h), paste0("x",k),paste0("y",k))
    print(c(paste0("x",h),paste0("y",h), paste0("x",k),paste0("y",k)))
    deto <- cbind(deto, tmp_h)
    if(k >= 4){
      for(h in c(2:floor(k/2))){
        h_x <- paste0("x",h)
        lh_x <- paste0("x",k+1-h)
        h_y <- paste0("y",h)
        lh_y <- paste0("y",k+1-h)
        
        tmp_h <- tmp %>%
          mutate(hx = tmp[,h_x],
                 hy = tmp[,h_y],
                 lhx = tmp[,lh_x],
                 lhy = tmp[,lh_y])%>%
          mutate(xh = ifelse(x > xx, hx, lhx)) %>%
          mutate(yh = ifelse(x > xx, hy, lhy)) %>%
          mutate(xlh = ifelse(x > xx, lhx, hx)) %>%
          mutate(ylh = ifelse(x > xx, lhy, hy))%>%
          select(xh,yh,xlh,ylh)
        deto <- deto %>% select(-c(h_x,h_y,lh_x,lh_y))
        colnames(tmp_h) <- c(h_x,h_y,lh_x,lh_y)
        print(c(h_x,h_y,lh_x,lh_y))
        deto <- cbind(deto,tmp_h)
      }
    }
    nodus <- rbind(nodus,deto)
    k <- k + 1
    more <- more[!is.na(more[,tag_x]),] # reduce what left for the next loop
  }

  
  return(nodus)
  #' Plot root anatomy
  #'
  #' This function plot the results of a 2D root cross section anatomy simulation
  #' @param sim the simulation objkect, returned by 'create_anatomy.R'
  #' @param col Parameter to choose for the coloring of the cells. Accepted arguments are 'type', 'area', 'dist', 'id_cell', "segment' and 'angle'. Default = 'type'
  #' @param leg Display the legend; Default= TRUE
  #' @param apo_bar Display apolastic barrier when col = "segment". 1 endodermal casparian strip, 2 fully suberized endodermis, 3 fully suberized endodermis and an exodermal casparian strip, and 4 exodermis and endodermis are fully suberized.
  #' @keywords root
  #' @export
  #' @import viridis
  #' @examples
  #' create_anatomy("PATH_TO_XLM_FILE")
  #'
  
  
  plot_anatomy <- function(sim=NULL,
                           col = "type",
                           leg = T,
                           apo_bar = 2){
    
    pl <- ggplot(sim$nodes) +
      geom_polygon(aes_string("x", "y", group="id_cell", fill=col), colour="white") +
      theme_classic() +
      coord_fixed() +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    if(!col %in% c("type", "cell_group")){
      pl <- pl + scale_fill_viridis()
    }
    
    if(!leg){
      pl <- pl + theme(legend.position="none")
    }
    if(col == "segment"){
      pl <- ggplot()+
        geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = sim$nodes)+
        theme_classic()+
        coord_fixed()+
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank())
      
      if(apo_bar != 0){
        if(apo_bar %in% c(2,4)){
          if(apo_bar == 4){
            pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                    data = sim$nodes%>%
                                      filter(type == "exodermis"))
            
          }
          pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                  data = sim$nodes%>%
                                    filter(type == "endodermis"))
        }
        if(apo_bar %in% c(1,3)){
          nod_apo <- sim$nodes
          cx <- unique(nod_apo$mx[nod_apo$id_cell == 1])
          cy <- unique(nod_apo$my[nod_apo$id_cell == 1])
          nod_apo <- nod_apo%>%
            mutate(slo = atan((y2-y1)/(x2-x1)),
                   SLO = atan((y1-cy)/(x1-cy)),
                   d = (slo-SLO))
          
          if(apo_bar == 1){
            pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                                      filter(type == "endodermis"
                                             ,d < 0.4 & d > - 0.4))
          }else{
            pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                    data = sim$nodes%>%
                                      filter(type == "endodermis"))+
              geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                             filter(type == "exodermis"
                                    ,d < 0.4 & d > - 0.4))
          }
        }
      }
    }
    
    return(pl)
  }
  #' @title Write the root anatomy as an XML file.
  #'
  #' The structure of the XL file matches the one of CellSet
  #' @param sim The simulation results
  #' @param path The path where to save the xml file
  #' @keywords root
  #' @export
  #' @examples
  #' write_anatomy_xml(sim, path = "current_root.xml" )
  #'
  
  write_anatomy_xml <- function(sim = NULL, path = NULL){
    
    if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
    if(is.null(path)) warning("No path found to save the XML file")
    
    
    if(length(sim$walls$x3) > 0){
      nodal <- sim$walls_nodes
    }else{
      nodal <- sim$nodes
    }
    nodal <- nodal%>%
      select(-id_group)
    
    cellgroups <- data.frame(id_group = c(1, 2, 3, 3, 4, 5, 13, 16, 12, 11, 4, 4, 6, 13),
                             type = c("exodermis", "epidermis", "endodermis", "passage_cell",  "cortex",
                                      "stele", "xylem", "pericycle", "companion_cell", "phloem",
                                      "inter_cellular_space", "aerenchyma", "cambium", "metaxylem"))
    
    xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
    xml <- paste0(xml, '<granardata>\n')
    
    # Write the Metadata
    xml <- paste0(xml, '\t<metadata>\n')
    xml <- paste0(xml, '\t\t<parameters>\n')
    xml <- paste0(xml,paste0('\t\t\t<parameter io="',sim$output$io,'" ',
                             'name="',sim$output$name,'" ',
                             'type="',sim$output$type,'" ',
                             'value="',sim$output$value,'"/>\n', collapse = ""))
    xml <- paste0(xml, '\t\t</parameters>\n')
    xml <- paste0(xml, '\t</metadata>\n')
    
    # Write the cells information
    xml <- paste0(xml, '\t<cells count="',nrow(sim$cells),'">\n')
    
    nodes_data <- merge(nodal, cellgroups, by="type")
    
    temp_wall <- ddply(nodes_data, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="',
                                                                                   paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                   '"/>\n'))
    xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                              '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                              '\t\t</cell>\n', collapse=""))
    xml <- paste0(xml, '\t</cells>\n')
    
    
    # Write the walls information
    xml <- paste0(xml, '\t<walls count="',length(unique(nodes_data$id_wall)),'">\n')
    
    walls <- sim$walls%>%
      select(starts_with("x"), starts_with("y"))
    col_nam <- colnames(walls)
    
    substr1(col_nam[nchar(col_nam) == 2], 1.5) <- "0"
    col_nam <- paste0(substr1(col_nam, -1), "_", substr1(col_nam, 1))
    colnames(walls) <- col_nam
    
    sorted_name <- sort(col_nam)
    walls <- walls%>%select(sorted_name)
    
    N <- max(parse_number(col_nam))
    begin <- tibble(tag1 = '\t\t<wall id="',
                    id_wall = sim$walls$id_wall-1,
                    tag2 = '" group="0" edgewall="false" >\n\t\t\t<points>\n')
    middle <- tibble(tag_x1 = '\t\t\t\t<point x="',
                     x1 = round(walls[,sorted_name[1]],6),
                     tag_y1 = '" y="',
                     y1 = round(walls[,sorted_name[2]],6),
                     tag_end1 = '"/>\n')
    for(k in 2:N){
      h <- k*2-1 # odd number
      tmp_coord <- walls[,sorted_name[c(h,h+1)]]
      tmp_middle <- tibble(tag_x = '\t\t\t\t<point x="',
                           x = round(tmp_coord[,1],6),
                           tag_y = '" y="',
                           y = round(tmp_coord[,2],6),
                           tag_end = '"/>\n')
      tmp_col_name <- colnames(tmp_middle)
      colnames(tmp_middle) <- paste0(t(tmp_col_name), k)
      middle <- cbind(middle, tmp_middle)
    }
    taged_walls <- cbind(begin,middle)%>%
      mutate(tag_ending = '\t\t\t</points>\n\t\t</wall>\n')
    xml <- paste0(xml, paste0(t(taged_walls), collapse = ""))
    xml <- paste0(xml, '\t</walls>\n')
    xml <- str_remove_all(xml, '\t\t\t\t<point x=\"NA\" y=\"NA\"/>\n')
    
    # Write the cell group informations
    print(cellgroups)
    xml <- paste0(xml, '\t<groups>\n')
    xml <- paste0(xml, '\t\t<cellgroups>\n')
    for(i in c(1:nrow(cellgroups))){
      xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id_group[i],'" name="',cellgroups$type[i],'" />\n')
    }
    xml <- paste0(xml, '\t\t</cellgroups>\n')
    xml <- paste0(xml, '\t\t<wallgroups>\n')
    xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
    xml <- paste0(xml, '\t\t</wallgroups>\n')
    xml <- paste0(xml, '\t</groups>\n')
    
    xml <- paste0(xml, '</granardata>')
    
    if(!is.null(path)){
      cat(xml, file = path)
      return(TRUE)
    }else{
      return(xml)
    }
    
    
  }
  
  substr1 <- function(x,y) {
    z <- sapply(strsplit(as.character(x),''),function(w) paste(na.omit(w[y]),collapse=''))
    dim(z) <- dim(x)
    return(z) }
  
  `substr1<-` <- function(x,y,value) {
    names(y) <- c(value,rep('',length(y)-length(value)))
    z <- sapply(strsplit(as.character(x),''),function(w) {
      v <- seq(w)
      names(v) <- w
      paste(names(sort(c(y,v[setdiff(v,y)]))),collapse='') })
    dim(z) <- dim(x)
    return(z) }
  
}



#' Write the root anatomy as an XML file.
#'
#' The structure of the XL file matches the one of CellSet
#' @param sim The simulation results
#' @param path The path where to save the xml file
#' @keywords root
#' @export
#' @examples
#' write_anatomy_xml()
#'

write_anatomy_xml <- function(sim = NULL, path = NULL){
  
  if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")
  
  
  if(length(sim$walls$x3) > 0){
    nodal <- sim$walls_nodes
  }else{
    nodal <- sim$nodes
  }
  nodal <- nodal%>%
    select(-id_group)
  
  cellgroups <- data.frame(id_group = c(1, 2, 3, 3, 4, 5, 13, 16, 12, 11, 4, 4, 6),
                           type = c("exodermis", "epidermis", "endodermis", "passage_cell",  "cortex", "stele", "xylem", "pericycle", "companion_cell", "phloem", "inter_cellular_space", "aerenchyma", "cambium"))
  
  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granardata>\n')
  
  # Write the Metadata
  xml <- paste0(xml, '\t<metadata>\n')
  xml <- paste0(xml, '\t\t<parameters>\n')
  xml <- paste0(xml,paste0('\t\t\t<parameter io="',sim$output$io,'" ',
                           'name="',sim$output$name,'" ',
                           'type="',sim$output$type,'" ',
                           'value="',sim$output$value,'"/>\n', collapse = ""))
  xml <- paste0(xml, '\t\t</parameters>\n')
  xml <- paste0(xml, '\t</metadata>\n')
  
  # Write the cells information
  xml <- paste0(xml, '\t<cells count="',nrow(sim$cells),'">\n')
  
  nodes_data <- merge(nodal, cellgroups, by="type")
  
  temp_wall <- ddply(nodes_data, .(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="',
                                                                                 paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                 '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')
  
  
  # Write the walls information
  xml <- paste0(xml, '\t<walls count="',length(unique(nodes_data$id_wall)),'">\n')
  
  walls <- sim$walls%>%
    select(starts_with("x"), starts_with("y"))
  col_nam <- colnames(walls)
  
  substr1(col_nam[nchar(col_nam) == 2], 1.5) <- "0"
  col_nam <- paste0(substr1(col_nam, -1), "_", substr1(col_nam, 1))
  colnames(walls) <- col_nam
  
  sorted_name <- sort(col_nam)
  walls <- walls%>%select(sorted_name)
  
  N <- max(parse_number(col_nam))
  begin <- tibble(tag1 = '\t\t<wall id="',
                  id_wall = sim$walls$id_wall-1,
                  tag2 = '" group="0" edgewall="false" >\n\t\t\t<points>\n')
  middle <- tibble(tag_x1 = '\t\t\t\t<point x="',
                   x1 = round(walls[,sorted_name[1]],6),
                   tag_y1 = '" y="',
                   y1 = round(walls[,sorted_name[2]],6),
                   tag_end1 = '"/>\n')
  for(k in 2:N){
    h <- k*2-1 # odd number
    tmp_coord <- walls[,sorted_name[c(h,h+1)]]
    tmp_middle <- tibble(tag_x = '\t\t\t\t<point x="',
                         x = round(tmp_coord[,1],6),
                         tag_y = '" y="',
                         y = round(tmp_coord[,2],6),
                         tag_end = '"/>\n')
    tmp_col_name <- colnames(tmp_middle)
    colnames(tmp_middle) <- paste0(t(tmp_col_name), k)
    middle <- cbind(middle, tmp_middle)
  }
  taged_walls <- cbind(begin,middle)%>%
    mutate(tag_ending = '\t\t\t</points>\n\t\t</wall>\n')
  xml <- paste0(xml, paste0(t(taged_walls), collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  xml <- str_remove_all(xml, '\t\t\t\t<point x=\"NA\" y=\"NA\"/>\n')
  
  # Write the cell group informations
  print(cellgroups)
  xml <- paste0(xml, '\t<groups>\n')
  xml <- paste0(xml, '\t\t<cellgroups>\n')
  for(i in c(1:nrow(cellgroups))){
    xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id_group[i],'" name="',cellgroups$type[i],'" />\n')
  }
  xml <- paste0(xml, '\t\t</cellgroups>\n')
  xml <- paste0(xml, '\t\t<wallgroups>\n')
  xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
  xml <- paste0(xml, '\t\t</wallgroups>\n')
  xml <- paste0(xml, '\t</groups>\n')
  
  xml <- paste0(xml, '</granardata>')
  
  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }
  
  
}


substr1 <- function(x,y) {
  z <- sapply(strsplit(as.character(x),''),function(w) paste(na.omit(w[y]),collapse=''))
  dim(z) <- dim(x)
  return(z) }

`substr1<-` <- function(x,y,value) {
  names(y) <- c(value,rep('',length(y)-length(value)))
  z <- sapply(strsplit(as.character(x),''),function(w) {
    v <- seq(w)
    names(v) <- w
    paste(names(sort(c(y,v[setdiff(v,y)]))),collapse='') })
  dim(z) <- dim(x)
  return(z) }

ApoSymp <- function(path = "MECHA/Projects/granar/out/Root/Project_Test/baseline/Macro_prop_1,0.txt"){
  
  coef_width_symplast=4/5
  mpercm=0.01
  dpersec=1/3600/24
  
  mecha_out_0 <- readLines(path)
  pos_dist <- which(mecha_out_0 =="Radial distance from stele centre (microns): ")
  pos_STF_up <- which(mecha_out_0 =="Standard Transmembrane uptake Fractions (%): ")
  pos_STF_rel <- which(mecha_out_0 =="Standard Transmembrane release Fractions (%): ")
  pos_end <- which(mecha_out_0 == "Scenario 1 ")
  pos_UP <- which(mecha_out_0 =="Stele, cortex, and epidermis uptake distribution cm^3/d: ")
  pos_UN <- which(mecha_out_0 =="Stele, cortex, and epidermis release distribution cm^3/d: ")
  pos_QXyl <- which(mecha_out_0 =="Xylem uptake distribution cm^3/d: ")
  pos_QPhlo <- which(mecha_out_0 =="Phloem uptake distribution cm^3/d: ")
  pos_disc <- which(mecha_out_0 =="Number of radial discretization boxes: ")
  
  q_tot <- parse_number( mecha_out_0[grepl("q_tot", mecha_out_0)]) # cm2 d-1
  c_height <- parse_number(unlist(str_split( mecha_out_0[grepl("height:", mecha_out_0)], " ")))
  c_height <- c_height[!is.na(c_height)]
  
  Q_tot = (q_tot*c_height)
  
  
  disc <- as.numeric(unlist(str_split(mecha_out_0[pos_disc+1], " ")))
  disco <- tibble(nlayer = disc[!is.na(disc)][-1], type = c("stele", "pericycle", "endodermis", "cortex", "exodermis", "epidermis"))
  
  
  dist <- parse_number(mecha_out_0[seq(pos_dist+1,pos_STF_up-2,1) ])
  STF_up <- parse_number(mecha_out_0[seq(pos_STF_up+1,pos_STF_rel-2,1) ])
  STF_down <- parse_number(mecha_out_0[seq(pos_STF_rel+1,pos_end-3,1) ])
  
  UP1 <- parse_number(mecha_out_0[seq(pos_UP[1]+1,pos_UN[1]-2,1) ])
  UN1 <- parse_number(mecha_out_0[seq(pos_UN[1]+1,pos_QXyl[1]-2,1) ])
  Qxyl1 <- parse_number(mecha_out_0[seq(pos_QXyl[1]+1,pos_QPhlo[1]-2,1) ])
  
  flux <- tibble(dist = dist, STF_up, STF_down, UP1,UN1, Qxyl1)
  
  Type <- NULL
  for (k in 1:nrow(disco) ){
    tmp <- rep(disco$type[k], disco$nlayer[k])
    Type <- c(Type, tmp)
  }
  flux$type <- Type
  
  flux$passage <- 0
  flux$passage[flux$type == "endodermis"] <- c(0,1,1,0)
  flux$intern <- 0
  flux$intern[flux$type == "endodermis"] <- c("in","in","out","out")
  
  flux_passage <- flux
  flux <- flux[flux$passage == 0 , ]
  
  flux$UP1[flux$intern == "in"] <- sum(flux_passage$UP1[flux_passage$intern == "in"])
  flux$UN1[flux$intern == "in"] <- sum(flux_passage$UN1[flux_passage$intern == "in"])
  flux$UP1[flux$intern == "out"] <- sum(flux_passage$UP1[flux_passage$intern == "out"])
  flux$UN1[flux$intern == "out"] <- sum(flux_passage$UN1[flux_passage$intern == "out"])
  
  endo <- flux$dist[flux$type == "endodermis"][1]
  flux <- flux %>%# filter(!duplicated(dist))%>%
    mutate(r_coord_layer = dist - endo)
  
  Flux <- NULL
  for(h in 2:nrow(flux)){
    tmp_m <- flux[c(h-1,h),]
    if(unique(tmp_m$type) %in% c("exodermis", "endodermis") & length(unique(tmp_m$type)) == 1){
      tmp <- tmp_m%>%mutate(label = "gate")
    }else{
      
      mid_wall <- mean(tmp_m$r_coord_layer)
      memb = (mid_wall-tmp_m$r_coord_layer[1])*coef_width_symplast
      tmp_m$label <- "mid_cell"
      
      if(tmp_m$Qxyl1[2] != 0){
        Addon <- tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                    label = "mid_wall",
                                    Qxyl1 = 0)
      }else{
        Addon <- NULL
      }
      
      tmp <- rbind(tmp_m[1,], # mid cell
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_out"), # wall
                   Addon,
                   tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                      label = "mid_wall"), # mid wall
                   
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_out"), # wall
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[2,]
      )
    }
    Flux <- rbind(Flux, tmp)
    
  }
  
  Flux <- Flux %>%mutate(r_coord_layer = -r_coord_layer)%>%arrange (r_coord_layer)
  
  Flux$type[Flux$label == "mid_wall"] <- "cell_wall"
  Flux$id <- 1:nrow(Flux)
  
  Flux$Sympl <- NA
  Flux$Apo <- NA
  
  Flux$Sympl[1:2] <- Flux$UP1[1:2]
  Flux$Sympl[3] <- Flux$Sympl[2]+Flux$UN1[3]
  Flux <- rbind(Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = 0,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = Q_tot,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0, Apo = Q_tot),
                Flux)
  
  com <- Flux[grepl("membrane", Flux$label),]
  
  for(i in 3:nrow(com)){
    if(i %% 2 == 0){
      if(com$type[i] %in% c("exodermis", "endodemris") ){
        com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]+com$UP1[i]
      }else{
        if (length(com$label[i][grepl("out", com$label[i])])>0){
          com$Sympl[i] <- com$Sympl[i-1]+com$UP1[i]
        }else{
          com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]
        }
      }
    }else{
      com$Sympl[i] <- com$Sympl[i-1]
    }
    Flux$Sympl[Flux$id == com$id[i]] <- com$Sympl[i]
  }
  
  
  cxyl <- Flux[grepl("mid_wall", Flux$label),]
  cxyl$Apo[1] <- Q_tot
  
  for(i in 2:nrow(cxyl)){
    cxyl$Apo[i] <- cxyl$Apo[i-1]+cxyl$Qxyl1[i]
    Flux$Apo[Flux$id == cxyl$id[i]] <- cxyl$Apo[i]
  }
  
  return(Flux)
}

plot_water_flux <- function(Flux, apobar = 1){
  Q_tot = Flux$Apo[2]
  r_0 <- Flux$r_coord_layer[1]
  
  colo <- c("Apoplastic compartment" = "burlywood3", "Symplastic compartment"= "dodgerblue2")
  
  pl <- Flux%>%
    filter(!is.na(Sympl) )%>%
    ggplot()+
    geom_polygon(aes(r_coord_layer, Apo, fill = "Apoplastic compartment"), data = Flux%>%filter(!is.na(Apo) ))+
    geom_polygon(aes(r_coord_layer, Sympl, fill = "Symplastic compartment") )+
    geom_line(aes(r_coord_layer, Sympl))+
    geom_vline(aes(xintercept = r_coord_layer), alpha = 0.2, linetype = 2, data = Flux%>%filter(label == "mid_wall"))+
    geom_vline(aes(xintercept = r_coord_layer), alpha = 0.5, linetype = 2, data = Flux%>%filter(type == "soil"))+
    geom_hline(yintercept = Q_tot)+
    geom_hline(yintercept = 0)+
    xlim(r_0-2, 90)+
    theme_classic()+
    ylab("Partition of radial water flow rates \n between compartments [cm3/d]")+
    xlab("Distance from the endodermis layer [um]")+
    labs(fill = "Compartment")+
    scale_fill_manual(values = colo)
  
  if(apobar == 1){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")
  }
  if(apobar == 2){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(id %in% c(48,57)))
  }
  if(apobar == 3){
    pl <- pl + geom_vline(aes(xintercept = 0), alpha = 0.6, size = 2, color = "red")+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(id %in% c(48,57)))+
      geom_vline(aes(xintercept = r_coord_layer), alpha = 0.6, size = 2, color = "red",
                 data = Flux%>%filter(label == "mid_cell",
                                      type == "exodermis"))
    
  }
  
  
  print(pl)
}


ApoSymp2 <- function(path = "MECHA/Projects/granar/out/Root/Project_Test/baseline/Macro_prop_1,0.txt"){
  
  coef_width_symplast=4/5
  mpercm=0.01
  dpersec=1/3600/24
  
  mecha_out_0 <- readLines(path)
  pos_dist <- which(mecha_out_0 =="Radial distance from stele centre (microns): ")
  pos_STF_up <- which(mecha_out_0 =="Standard Transmembrane uptake Fractions (%): ")
  pos_STF_rel <- which(mecha_out_0 =="Standard Transmembrane release Fractions (%): ")
  
  pos_end <- length(mecha_out_0)
  #pos_UP <- which(mecha_out_0 =="Stele, cortex, and epidermis uptake distribution cm^3/d: ")
  #pos_UN <- which(mecha_out_0 =="Stele, cortex, and epidermis release distribution cm^3/d: ")
  #pos_QXyl <- which(mecha_out_0 =="Xylem uptake distribution cm^3/d: ")
  #pos_QPhlo <- which(mecha_out_0 =="Phloem uptake distribution cm^3/d: ")
  pos_disc <- which(mecha_out_0 =="Number of radial discretization boxes: ")
  
  q_tot <- 1
  c_height <- parse_number(unlist(str_split( mecha_out_0[grepl("height:", mecha_out_0)], " ")))
  c_height <- c_height[!is.na(c_height)]
  
  Q_tot = (q_tot*c_height)
  
  
  disc <- as.numeric(unlist(str_split(mecha_out_0[pos_disc+1], " ")))
  disco <- tibble(nlayer = disc[!is.na(disc)][-1], type = c("stele", "pericycle", "endodermis", "cortex", "exodermis", "epidermis"))
  
  
  dist <- parse_number(mecha_out_0[seq(pos_dist+1,pos_STF_up-2,1) ])
  STF_up <- parse_number(mecha_out_0[seq(pos_STF_up+1,pos_STF_rel-2,1) ])
  STF_down <- parse_number(mecha_out_0[seq(pos_STF_rel+1,pos_end,1) ])
  
  UP1 <- STF_up*Q_tot
  UN1 <- STF_down*Q_tot
  Qxyl1 <- 0
  
  flux <- tibble(dist = dist, STF_up, STF_down, UP1,UN1, Qxyl1)
  
  Type <- NULL
  for (k in 1:nrow(disco) ){
    tmp <- rep(disco$type[k], disco$nlayer[k])
    Type <- c(Type, tmp)
  }
  flux$type <- Type
  
  flux$passage <- 0
  flux$passage[flux$type == "endodermis"] <- c(0,1,1,0)
  flux$intern <- 0
  flux$intern[flux$type == "endodermis"] <- c("in","in","out","out")
  
  flux_passage <- flux
  flux <- flux[flux$passage == 0 , ]
  
  flux$UP1[flux$intern == "in"] <- sum(flux_passage$UP1[flux_passage$intern == "in"])
  flux$UN1[flux$intern == "in"] <- sum(flux_passage$UN1[flux_passage$intern == "in"])
  flux$UP1[flux$intern == "out"] <- sum(flux_passage$UP1[flux_passage$intern == "out"])
  flux$UN1[flux$intern == "out"] <- sum(flux_passage$UN1[flux_passage$intern == "out"])
  
  endo <- flux$dist[flux$type == "endodermis"][1]
  flux <- flux %>%# filter(!duplicated(dist))%>%
    mutate(r_coord_layer = dist - endo)
  
  Flux <- NULL
  for(h in 2:nrow(flux)){
    tmp_m <- flux[c(h-1,h),]
    if(unique(tmp_m$type) %in% c("exodermis", "endodermis") & length(unique(tmp_m$type)) == 1){
      tmp <- tmp_m%>%mutate(label = "gate")
    }else{
      
      mid_wall <- mean(tmp_m$r_coord_layer)
      memb = (mid_wall-tmp_m$r_coord_layer[1])*coef_width_symplast
      tmp_m$label <- "mid_cell"
      
      if(tmp_m$Qxyl1[2] != 0){
        Addon <- tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                    label = "mid_wall",
                                    Qxyl1 = 0)
      }else{
        Addon <- NULL
      }
      
      tmp <- rbind(tmp_m[1,], # mid cell
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[1,]%>%mutate(r_coord_layer = r_coord_layer+memb,
                                      label = "membrane_out"), # wall
                   Addon,
                   tmp_m[2,]%>%mutate(r_coord_layer = mid_wall,
                                      label = "mid_wall"), # mid wall
                   
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_out"), # wall
                   tmp_m[2,]%>%mutate(r_coord_layer = r_coord_layer-memb,
                                      label = "membrane_in"), # membrane
                   tmp_m[2,]
      )
    }
    Flux <- rbind(Flux, tmp)
    
  }
  
  Flux <- Flux %>%mutate(r_coord_layer = -r_coord_layer)%>%arrange (r_coord_layer)
  
  Flux$type[Flux$label == "mid_wall"] <- "cell_wall"
  Flux$id <- 1:nrow(Flux)
  
  Flux$Sympl <- NA
  Flux$Apo <- NA
  
  Flux$Sympl[1:2] <- Flux$UP1[1:2]
  Flux$Sympl[3] <- Flux$Sympl[2]+Flux$UN1[3]
  Flux <- rbind(Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = 0,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0,
                                  Apo = Q_tot,
                                  type = "soil",
                                  r_coord_layer = r_coord_layer-diff(Flux$r_coord_layer[c(1,3)])/5),
                Flux[1,]%>%mutate(Sympl = 0, Apo = Q_tot),
                Flux)
  
  com <- Flux[grepl("membrane", Flux$label),]
  
  for(i in 3:nrow(com)){
    if(i %% 2 == 0){
      if(com$type[i] %in% c("exodermis", "endodemris") ){
        com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]+com$UP1[i]
      }else{
        if (length(com$label[i][grepl("out", com$label[i])])>0){
          com$Sympl[i] <- com$Sympl[i-1]+com$UP1[i]
        }else{
          com$Sympl[i] <- com$Sympl[i-1]+com$UN1[i]
        }
      }
    }else{
      com$Sympl[i] <- com$Sympl[i-1]
    }
    Flux$Sympl[Flux$id == com$id[i]] <- com$Sympl[i]
  }
  
  
  cxyl <- Flux[grepl("mid_wall", Flux$label),]
  cxyl$Apo[1] <- Q_tot
  
  for(i in 2:nrow(cxyl)){
    cxyl$Apo[i] <- cxyl$Apo[i-1]+cxyl$Qxyl1[i]
    Flux$Apo[Flux$id == cxyl$id[i]] <- cxyl$Apo[i]
  }
  
  return(Flux)
}
#' Plot root anatomy
#'
#' This function plot the results of a 2D root cross section anatomy simulation
#' @param sim the simulation objkect, returned by 'create_anatomy.R'
#' @param col Parameter to choose for the coloring of the cells. Accepted arguments are 'type', 'area', 'dist', 'id_cell', "segment' and 'angle'. Default = 'type'
#' @param leg Display the legend; Default= TRUE
#' @param apo_bar Display apolastic barrier when col = "segment". 1 endodermal casparian strip, 2 fully suberized endodermis, 3 fully suberized endodermis and an exodermal casparian strip, and 4 exodermis and endodermis are fully suberized.
#' @keywords root
#' @export
#' @import viridis
#' @examples
#' create_anatomy("PATH_TO_XLM_FILE")
#'


plot_anatomy <- function(sim=NULL,
                         col = "type",
                         leg = T,
                         apo_bar = 2){
  
  pl <- ggplot(sim$nodes) +
    geom_polygon(aes_string("x", "y", group="id_cell", fill=col), colour="white") +
    theme_classic() +
    coord_fixed() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  if(!col %in% c("type", "cell_group")){
    pl <- pl + scale_fill_viridis()
  }
  
  if(!leg){
    pl <- pl + theme(legend.position="none")
  }
  if(col == "segment"){
    pl <- ggplot()+
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = sim$nodes)+
      theme_classic()+
      coord_fixed()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    if(apo_bar != 0){
      if(apo_bar %in% c(2,4)){
        if(apo_bar == 4){
          pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                  data = sim$nodes%>%
                                    filter(type == "exodermis"))
          
        }
        pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                data = sim$nodes%>%
                                  filter(type == "endodermis"))
      }
      if(apo_bar %in% c(1,3)){
        nod_apo <- sim$nodes
        cx <- unique(nod_apo$mx[nod_apo$id_cell == 1])
        cy <- unique(nod_apo$my[nod_apo$id_cell == 1])
        nod_apo <- nod_apo%>%
          mutate(slo = atan((y2-y1)/(x2-x1)),
                 SLO = atan((y1-cy)/(x1-cy)),
                 d = (slo-SLO))
        
        if(apo_bar == 1){
          pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                                    filter(type == "endodermis"
                                           ,d < 0.4 & d > - 0.4))
        }else{
          pl <- pl + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1, colour = "red",
                                  data = sim$nodes%>%
                                    filter(type == "endodermis"))+
            geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), size = 1.2, colour = "red", data = nod_apo%>%
                           filter(type == "exodermis"
                                  ,d < 0.4 & d > - 0.4))
        }
      }
    }
  }
  
  return(pl)
}


