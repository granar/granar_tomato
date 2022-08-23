#' @title Write the parameters of the param dataframe as an XML file.
#'
#' The structure of the XL file matches the default input files of GRANAR
#' @param param The param dataframe
#' @param path The path where to save the xml file
#' @keywords root
#' @export
#' @examples
#' write_param_xml(params, path = "params.xml" )
#'
write_param_xml2 <- function(params, path){
  
  if(is.null(params)) warning("No params found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")
  
  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granar>\n')
  
  
  for(i in unique(params$name)){
    
    xml <- paste0(xml, '<',i)
    for(j in params$type[params$name == i]){
      xml <- paste0(xml, ' ',j, '="',params$value[params$name == i & params$type == j],'"')
    }
    xml <- paste0(xml, '/>\n')
  }
  
  xml <- paste0(xml, '</granar>')
  
  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }
}
