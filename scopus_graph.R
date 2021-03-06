tos_scopus <- function(file) {
  
  source("readFiles.R")
  
  giant.component <- function(graph) {
    cl <- clusters(graph)
    induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
  }
  
  data_raw_scopus <- 
    readScopus(file) %>% 
    separate_rows(CR, 
                  sep = "; ") 
  
  data_raw_edgelist_scopus <- 
    data_raw_scopus %>% 
    mutate(lastname = sub("\\., .*", "", CR),
           lastname = sub(",", "", lastname),
           lastname = sub("\\.", "", lastname),# extracting lastnames
           year = str_extract(CR, "[0-9]+")) %>%   # extracting year
    filter(!grepl(pattern = "[():[:digit:]]", lastname),
           str_length(year) == 4) %>% 
    mutate(id_scopus = paste0(lastname, ", ", year, ","),
           SR = str_trim(SR)) %>% 
    select(SR, id_scopus)
  
  graph_raw <- 
    graph.data.frame(data_raw_edgelist_scopus, 
                            directed = TRUE) %>% 
    simplify()
  
  graph <- 
    delete.vertices(graph_raw, 
                    which(degree(graph_raw, mode = "in") == 1 & 
                            degree(graph_raw, mode = "out") == 0)) %>% 
    giant.component()
    
  network.metrics <- tibble(
    id = V(graph)$name,
    indegree = degree(graph, mode = "in"),
    outdegree = degree(graph, mode = "out"),
    bet = betweenness(graph))
    
    
  
}
