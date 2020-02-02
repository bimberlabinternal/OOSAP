#' @title plot STRINGdb networks
#'
#' @description Takes Differential Expression table and plots STRINGdb networks 
#' @param DEtable A DE table
#' @param numHits The num of mapped hits to plot
#' @param refSpeciesNum The dataset (see STRINGdb docs) to use as a reference; 9606=Human, ?=Rhesus Macaque
#' @return The PNGs of network plots
#' @keywords STRINGdb
#' @import STRINGdb
#' @export
#' @importFrom

##TODO: save pngs
  ##change db version

runSTRINGdb <- function(DEtable, numHits, refSpeciesNum){
  string_db <- STRINGdb$new(version="10", species=refSpeciesNum, score_threshold=0, input_directory="")
  
  DEtable.split <- split(DEtable, DEtable$cluster)
  
  for (i in as.vector(names(DEtable.split))){
    tryCatch({
      clusterTable <- DEtable.split[[i]]
      
      cluster.map <- string_db$map(clusterTable, "gene", removeUnmappedRows = FALSE)
      hits <- cluster.map$STRING_id[1:numHits]
      
      ##payload mechanism for upregulated vs downregulated genes:
      ##adds a color column for up vs downregulated genes
      cluster.color <- string_db$add_diff_exp_color(cluster.map, logFcColStr="avg_logFC")
      # post payload information to the STRING server
      payload_id <- string_db$post_payload(cluster.color$STRING_id, colors=cluster.color$color )
      string_db$plot_network(hits, payload_id=payload_id)
      
      ##clustering/community algorithms: ”fastgreedy”, ”walktrap”, ”spinglass”, ”edge.betweenness”.
      networkClustersList <- string_db$get_clusters(cluster.map$STRING_id[1:600], algorithm = "fastgreedy")
      par(mfrow=c(2,2))
      for(i in seq(1:length(networkClustersList))){
        string_db$plot_network(networkClustersList[[i]], payload_id=payload_id)
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
