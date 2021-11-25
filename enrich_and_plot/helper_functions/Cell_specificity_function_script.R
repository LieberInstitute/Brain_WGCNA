##################################################################################################
### CELL-SPECIFICITY geneset enrichment analysis with limma ######################################
##################################################################################################

Cell_Specificity_function = function(Net, network.name, genemap, specificity){
  
  library(limma)
  
  getResults <- function(specificity) {
    ### GSEA
    
    results = sapply(colnames(specificity), function(v){ 
      sapply(candidate.list, function(m){ 
        geneSetTest(index = toupper(rownames(specificity)) %in% toupper(m), statistics = specificity[,v], 
                    alternative = "up", type= "t", ranks.only = T, nsim=9999) 
      })
    })
    
    results[results==0] <- 1e-200
    
    nn = nrow(results)*ncol(results)
    #results.fdr  = sapply(colnames(results), function(x) p.adjust(results[,x],method = "fdr"       , n = nn))
    results.bonf = sapply(colnames(results), function(x) p.adjust(results[,x],method = "bonferroni", n = nn))
    
    colnames(results)      = paste0("pvals..",colnames(results))
    colnames(results.bonf) = paste0("bonf.." ,colnames(results.bonf))
    
    return(cbind(results,results.bonf))
  }
  
  #####################
  ####### MAIN ########          
  
  ### load network data
  Net = Net[!names(Net) %in% "grey"]
  candidate.list = sapply(Net,function(x){
    genemap$Symbol[genemap$ensemblID %in% x]
  })
  
  out = getResults(specificity)  
  
}

