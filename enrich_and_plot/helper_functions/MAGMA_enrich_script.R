#################################################################
#          Function for MAGMA enrichment processing             #
#                Gene Level Analysis steps                      #  
#################################################################


#### MAGMA enrich ####
MAGMA_enrich = function(Net,pathology, Net_name,genemap,target.directory = getwd(),output.directory,magma.directory = NULL){
  require(future);require(furrr);require(purrr)
  #print(target.directory)
  
  # sumstats = grep("\\.RData|\\.BARI|\\.LIBD|MAGMA|H\\-MAGMA|clumped|PGC2|PGC2\\.clozuk|PGC3\\.eur",
  #                         list.files(target.directory),value = T,invert = T)
  sumstats = grep("\\.use$", list.files(target.directory),value = T)#; sumstats = gsub("\\.use","",sumstats)
  
  ## define kbp windows ########
  kbp.windows = c("0", "20", "50", "100", "150", "200", "250", "500")
  kbp.windows = c("35.10","100")
  
  ######## Gene-level analysis step preprocess ####
  purrr::map(c("protein.coding","all.biotypes"),~{
    
    biotype = .x
    
    ### obtain genes in the network and subset them based on gene_type
    if(biotype == "protein.coding"){
      gene_least = unique(unlist(Net)) 
      gene_least = genemap$ensemblID[with(genemap,ensemblID %in% gene_least & gene_type %in% "protein_coding")]
      Net = sapply(names(Net),function(x)intersect(Net[[x]],gene_least),USE.NAMES = T)
      if(any(lengths(Net)<10)) {warning(paste("Network ", Net_name," has PC module with less than 10 genes"))}
    }
    
    ### create gene set info file (for each gene obtain specific component)
    comp = rep(names(Net), times = sapply(Net, length))
    geneList_new      = as.data.frame(unlist(Net), stringsAsFactors = F)
    geneList_new$comp = comp
    geneList_new      = geneList_new[,c(2,1)]
    #View(geneList_new); print(paste0(output.directory,Net_name,"_",biotype,".gene.set"))
    
    write.table(geneList_new, file = paste0(output.directory,Net_name,"_",biotype,".gene.set"), row.names = F, col.names = F, quote = F, sep = "\t")
    print(paste0("Writing done for: ",output.directory,Net_name,"_",biotype,".gene.set"))
    
    map(kbp.windows,~{
      kbp = .x
      
      ### run Magma
      purrr::map(sumstats,~{
        pgc = .x
        magma = "magma"
        gene_results = paste0(output.directory,"/all.biotypes","_",pgc,"_",kbp,"kbp_gene.analysis.genes.raw") ### output of gene analysis already performed ###Select all.biotype .raw file only and then subset based on PC geneset
        #rr_gene_results = readLines(gene_results)
        #View(rr_gene_results)
        
       
        gene_set = paste0(output.directory,Net_name,"_",biotype,".gene.set") ### gene set info file created previously
        rr_gene_set = read.table(gene_set)
        writeLines(unique(rr_gene_set[,2]),paste0(output.directory,paste0(Net_name,"_",biotype,"_genes.txt")))
        
        gene_include = paste0(output.directory,Net_name,"_",biotype,"_genes.txt") ### output of gene analysis already performed
        
        out = paste0(output.directory,Net_name,"_",biotype,"_",pgc,"_",kbp,"kbp_gene-level.analysis")
        #cmd = paste0(magma.directory,magma," --gene-results ",gene_results," --set-annot ",gene_set," col=2,1 --out ",out)
        cmd = paste0(magma.directory,"/",magma," --gene-results ",gene_results," --set-annot ",gene_set," col=2,1 --out ",out, " --settings gene-include=", gene_include)
    
        ##### run analysis
        if (target.directory == "C:/Users/mpariha1/Documents/SCZ/") {system(cmd)}
        
        
      })
    })
  })
  
  
  future::plan("sequential")
  print("I am here")
  print(target.directory)
  
  rm(list = setdiff(ls(),c("Net","Net_name","genemap","target.directory","output.directory","sumstats","kbp.windows")))
  
  ##### organize output file #####
  PGC_MAGMA_out = list()
  for(biotype in c("protein.coding","all.biotypes")){
    for (kbp in kbp.windows){
      for (pgc in sumstats){
        ### load output files
        
        out = paste0(output.directory,Net_name,"_",biotype,"_",pgc,"_",kbp,"kbp_gene-level.analysis.gsa.out.txt")
        print(paste0("Reading file:  ",out))
       
        PGC_out <- read.table(out,
                              header=TRUE, quote="\"", stringsAsFactors=FALSE)
        
        ### set pvalue = 1 for missing components
        miss_component = setdiff(names(Net),PGC_out$VARIABLE)
        if(length(miss_component)!=0){
          for (i in miss_component){
            miss_comp = PGC_out[1,]
            miss_comp$VARIABLE = i
            miss_comp$P = 1
            PGC_out = rbind(PGC_out,miss_comp)
          }
        }
        ### adjust with fdr and bonf correction
        PGC_out$P_fdr = p.adjust(PGC_out$P,method = "fdr")
        PGC_out$P_bonf = p.adjust(PGC_out$P,method = "bonferroni")
        ### order dataframe based on components name
        PGC_out$VARIABLE = factor(PGC_out$VARIABLE,levels = names(Net))
        PGC_out = PGC_out[order(PGC_out$VARIABLE),]
        PGC_MAGMA_out[[biotype]][[paste0(kbp,"kbp")]][[gsub(".use","",pgc)]] = PGC_out
      }
    }
  }
  ### save ###
  print(paste0(output.directory,Net_name,"_MAGMA.RData"))
  save(PGC_MAGMA_out,file = paste0(output.directory,Net_name,"_MAGMA.RData"))
  
  
  ##### plot ####
  require(grDevices)
  pdf(paste0(output.directory,Net_name,"_MAGMA.pdf"),width = 18, useDingbats = FALSE)
  par(mar = c(8,5,5,2), lwd = 96/72, ps = 8, las = 2)
  # print(warnings())
  # print(names(PGC_MAGMA_out))
  # View(PGC_MAGMA_out)
  for(biotype in names(PGC_MAGMA_out)){
    for(kbp in names(PGC_MAGMA_out[[biotype]])){
      for(pgc in names(PGC_MAGMA_out[[biotype]][[kbp]])){
        
        signif = which(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"P_fdr"] < 0.05)
        print(signif)
        if (length(signif) == 0){
          
          barplot(height = -log10(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"P_fdr"]),
                  col = "grey", names.arg = PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"VARIABLE"], space = 0,
                  main = paste0('MAGMA enrichment of ', Net_name,
                                '\n',biotype,' ',kbp,' ',pgc),
                  ylab = '-log10(fdr-adjusted p-value)', cex.lab = 3, cex.axis = 2, ylim = c(0,max(-log10(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"P_fdr"]))+1),las = 2)
          abline(h = -log10(0.05), lwd = 2, lty = 2, col = 'darkred')
          text(x = 5, y = -log10(.05) +.3, labels = 'FDR=0.05', col = 'darkred')
        } else {
          
          bp = barplot(height = -log10(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"P_fdr"]),
                       col = "grey", names.arg = PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"VARIABLE"], space = 0,
                       main = paste0('MAGMA enrichment of ', Net_name,
                                     '\n',biotype,' ',kbp,' ',pgc),
                       ylab = '-log10(fdr-adjusted p-value)', cex.lab = 3, cex.axis = 2, ylim = c(0,max(-log10(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][,"P_fdr"]))+2),las = 2)
          abline(h = -log10(0.05), lwd = 2, lty = 2, col = 'darkred')
          text(x = 5, y = -log10(.05) +.3, labels = 'FDR=0.05', col = 'darkred')
          text(x = bp[signif], y = -log10(PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][signif,"P_fdr"])+1.5,
               labels = PGC_MAGMA_out[[biotype]][[kbp]][[pgc]][signif,"NGENES"], cex = 1.5)
          
        }
      }
    }
  }
  
  dev.off()
  print(paste0(output.directory,Net_name,"_MAGMA.PDF"))
}
