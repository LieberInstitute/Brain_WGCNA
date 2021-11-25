#### MAGMA enrich preprocess ####
MAGMA_enrich_preprocess = function(genemap,pathology,target.directory = paste0(getwd(),"/"),output.directory = paste0(getwd(),"/"), magma.directory = paste0(getwd(),"/")){
  
  require(future);require(furrr);require(purrr)
  geneMap_fun = function(genes, ensembl = "37")
  {
    ### function to obtain gene location coordinates
    require(biomaRt)
    if(ensembl == "37"){
      ensembl.v = useMart("ENSEMBL_MART_ENSEMBL", 
                          dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
    }else if(ensembl == "38"){
      ensembl.v = useMart("ENSEMBL_MART_ENSEMBL",  
                          dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
    }else stop("incorrect ensembl version")
    Attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","strand")
    Filters = "ensembl_gene_id"
    df = getBM(Attributes,Filters,values = genes,mart = ensembl.v)
    df = df[df$chromosome_name %in% c(1:22,"X","Y"),]
    df$chromosome_name = factor(df$chromosome_name,levels = c(1:22,"X","Y"))
    df = df[order(df$chromosome_name,df$start_position,df$end_position),]
    df = df[!duplicated(df$ensembl_gene_id),]
    df$strand = ifelse(df$strand < 0, "-","+")
    return(df)
  }
  
  
  ### Main

  ## define kbp windows ########
  kbp.windows = c("0", "20", "50", "100", "150", "200", "250", "500")
#  kbp.windows = c("35.10","100")
  
  if(grepl(paste0("annotation\\.step\\.genes\\.annot$"),list.files(output.directory))){ 
    print("Annotation file(s) already present")} else {
      print("Annotation file will be generated")
      
      ### create gene location file first
      gene_least_protein.coding = geneMap_fun(genemap$ensemblID[genemap$gene_type %in% "protein_coding"],ensembl = "38")
      write.table(gene_least_protein.coding,file = paste0(output.directory,"protein.coding.gene.loc"), row.names = F, col.names = F, quote = F, sep = "\t")
      
      gene_least_all.biotypes = geneMap_fun(genemap$ensemblID,ensembl = "38")
      write.table(gene_least_all.biotypes,file = paste0(output.directory,"all.biotypes.gene.loc"), row.names = F, col.names = F, quote = F, sep = "\t")
      
      ######## Annotation step ####
      ### run Magma
      biotype = "all.biotypes"
      map(kbp.windows,~{
        kbp = .x
        magma = "magma"
        snp_loc  = paste0(target.directory,"g1000_eur.bim") ### snpmap of 1000 genome
        gene_loc = paste0(output.directory,biotype,".gene.loc") ### geneMap created previously
        if(kbp == "35.10"){window = paste0("window=",gsub("\\.",",",kbp))} else window = paste0("window=",kbp)
        out = paste0(output.directory,biotype,"_",kbp,"kbp_annotation.step")
        cmd = paste0(magma.directory,magma," --annotate ",window," --snp-loc ",snp_loc," --gene-loc ",gene_loc," --out ",out)
        
        if (!paste0(biotype,"_",kbp,"kbp_annotation.step.genes.annot.txt") %in% list.files(output.directory)){
          ##### run analysis
          system(cmd)
        }
      })
    }
  
  ######## Gene analysis step ####
  # sumstats = grep("\\.RData|\\.BARI|\\.LIBD|MAGMA|H\\-MAGMA|clumped|PGC2|PGC2\\.clozuk|PGC\\.eur",
  #                 list.files(paste0(target.directory,pathology)),value = T,invert = T)
  sumstats = grep("\\.use$", list.files(paste0(target.directory,"/",pathology)),value = T)#; sumstats = gsub("\\.use","",sumstats)
  
  biotype = "all.biotypes"
  purrr::walk(kbp.windows,~{
    kbp = .x; print(.x)
    purrr::walk(sumstats,~{
      pgc = .x; print(.x)
      
      ### run Magma
      magma      = "magma"
      gene_annot = paste0(output.directory,biotype,"_",kbp,"kbp_annotation.step.genes.annot") ### output of annotation step
      pval       = paste0(target.directory,pathology,"/",pgc," ncol=N") ### summary statistic of PGC
      out        = paste0(output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis")
      final.merged.file = paste0(biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.out.txt")
      
      if (basename(final.merged.file) %in% list.files(paste0(output.directory,pathology))){
        print(paste0("Merged final file already present: ", final.merged.file))
      } else {
        
        future::plan(multiprocess, workers = 3, verbose = T)
        #future(sequential)
        dir.file = future_map_chr(1:22,~{
          chr = .x
          
          bfile = paste0(target.directory,"g1000_eur_splitted/g1000_eur.chr",chr) ### 1000 genome genotype file for LD estimation
          cmd = paste0(magma.directory,magma," --bfile ",bfile," --gene-annot ",gene_annot," --pval ",pval," --batch ",chr," chr --out ",out)
          
          if ((pathology == "RA" & chr == 6)|(pathology == "UC" & chr == 6)) {
            cmd = paste0(cmd, " --gene-settings adap-permp=10000,10")
          }
          
          ##### run analysis
          raw.file = paste0(out, ".batch", chr, "_chr.genes.raw")
          out.file = paste0(out, ".batch", chr, "_chr.genes.out")
          
          if (!(basename(raw.file) %in% list.files(paste0(output.directory,pathology)) & basename(paste0(out.file,".txt")) %in% list.files(paste0(output.directory,pathology)))){
            print("Splitted gene analysis file not present")
            system(cmd)
          } else {
            print(paste0("Splitted gene analysis file already present: ", out.file))
          }
          return(gsub(".genes.raw$","",raw.file))
        })
        future::plan("sequential")
        
        #Test output
        iwalk(dir.file %>% magrittr::set_names(strsplit2(.,"batch")[, 2]), ~{
          print(.y)
          raw_files = readLines(paste0(.x, ".genes.raw"))
          out_files = readLines(paste0(.x, ".genes.out.txt"))
          append. = F
          if (.y != "1_chr"){ 
            raw_files = raw_files[3:length(raw_files)]
            out_files = out_files[2:length(out_files)]
            append. = T
          }
          write(raw_files, file =  paste0(output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.raw"), append = append.)
          write(out_files, file =  paste0(output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.out"), append = append.)
        })
        # ### merge MAGMA batch output
        # cmd.1 = paste0(paste0("{ cat ",dir.file[1],".genes.raw;"),paste0(" sed '1,2d' ",dir.file[-1],collapse = ".genes.raw;"),".genes.raw; } > ",
        #                output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.raw")
        # cmd.2 = paste0(paste0("{ cat ",dir.file[1],".genes.out;"),paste0(" sed '1d' ",dir.file[-1],collapse = ".genes.out;"),".genes.out; } > ",
        #                output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.out")
        # 
        # system(cmd.1)
        # print(paste0("Writing merged file: ",output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.raw"))
        # system(cmd.2)
        # print(paste0("Writing merged file: ",output.directory,pathology,"/",biotype,"_",pgc,"_",kbp,"kbp_gene.analysis.genes.out"))
        
      }
    })
  })
}



