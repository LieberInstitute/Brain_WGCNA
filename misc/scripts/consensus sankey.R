con_ensembl = c(
  "ENSG00000101204"
  ,"ENSG00000119522"
  ,"ENSG00000162104"
  ,"ENSG00000173175"
  ,"ENSG00000196588"
  ,"ENSG00000198010"
  ,"ENSG00000215012"
  ,"ENSG00000002746"
  ,"ENSG00000105649"
  ,"ENSG00000174437"
  ,"ENSG00000175931"
  ,"ENSG00000184156"
  ,"ENSG00000198929"
  ,"ENSG00000084731"
  ,"ENSG00000130226"
  ,"ENSG00000139767"
  ,"ENSG00000171435"
  ,"ENSG00000178235"
  ,"ENSG00000183715"
  ,"ENSG00000241973"
  ,"ENSG00000110427"
  ,"ENSG00000133958"
  ,"ENSG00000139182"
  ,"ENSG00000144406"
  ,"ENSG00000145934"
  ,"ENSG00000147724"
  ,"ENSG00000155093"
  ,"ENSG00000242732"
)

load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/enrichment.shorter.RData")

names(gm_AB)[grep("Gandal2018$"     ,names(gm_AB))] = "Gandal2018a"
names(gm_AB)[grep("Gandal2018PE$"   ,names(gm_AB))] = "Gandal2018b"
names(gm_AB)[grep("Gandal2018PE_cs$",names(gm_AB))] = "Gandal2018b_cs"
Hartl2021_gm <- readRDS("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/data/Hartl2021-gene-module-list.rds")
gm_AB$Hartl2021_BRNCTX =Hartl2021_gm$BRNCTX

conmodIDs =  function(gm, con) {
  imap_dfr(con, ~ {
    con_gene = ..1
    imap_dfr(gm, ~ {
      net = ..1
      net_name = ..2
      con_ind = map_lgl(net, ~ {
        if(sum(con_gene %in% ..1)>0) {
          out = TRUE
        } else {out = FALSE}
      })
      mod_con = names(net[con_ind])
      out = c(Network = ..2, Gene = con_gene, Module = mod_con)
    })
  })
}

gm_AB_DLPFC = gm_AB[c("dlpfc__.1.6","dlpfc__6.25","dlpfc__25.50","dlpfc__50.100")] %>% 
  setNames(c("DLPFC Perinatal","DLPFC Juvenile","DLPFC Adult","DLPFC Older Adult"))
gm_AB_con_modulIDs = conmodIDs(gm_AB_DLPFC, con_ensembl)
gm_AB_con_modulIDs <- pivot_wider(as.data.frame(gm_AB_con_modulIDs), names_from = "Network", values_from = "Module") 
gm_AB_con_modulIDs[2:5] = imap_dfc(gm_AB_con_modulIDs[2:5], ~ {
  paste0(..2, ">", ..1)
})
gm_AB_con_modulIDs = rbind(
  gm_AB_con_modulIDs[1:3] %>% setNames(.,c("hitgenes.intersect","source", "target")),
  gm_AB_con_modulIDs[c(1,3:4)] %>% setNames(.,c("hitgenes.intersect","source", "target")),
  gm_AB_con_modulIDs[c(1,4:5)] %>% setNames(.,c("hitgenes.intersect","source", "target"))
)
gm_AB_con_modulIDs = cbind(gm_AB_con_modulIDs,map_dfr(gm_AB_con_modulIDs$source, ~ {
  # edges.df = as.data.frame(edges.df)
  # print(..1)
  gene_sourceindex = unique(edges.df[edges.df$source ==..1 , "source.index" ])
  c(source = ..1,  gene_sourceindex)
})[2])
gm_AB_con_modulIDs = cbind(gm_AB_con_modulIDs,map_dfr(gm_AB_con_modulIDs$target, ~ {
  # edges.df = as.data.frame(edges.df)
  # print(..1)
  gene_sourceindex = unique(edges.df[edges.df$target ==..1 , "target.index" ])
  c(target = ..1, gene_sourceindex)
})[2] )
gm_AB_con_modulIDs$source.module = gsub(".*>", "", gm_AB_con_modulIDs$source)
gm_AB_con_modulIDs$target.module = gsub(".*>", "", gm_AB_con_modulIDs$target)
gm_AB_con_modulIDs$transitions = gm_AB_con_modulIDs$source
gm_AB_con_modulIDs$transitions[grepl("Perinatal",gm_AB_con_modulIDs$source)] = "t1"
gm_AB_con_modulIDs$transitions[grepl("Juvenile",gm_AB_con_modulIDs$source)] = "t2"
gm_AB_con_modulIDs$transitions[grepl("Adult",gm_AB_con_modulIDs$source)] = "t3"
gm_AB_con_modulIDs$network = "DLPFC"
gm_AB_con_modulIDs$hitgenes.intersect.count = 20
gm_AB_con_modulIDs = gm_AB_con_modulIDs[names(edges.df)]

new.names = read_csv("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/madhur data/old-to-new Network names all.csv") %>% tibble::deframe()
new.names[grepl("Gandal2018$"     ,names(new.names))] = "Gandal2018a"
new.names[grepl("Gandal2018PE$"   ,names(new.names))] = "Gandal2018b"
new.names[grepl("Gandal2018PE_cs$",names(new.names))] = "Gandal2018b_cs"
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))
new.names = c(new.names,c("Gandal2018a"="Gandal2018a","Gandal2018b"="Gandal2018b","Gandal2018b_cs"="Gandal2018b_cs"))
new.names = c(new.names,c("Hartl2021_BRNCTX"="Hartl2021_BRNCTX"))

load("C:/Users/Proprietario/OneDrive - Università degli Studi di Bari/R Projects/wgcna age/data/sankey_raw_data.RData")
##### Code for sankey plots (hitgenes) -----
make_sankeys = function(all.files0 = all.files, SCZ_risk_modules_df0 = SCZ_risk_modules_df, new.names0 = new.names){
  
  sankey = all.files0 %>% dplyr::select(network, enrichment, modules, module.length, hitgenes..kbp_200) %>% 
    mutate(new.network        = new.names0[network],
           new.network.module = paste0(new.network,">",modules),
           hitgenes..kbp_200  = hitgenes..kbp_200 %>% set_names(new.network.module))                   %>%
    filter(enrichment %in% "PGC3.all.biotypes" & grepl(" ",new.network))                               %>%
    split(.,.$new.network)
  
  sankey.hitgenes    = map(sankey,"hitgenes..kbp_200") %>% flatten
  
  sankey.pairs = list(
    DLPFC = list(
      t1 = crossing(source = sankey$`DLPFC Perinatal`$new.network.module, target = sankey$`DLPFC Juvenile`$new.network.module)  ,
      t2 = crossing(source = sankey$`DLPFC Juvenile`$new.network.module , target = sankey$`DLPFC Adult`$new.network.module)     ,
      t3 = crossing(source = sankey$`DLPFC Adult`$new.network.module    , target = sankey$`DLPFC Older Adult`$new.network.module)
    ),
    
    HP = list(
      t1 = crossing(source = sankey$`HP Perinatal`$new.network.module   , target = sankey$`HP Juvenile`$new.network.module)  ,
      t2 = crossing(source = sankey$`HP Juvenile`$new.network.module    , target = sankey$`HP Adult`$new.network.module)     ,
      t3 = crossing(source = sankey$`HP Adult`$new.network.module       , target = sankey$`HP Older Adult`$new.network.module)
    ),
    
    CN = list(
      t1 = crossing(source = sankey$`CN Juvenile`$new.network.module    , target = sankey$`CN Adult`$new.network.module)     ,
      t2 = crossing(source = sankey$`CN Adult`$new.network.module       , target = sankey$`CN Older Adult`$new.network.module)
    )
  ) %>% map_dfr(.id = "network",~ bind_rows(.x,.id = "transitions")) %>%
    rowwise()                                                        %>%
    mutate(hitgenes.intersect       = {list(intersect(sankey.hitgenes[[source]],sankey.hitgenes[[target]]))}) %>%
    ungroup()                                                                                                 %>%
    mutate(hitgenes.intersect.count = lengths(hitgenes.intersect))                                            %>%
    filter(hitgenes.intersect.count >0)
  
  sankey.nodes.list  = data.frame(node = unique(c(sankey.pairs$source,sankey.pairs$target))) %>% separate(node, into = c("new.network","module"),sep =">",remove = F)
  
  ##Plot sankeys for three tissues
  # walk(c("DLPFC"
  #        # ,"HP","CN"
  #        ),~{
    tissue = "DLPFC"
    
    nodes.df = sankey.nodes.list %>% filter(grepl(tissue, node)) %>%
      mutate(module2 = ifelse(node %in% SCZ_risk_modules_df0$new.network.module, toupper(module),"")) %>%
      .[c(grep(paste0(tissue," Perinatal"),.$node),grep(paste0(tissue," Juvenile"),.$node),grep(paste0(tissue," Adult"),.$node),grep(paste0(tissue," Older Adult"),.$node)),]
    
    nodes.df$module2[nodes.df$module2 == "BLACK "] = "fsfsf "
    
    edges.df = sankey.pairs %>% filter(network %in% tissue) %>%
      mutate(source.index  = match(source, nodes.df$node)-1 ,
             target.index  = match(target, nodes.df$node)-1,
             source.module = strsplit2(source,">")[,2],
             target.module = strsplit2(target,">")[,2],)
    
    library(sankeyD3)
    
    mod.colors = paste0('"',unique(tolower(nodes.df$module)),'"', collapse = ',')
    my_color = paste0("d3.scaleOrdinal() .domain([", mod.colors, "]) .range([", mod.colors, "])")
    
    gm_AB_con_modulIDs$hitgenes.intersect.count = 20000
    edges.df_withcon = rbind(edges.df,gm_AB_con_modulIDs)
    # edges.df =gm_AB_con_modulIDs
    
    sankey1 = sankeyD3::sankeyNetwork(Links          = edges.df_withcon, 
                                      Source         = "source.index", 
                                      Target         = "target.index", 
                                      Value          = "hitgenes.intersect.count", 
                                      # LinkGroup      = "source.module", 
                                      linkGradient   = T, 
                                      linkColor      = "darkcyan",
                                      
                                      Nodes          = nodes.df, 
                                      NodeGroup      = "module", 
                                      NodeID         = "module2",
                                      showNodeValues = F, 
                                      #NodeValue = "totalgenes.count",
                                      nodeCornerRadius = 5,
                                      nodeLabelMargin  = 6,
                                      nodeStrokeWidth  = 1, 
                                      units            = "genes", 
                                      #NodePosX = "col.x.index",
                                      #nodePadding = 13, #scaleNodeBreadthsByString = T, nodeShadow = F,
                                      nodeWidth = 20, curvature = 0.5, linkOpacity = 0.4,
                                      # NodeFontSize = 50,
                                      fontSize            = 50, 
                                      iterations          = 10, 
                                      colourScale         = my_color, 
                                      # title               = paste0("Consensus genes"), 
                                      margin              = list(top = 20, bottom = 40,left = 75, right = 55), 
                                      numberFormat        = ",d",
                                      align               = "center", 
                                      dragX               = T, 
                                      dragY               = T, 
                                      zoom                = T, 
                                      highlightChildLinks = T, 
                                      doubleclickTogglesChildren = F, 
                                      linkType            = "path1", 
                                      xAxisDomain         = as.list(unique(nodes.df$new.network)),
                                      yOrderComparator = JS("function(a,b) {if (a.name == 'grey'||a.name == 'GREY') {return  -1;}; return a.y - b.y; }"),
                                      xScalingFactor = .95
    )
    
    sankey1 = htmlwidgets::onRender(sankey1, 'function(el){
                                alltext   = el.querySelector("g.x.axis");
                                alltext.setAttribute("font-size","19");
                                    }')
    
    sankey1$x$links$test = paste0(edges.df$source.module," -> ",edges.df$hitgenes.intersect.count," genes -> ",edges.df$target.module)
    
    sankey1 <- htmlwidgets::onRender(
      sankey1,
      '
  function(el, x) {
  d3.selectAll(".link").select("*").text(function(d) { return d.test ; });
  }
  '
    )
    print(sankey1)
    sankey1 %>% saveNetwork(file = paste0("Consensusflow","__new.html"),selfcontained = T)
  # })
  
    sankey2 = htmlwidgets::onRender(sankey1, 'function(el){
                                nodestext = el.querySelectorAll("text.node-text");
                                nodevalue = el.querySelectorAll("text.node-number");
                                link = el.getElementsByClassName("link");
                                
                                //nodestext.forEach(d =>{
                                nodestext.forEach(function(d,i){
                                    ff1 = "0px";ff2 = "0px"
                                    if (d.textContent == d.textContent.toUpperCase()) {ff1 = "25px";ff2 = "0px";}
                                    d.style.fontSize = ff1;
                                    d.style.textShadow = "1px 1px 3px " + d.textContent;
                                    nodevalue[i].style.fontSize = ff2;
                                })
                                
                                //allnodes  = d3.selectAll(".node-number");
                                //allnodes.text(d => {return(d.newtest);});
                                
                                alledges  = d3.selectAll(".link").select("*");
                                alledges.text(d => {
                                ans = d.source.name + " -> "+d.target.name +" (" + d.value +" genes) \t" + d.test;
                                return(ans);
                                });
                                
                                alltext   = el.querySelector("g.x.axis");
                                alltext.setAttribute("font-size","22");
                                
                                
                                ////nodevalue.forEach(d => {
                                ////if (d.textContent !== d.textContent.toUpperCase()) {ff2 = "0px";}
                                ////d.style.fontSize = ff2;
                                ////})
                                    }')
    # sankey2 = htmlwidgets::onRender(sankey2, 'function(el){
    #                             var list1,list2, index;
    #                             list1 = el.querySelectorAll("text.node-text");
    #                             list2 = el.querySelectorAll("text.node-number");
    #                             for (index = 0; index < list1.length; ++index) {
    #                                  ff1 = "0px";ff2 = "0px";
    #                                  if (list1[index].textContent == list1[index].textContent.toUpperCase()) {ff1 = "21px";ff2 = "0px"}
    #                                  ss1 = "cursor: move; fill: inherit; font-size: " + ff1 + "; font-family: inherit; text-shadow: 1px 1px 3px " + list1[index].textContent + ";";
    #                                  ss2 = "cursor: move; fill: inherit; font-size: " + ff2 + "; font-family: inherit;";
    # 
    #                                  list1[index].setAttribute("style",ss1);
    #                                  list2[index].setAttribute("style",ss2);
    #                                  }
    #                             }')
    
    #sankey2 = htmlwidgets::onRender(sankey2, 'function(el) {el.querySelector("g.x.axis").setAttribute("font-size","19")}')
    
    #sankey2 = htmlwidgets::onRender(sankey2, '$.each( getElementsByTagName("text.node-text"), function( index, value ){value.style.textShadow = "3px 1px 5px white";});')
    #htmlwidgets::onStaticRenderComplete('$.each( document.getElementsByTagName("svg"), function( index, value ){value.setAttribute("viewBox", null);});')
    
    #sankey2$x$links$test = unlist(map(edge.table.hitgenes$hitgenes.intersection.symbol, ~paste0(.x, collapse = ", ")))
    #sankey1$x$nodes$test = unlist(map(edge.table.allgenes$allgenes.intersection, ~paste0(.x, collapse = ", ")))
    
    #   sankey2 <- htmlwidgets::onRender(
    #     sankey2,
    #     '
    # function(el, x) {
    # d3.selectAll(".link").select("*").text(function(d) { return d.source.name + " -> "+d.target.name +" (" + d.value +" genes) \t" + d.test ; });
    # }
    # '
    #   )
    
    sankey2 %>% saveNetwork(file = paste0("DLPFC","Consensus","__new.html"),selfcontained = F)  
    #browser()
    
    library(webshot)
    webshot2::webshot(url = paste0("DLPFC","Consensus","__new.html"),file =  paste0("DLPFC","Consensus","__new.png"),
                      cliprect = NULL, expand = -5,delay = 1, vwidth = 25000, vheight = 16000,
                      # eval =  "casper.then(function() {
                      #   $(window).load(function() {
                      #     //dom not only ready, but everything is loaded
                      #     this.wait(500)
                      #   });});"
    )
    
    
    
  rm(sankey, sankey.hitgenes, sankey.pairs, sankey.nodes.list) 
}

make_sankeys(all.files, SCZ_risk_modules_df, new.names) #Makes sankey plots for DLPFC, CN, HIPPO hitgenes
