load("C:/Users/Proprietario/Desktop/OneDrive - Universit? degli Studi di Bari/age.WGCNA.paper/age.WGCNA.paper_moduleList_grch38.RData")
jaccard = function(x,y){ length(intersect(x,y))/length(union(x,y)) }
caudate = moduleList[2:4]
dlpfc = moduleList[7:10]
hippo = moduleList[12:15]

caudate.all = lapply(caudate,function(x) x$grey)
dlpfc.all = lapply(dlpfc,function(x) x$grey)
hippo.all = lapply(hippo,function(x) x$grey)

caudate.scz = lapply(caudate.all,function(x) intersect(x,PGC3.all.biotypes$kbp_200))
dlpfc.scz = lapply(dlpfc.all,function(x) intersect(x,PGC3.all.biotypes$kbp_200))
hippo.scz = lapply(hippo.all,function(x) intersect(x,PGC3.all.biotypes$kbp_200))


# signif(jaccard(caudate.scz$caudate__.1.25,caudate.scz$caudate__25.50),2)

caudate.null = lapply(1:2,function(x){
  
  universe = unique(unlist(caudate.all))
  module = caudate.all[[x]]
  if(x == 1) { to.test = caudate.all[[2]]; age = "juvenile-adult" } else { to.test = caudate.all[[3]]; age = "adult-older adult" }
  pb2 <- progress_estimated(10000)
  module.null = purrr::map_dbl(1:10000,~{
    pb2$tick()$print()
    null = sample(universe, size = length(module), replace = F)
    set.seed(.x)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("caudate",10000),
             age_window = rep(age,10000),
             type = rep("ALL",10000),
             JI = module.null)
})

caudate.null.scz = lapply(1:2,function(x){
  
  universe = unique(unlist(caudate.scz))
  module = caudate.scz[[x]]
  if(x == 1) { to.test = caudate.scz[[2]]; age = "juvenile-adult" } else { to.test = caudate.scz[[3]]; age = "adult-older adult" }
  module.null = purrr::map_dbl(1:10000,~{
    null = sample(universe, size = length(module), replace = F)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("caudate",10000),
             age_window = rep(age,10000),
             type = rep("SCZ",10000),
             JI = module.null)
})
Caudate.null = as.data.frame( rbind( Reduce(rbind,caudate.null), Reduce(rbind,caudate.null.scz) ) )

dlpfc.null = lapply(1:3,function(x){
  
  universe = unique(unlist(dlpfc.all))
  module = dlpfc.all[[x]]
  if(x == 1) { to.test = dlpfc.all[[2]]; age = "perinatal-juvenile" } else if (x == 2) { to.test = dlpfc.all[[3]]; age = "juvenile-adult" }
  else { to.test = dlpfc.all[[4]]; age = "adult-older adult"}
  module.null = purrr::map_dbl(1:10000,~{
    null = sample(universe, size = length(module), replace = F)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("dlpfc",10000),
             age_window = rep(age,10000),
             type = rep("ALL",10000),
             JI = module.null)
})

dlpfc.null.scz = lapply(1:3,function(x){
  
  universe = unique(unlist(dlpfc.scz))
  module = dlpfc.scz[[x]]
  if(x == 1) { to.test = dlpfc.scz[[2]]; age = "perinatal-juvenile" } else if (x == 2) { to.test = dlpfc.scz[[3]]; age = "juvenile-adult" }
  else { to.test = dlpfc.scz[[4]]; age = "adult-older adult"}
  module.null = purrr::map_dbl(1:10000,~{
    null = sample(universe, size = length(module), replace = F)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("dlpfc",10000),
             age_window = rep(age,10000),
             type = rep("SCZ",10000),
             JI = module.null)
})
Dlpfc.null = as.data.frame( rbind( Reduce(rbind,dlpfc.null), Reduce(rbind,dlpfc.null.scz) ) )

hippo.null = lapply(1:3,function(x){
  
  universe = unique(unlist(hippo.all))
  module = hippo.all[[x]]
  if(x == 1) { to.test = hippo.all[[2]]; age = "perinatal-juvenile" } else if (x == 2) { to.test = hippo.all[[3]]; age = "juvenile-adult" }
  else { to.test = hippo.all[[4]]; age = "adult-older adult"}
  module.null = purrr::map_dbl(1:10000,~{
    null = sample(universe, size = length(module), replace = F)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("hippocampus",10000),
             age_window = rep(age,10000),
             type = rep("ALL",10000),
             JI = module.null)
})

hippo.null.scz = lapply(1:3,function(x){
  
  universe = unique(unlist(hippo.scz))
  module = hippo.scz[[x]]
  if(x == 1) { to.test = hippo.scz[[2]]; age = "perinatal-juvenile" } else if (x == 2) { to.test = hippo.scz[[3]]; age = "juvenile-adult" }
  else { to.test = hippo.scz[[4]]; age = "adult-older adult"}
  module.null = purrr::map_dbl(1:10000,~{
    null = sample(universe, size = length(module), replace = F)
    signif(jaccard(x = null, y = to.test),3)
  })
  data.frame(tissue = rep("hippocampus",10000),
             age_window = rep(age,10000),
             type = rep("SCZ",10000),
             JI = module.null)
})
Hippo.null = as.data.frame( rbind( Reduce(rbind,hippo.null), Reduce(rbind,hippo.null.scz) ) )


