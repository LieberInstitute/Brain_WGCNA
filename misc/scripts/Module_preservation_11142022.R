####################################################################################################
#                      Script to  test module preservation using WGCNA                             #
#            Preserving DLPFC Perinatal modules in DLPFC Prenatal networks                         #
#                                  Author: Eugenia Radulescu                                       #
####################################################################################################


####### Load the packages #######

library(WGCNA)

options(stringsAsFactors=FALSE)

####### Module preservation ##########


# - Preserve DLPFC -1_6 into DLPFC -1_0; mp1
# - Preserve DLPFC -1_0 into DLPFC 6_25; mp2
# - Preserve DLPFC -1_0 into Hippo -1_0; mp3
# - Preserve Hippo -1_6 into Hippo -1_0; mp4

# mp1:

load("C:/Users/eradule1/OneDrive - Johns Hopkins/Data_challenge_project/Giulio_paper/Module_preservation_11092022/new_age_groups.allData.RData")

dlpfc.O.gr.1_0<-assays(rse[[1]])

dlpfc.O.gr.1_0<-t(dlpfc.O.gr.1_0$ranknorm)

dlpfc.O.gr.1_6<-assays(rse[[2]])

dlpfc.O.gr.1_6<-t(dlpfc.O.gr.1_6$ranknorm)

dim(dlpfc.O.gr.1_0)
dim(dlpfc.O.gr.1_6)



# Create the multiExpr list:

nSets=2;
setLabels=c("dlpfc.O.gr.1_6","dlpfc.O.gr.1_0");
shortLabels=c("Gr_1_6","Gr_1_0");
multiExpr=vector(mode="list",length=nSets);

multiExpr[[1]]=list(data=dlpfc.O.gr.1_6);
names(multiExpr[[1]]$data)=names(dlpfc.O.gr.1_6);
rownames(multiExpr[[1]]$data)=dimnames(dlpfc.O.gr.1_6)[[1]]

multiExpr[[2]]=list(data=dlpfc.O.gr.1_0);
names(multiExpr[[2]]$data)=names(dlpfc.O.gr.1_0);
rownames(multiExpr[[2]]$data)=dimnames(dlpfc.O.gr.1_0)[[1]]

exprSize=checkSets(multiExpr)

# Create the multiColor list (1=module colors for reference; 2=module colors for test)

colors=dlpfc_gm_flat$`dlpfc.O.gr__-1_0`
colors_1_6=dlpfc_gm_flat$`dlpfc.O.gr__-1_6`

multiColor = list("Gr_1_6" = colors_1_6, "Gr_1_0" = colors);

identical(multiColor[[1]],dlpfc_gm_flat$`dlpfc.O.gr__-1_6`)
identical(colnames(multiExpr[[1]]$data),colnames(multiExpr[[2]]$data))


save(multiExpr,multiColor,file="Input1_11162022.RData")

# - Preserve DLPFC -1_0 into DLPFC 6_25; mp2

load("C:/Users/eradule1/OneDrive - Johns Hopkins/Data_challenge_project/Giulio_paper/Module_preservation_11092022/new_age_groups.allData.RData")


dlpfc.O.gr.6_25<-assays(rse[[3]])

dlpfc.O.gr.6_25<-t(dlpfc.O.gr.6_25$ranknorm)

dlpfc.O.gr.1_0<-assays(rse[[1]])

dlpfc.O.gr.1_0<-t(dlpfc.O.gr.1_0$ranknorm)


dim(dlpfc.O.gr.1_0)
dim(dlpfc.O.gr.6_25)



# Create the multiExpr list:

nSets=2;
setLabels=c("dlpfc.O.gr.1_0","dlpfc.O.gr.6_25");
shortLabels=c("Gr_1_0","Gr_6_25");
multiExpr=vector(mode="list",length=nSets);

multiExpr[[1]]=list(data=dlpfc.O.gr.1_0);
names(multiExpr[[1]]$data)=names(dlpfc.O.gr.1_0);
rownames(multiExpr[[1]]$data)=dimnames(dlpfc.O.gr.1_0)[[1]]

multiExpr[[2]]=list(data=dlpfc.O.gr.6_25);
names(multiExpr[[2]]$data)=names(dlpfc.O.gr.6_25);
rownames(multiExpr[[2]]$data)=dimnames(dlpfc.O.gr.6_25)[[1]]

exprSize=checkSets(multiExpr)

# Create the multiColor list (1=module colors for reference; 2=module colors for test)

colors=dlpfc_gm_flat$dlpfc.O.gr__6_25
colors_1_0=dlpfc_gm_flat$`dlpfc.O.gr__-1_0`



multiColor = list("Gr_1_0" = colors_1_0, "Gr_6_25" = colors);

identical(multiColor[[1]],dlpfc_gm_flat$`dlpfc.O.gr__-1_0`)
identical(colnames(multiExpr[[1]]$data),colnames(multiExpr[[2]]$data))


save(multiExpr,multiColor,file="Input2_11172022.RData")

# - Preserve DLPFC -1_0 into Hippo -1_0; mp3


load("C:/Users/eradule1/OneDrive - Johns Hopkins/Data_challenge_project/Giulio_paper/Module_preservation_11092022/new_age_groups.allData.RData")


dlpfc.O.gr.1_0<-assays(rse[[1]])

dlpfc.O.gr.1_0<-dlpfc.O.gr.1_0$ranknorm

Hippo.1_0<- assays(rse[[4]])

Hippo.1_0<-Hippo.1_0$ranknorm

dim(dlpfc.O.gr.1_0)
dim(Hippo.1_0)

rowID_DLPFC<-rownames(dlpfc.O.gr.1_0)
rowID_hippo<-rownames(Hippo.1_0)

commonProbesA = intersect (rownames(dlpfc.O.gr.1_0),rownames(Hippo.1_0))

DLPFC_index<-rowID_DLPFC %in% commonProbesA 

hippo_index<-rowID_hippo %in% commonProbesA


dlpfc_gm_flat<-dlpfc_gm_flat[DLPFC_index,]
hippo_gm_flat<-hippo_gm_flat[hippo_index,]

dlpfc.O.gr.1_0 = t(dlpfc.O.gr.1_0[DLPFC_index,])
Hippo.1_0 = t(Hippo.1_0[hippo_index,])



colors=hippo_gm_flat$`hippo.gr__-1_0`
colors_1_0=dlpfc_gm_flat$`dlpfc.O.gr__-1_0`

identical(commonProbesA,colnames(dlpfc.O.gr.1_0))
identical(commonProbesA,colnames(Hippo.1_0))
identical(dlpfc_gm_flat$GENECODE,colnames(dlpfc.O.gr.1_0))

# Create the multiExpr list:

nSets=2;
setLabels=c("dlpfc.O.gr.1_0","Hippo.1_0");
shortLabels=c("Gr_1_0","Gr_Hippo.1_0");
multiExpr=vector(mode="list",length=nSets);

multiExpr[[1]]=list(data=dlpfc.O.gr.1_0);
names(multiExpr[[1]]$data)=names(dlpfc.O.gr.1_0);
rownames(multiExpr[[1]]$data)=dimnames(dlpfc.O.gr.1_0)[[1]]

multiExpr[[2]]=list(data=Hippo.1_0);
names(multiExpr[[2]]$data)=names(Hippo.1_0);
rownames(multiExpr[[2]]$data)=dimnames(Hippo.1_0)[[1]]

exprSize=checkSets(multiExpr)



# Create the multiColor list (1=module colors for reference; 2=module colors for test)


multiColor = list("Gr_1_0" = colors_1_0, "Gr_Hippo.1_0" = colors);


save(multiExpr,multiColor,file="Input3_11212022.RData")

# - Preserve Hippo -1_6 into Hippo -1_0;

load("C:/Users/eradule1/OneDrive - Johns Hopkins/Data_challenge_project/Giulio_paper/Module_preservation_11092022/new_age_groups.allData.RData")


Hippo.1_6<- assays(rse[[5]])

Hippo.1_6<-Hippo.1_6$ranknorm

Hippo.1_0<- assays(rse[[4]])

Hippo.1_0<-Hippo.1_0$ranknorm

dim(Hippo.1_6)
dim(Hippo.1_0)


Hippo.1_6 = t(Hippo.1_6)
Hippo.1_0 = t(Hippo.1_0)



colors=hippo_gm_flat$`hippo.gr__-1_0`
colors_1_6=hippo_gm_flat$`hippo.gr__-1_6`


# Create the multiExpr list:

nSets=2;
setLabels=c("Hippo.1_6","Hippo.1_0");
shortLabels=c("Hippo.1_6","Hippo.1_0");
multiExpr=vector(mode="list",length=nSets);

multiExpr[[1]]=list(data=Hippo.1_6);
names(multiExpr[[1]]$data)=names(Hippo.1_6);
rownames(multiExpr[[1]]$data)=dimnames(Hippo.1_6)[[1]]

multiExpr[[2]]=list(data=Hippo.1_0);
names(multiExpr[[2]]$data)=names(Hippo.1_0);
rownames(multiExpr[[2]]$data)=dimnames(Hippo.1_0)[[1]]

exprSize=checkSets(multiExpr)



# Create the multiColor list (1=module colors for reference; 2=module colors for test)


multiColor = list("Hippo_1_6" = colors_1_6, "Hippo.1_0" = colors);


save(multiExpr,multiColor,file="Input4_11262022.RData")

##### Calculate the preservation statistics ##########

# Load the inputs;

mp<-c("mp1","mp2","mp3","mp4")

# Apply modulePreservation():

mp<-modulePreservation(
  multiData=multiExpr,
  multiColor=multiColor,
  multiWeights = NULL,
  dataIsExpr = TRUE,
  networkType = "signed hybrid", 
  corFnc = "cor",
  corOptions = "use = 'p'",
  referenceNetworks = 1, 
  testNetworks = NULL,
  nPermutations = 100, 
  includekMEallInSummary = FALSE,
  restrictSummaryForGeneralNetworks = TRUE,
  calculateQvalue = FALSE,
  randomSeed = 12345, 
  maxGoldModuleSize = 1000, 
  maxModuleSize = 1000, 
  quickCor = 1, 
  ccTupletSize = 2, 
  calculateCor.kIMall = FALSE,
  calculateClusterCoeff = FALSE,
  useInterpolation = FALSE, 
  checkData = TRUE, 
  greyName = NULL, 
  goldName = NULL,
  savePermutedStatistics = TRUE, 
  loadPermutedStatistics = FALSE, 
  permutedStatisticsFile = "permutedStats-actualModules.RData", 
  plotInterpolation = FALSE, 
  discardInvalidOutput = TRUE,
  parallelCalculation = FALSE,
  verbose = 3, indent = 0)


###### Create plots #######



stats = mp$preservation$Z$ref.Set_1$inColumnsAlsoPresentIn.Set_2
#stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2

stats[order(-stats[,2]),c(1:2)]

# Analyze and plot the preservation results

ref=1;
test=2;

statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

## We look at the main output: the preservation medianRank and Zsummary statistics.

## Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

## Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

## leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
## Text labels for points
text = modColors[plotMods];

#sizeGrWindow(8,7);
par(mfrow=c(1,2))
#par(mar=c(5, 7.4, 2.7, 1)+0.3);
plotData1=mp$preservation$observed[[ref]][[test]][, 2]
Lim=c(1:29)
mains1 = "Preservation Median rank"
min = min(Lim, na.rm = TRUE);
max = max(Lim, na.rm = TRUE);
{if (min > -max/10) min = -max/10
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))}

par(mar = c(4.5,4.5,2.5,1))
#par(mar=c(9,9,2.5,1))
#par(mar=c(5, 7.4, 2.7, 1)+0.3);
plot(moduleSizes[plotMods], plotData1[plotMods], col = 1, bg = modColors[plotMods], pch = 21,
     main = mains1,
     cex = 1.5,
     ylab = mains1, xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000),cex.lab = 1.2, cex.axis = 1, cex.main =1.2)

plotData2=mp$preservation$Z[[ref]][[test]][,2]
mains2 = "Preservation Zsummary";
max=max(plotData2)
min=min(plotData2)
if (min > -max/10) min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))

plot(moduleSizes[plotMods], plotData2[plotMods], col = 1, bg = modColors[plotMods], pch = 21,
     main = mains2,
     cex = 1.5,
     ylab = mains2, xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000),cex.lab = 1.2, cex.axis = 1, cex.main =1.2)

labelPoints(moduleSizes[plotMods],plotData2[plotMods],text, cex = 0.7, offs = 0.08);
abline(h=0)
abline(h=2, col = "red", lty = 2)
abline(h=7,col="red",lty=2)
abline(h=10, col = "blue", lty = 2)

