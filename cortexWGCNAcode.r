## loadfile       
library(WGCNA)
options(stringsAsFactors = FALSE)
setwd("...")
load(file = "cortexdatagendereffectremoved.RData")

## determine thresholding power for analysis       
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red")
abline(h=0.80,col="red")



## find gene coexpression modules            
net = blockwiseModules(datExpr, power = 14,
                       TOMType = "signed", minModuleSize = 150,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3, maxBlockSize = 10000)



##dendrogram                              
table(net$colors) 
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]





## associate modules with traits
# Define numbers of genes and samples   
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels         
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# plot correlations and their p-values
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot   
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(traits),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.75,zlim = c(-1,1),main = paste("Module-trait relationships"),legendLabel = "Correlation")

##extract condition     
Condition = as.data.frame(traits$WT)
names(Condition) = "Condition"


## calculate                                        
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, Condition, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Condition), sep="")
names(GSPvalue) = paste("p.GS.", names(Condition), sep="")


## GS and MM correlation for sig modules      
module = "purple"
column = match(module, modNames)
moduleGenes = moduleColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for Genotype",main = paste("Module membership vs. gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
module = "green"
column = match(module, modNames)
moduleGenes = moduleColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for FASD",main = paste("Module membership vs. gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

## extract green genes 
green <- names(datExpr)[moduleColors=="green"]

module = "purple"
column = match(module, modNames)
moduleGenes = moduleColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for Genotype",main = paste("Module membership vs. gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


## extract purple genes   
purple <- names(datExpr)[moduleColors=="purple"]


##enrichGO （）                
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
universe <- rownames(vsddf)
gopurple <- enrichGO(gene = purple, OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',ont = "All",pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = universe)
gopurple2 <- clusterProfiler::simplify(gopurple)
gogreen <- enrichGO(gene = green, OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',ont = "All",pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = universe)

##try enrichGO() without universe  
gogreen <- enrichGO(gene = green, OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',ont = "All",pAdjustMethod = "BH", pvalueCutoff = 0.05)
gopurple <- enrichGO(gene = purple, OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL',ont = "All",pAdjustMethod = "BH", pvalueCutoff = 0.05)


##simplify                 
gopurple2 <- clusterProfiler::simplify(gopurple)
gogreen2 <- clusterProfiler::simplify(gogreen)

##visualize             
plotgreen <- pairwise_termsim(gogreen2)
treeplot(plotgreen)
plotpurple <- pairwise_termsim(gopurple2)
treeplot(plotpurple)
