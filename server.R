suppressMessages(library (shiny))

source(paste(getwd(),'global.R',sep="/"))
options(shiny.maxRequestSize = 30*1024^2)
server <- function(input, output,session){

##### data input
dataInput <- reactiveValues(
	dataCount = read.delim(paste(getwd(),"/data/public_data/raw_counts.csv",sep="/"), header=T, row.names=1),
	dataMeta = read.delim(paste(getwd(),"/data/public_data/sample_info.txt",sep="/"), header=T) 
    )

myValues <- eventReactive(input$upload, {
	withProgress(message = 'Submit:', value = 0,{
	incProgress(.1, detail = paste("Uploading dataset"))
    if (is.null(input$countFile) & is.null(input$metaTab) ) {
      dataInput$dataCount <- read.delim(paste(getwd(),"/data/public_data/raw_counts.csv",sep="/"), header=T, sep=",", row.names=1)
      dataInput$dataMeta <- read.delim(paste(getwd(),"/data/public_data/sample_info.txt",sep=""), header=T, sep="\t")
    } else if (!is.null(input$countFile) & is.null(input$metaTab) ) {
      dataInput$dataCount <- read.delim(paste(getwd(),"/data/public_data/raw_counts.csv",sep="/"), header=T, sep=",", row.names=1)
      dataInput$dataMeta <- NULL
    } else if (is.null(input$countFile) & !is.null(input$metaTab) ) {
      dataInput$dataCount <- NULL
      dataInput$dataMeta <- read.delim(paste(getwd(),"/data/public_data/sample_info.txt",sep=""), header=T, sep="\t")
    } else if (!is.null(input$countFile) & !is.null(input$metaTab)) {
      dataInput$dataCount <- read.delim(input$countFile$datapath, header=T,sep=',', row.names=1)
      dataInput$dataMeta <- read.delim(input$metaTab$datapath, sep='\t',header=T)
    }   
    list(meta=(dataInput$dataMeta), count=(dataInput$dataCount) )
  })})
  
#### ENter geo ID
geovalues <- eventReactive(input$uploadGEO, {
	withProgress(message = 'Submit:', value = 0,{
	incProgress(.1, detail = paste("Fetching dataset"))
	
	gset <- getGEO(input$geodataset,GSEMatrix=TRUE, AnnotGPL=TRUE)
	if (length(gset) > 1) idx <- grep(gset[[attr(gset, "names")]]@annotation, attr(gset, "names")) else idx <- 1
	gset <- gset[[idx]]
	gset
	})
	})	
	
### Enter SRP ID

SRPvalues <- eventReactive(input$uploadSRP,{
withProgress(message = 'Submit:', value = 0,{
incProgress(.1, detail = paste("Fetching dataset"))
url <- recount::download_study(input$srpdataset)
load(file.path(input$srpdataset, "rse_gene.Rdata"))
rse <- recount::scale_counts(rse_gene)
y <- SE2DGEList(rse)
list(y=y,rse_gene=rse_gene)
})
})
	
 ####counts 
  
  datareactive <- reactive ({              
    org.counts <- dataInput$dataCount
    metadata <- dataInput$dataMeta
	group <- metadata[,2]
	samples <- metadata[,1]
	org.count <- DGEList(counts=org.counts, group=group, samples= samples,
						lib.size = colSums(org.counts),
						norm.factors = calcNormFactors(org.counts),)
	if (input$geneinfo == 'SYMBOL'){
	Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(org.count$counts),
                                columns=c("ENTREZID"),
                                keytype="SYMBOL",
                                multiVals="first")
	}else if (input$geneinfo == 'ENSEMBL'){
	Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(org.count$counts),
                                columns=c("ENTREZID","SYMBOL"),
                                keytype="ENSEMBL",
                                multiVals="first")
	}else if (input$geneinfo == 'ENTREZID'){
	Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(org.count$counts),
                                columns=c("SYMBOL"),
                                keytype="ENTREZID",
                                multiVals="first")
	}
  Hs_ann <- Hs_ann[!duplicated(Hs_ann[,1]),]
  org.count$genes <- Hs_ann
    org.count
  })
  
  
 
#### Data summary

output$sampleInfo <- renderTable(
	{
	meta <- myValues()$meta
	meta
	})

#### Data Summary of GEO ID

sampleGeo <- reactive({
if (input$uploadGEO){
	meta <- pData(phenoData(geovalues()))[,c(1,2,8,grep(":ch1",names(pData(phenoData(geovalues())))))]
names(meta)= gsub(":ch1", "", names(meta))
names(meta)= gsub("_ch1", "", names(meta))
as.data.frame(meta, rownames=T)
}
})


output$sampleInfo2 <- renderTable({	
sampleGeo()
})

#### Data Summary of SRP ID

sampleSrp <- reactive({
if (input$uploadSRP){
	geochar <- lapply(split(colData(SRPvalues()$rse_gene), seq_len(nrow(colData(SRPvalues()$rse_gene)))), geo_characteristics)
	geochar1 =dplyr::bind_rows(geochar)
	geochar1
}
})

output$sampleInfo3 <- renderTable({	
sampleSrp()
})

sampleTableSRP <- reactive({
sampleTable1 <- data.frame(condition=sampleSrp()[input$colfactor_SRP])
colnames(sampleTable1)='condition'
group <- factor(sampleTable1$condition)
groups <- make.names(c(gsub(" ", "_", levels(group))))
levels(group) <- groups
y2= SRPvalues()$y
y2$group=group
list(sampleTable1=sampleTable1, group=y2$group, groups=groups)
})

####Visualization of GEO data

geovisualize <- reactive({

if (input$groupassignmentgeo){
geofactor1 <- unlist(strsplit(input$factor1geo, ","))
geofactor2 <- unlist(strsplit(input$factor2geo, ","))
grpwords= c(geofactor1, geofactor2)
pattern <- paste(grpwords, collapse = "|")
sml= c(rep(input$namgeo1, length(grep(input$factor1geo,sampleGeo()[[input$colfactor]]))),rep(input$namgeo2, length(grep(input$factor2geo,sampleGeo()[[input$colfactor]]))))
sel=grep(pattern,sampleGeo()[[input$colfactor]])
groups <- make.names(c(input$namgeo1,input$namgeo2))
}else if (input$groupassignment){
geofactor1 <- unlist(strsplit(input$factor1, ","))
geofactor2 <- unlist(strsplit(input$factor2, ","))

sml <-c(rep(input$nam1, length(geofactor1)),rep(input$nam2, length(geofactor2)))
#sel = which(sampleGeo()[[input$colfactor]]==input$factor1|sampleGeo()[[input$colfactor]]==input$factor2)
sel = which(sampleGeo()[["geo_accession"]] %in% geofactor1 | sampleGeo()[["geo_accession"]]  %in% geofactor2)
groups <- make.names(c(input$nam1,input$nam2))
}
#sml <- sml[sel]
geovalues2 <- geovalues()[ ,sel]
gs <- factor(sml)
print (gs)
levels(gs) <- groups
#gs=factor(sampleGeo()[sel, "geo_accession"])
# groups <- make.names(c(gsub(" ", "_", levels(gs))))

geovalues2$group <- gs
ex <- exprs(geovalues2)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(geovalues2) <- log2(ex) }
ex <- exprs(geovalues2)
ord <- order(gs)  # order samples by group
list(ex=ex, ord=ord,gs=gs, groups=groups, geovalues2=geovalues2)
})

boxplotgeo <- eventReactive(c(input$groupassignment,input$groupassignmentgeo),ignoreInit = T,{
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title1 <- paste (input$geodataset, "/", annotation(geovalues()), " value distribution", sep ="")
boxplot(geovisualize()$ex, outline=FALSE, main=title1,las=2, col=geovisualize()$gs)
legend("topleft", geovisualize()$groups, fill=palette(), bty="n")
})

output$boxplot1 <-  renderPlot({
boxplotgeo()
})

expressiondensityplotreactive <- eventReactive(c(input$groupassignment,input$groupassignmentgeo),ignoreInit = T,{
title1 <- paste (input$geodataset, "/", annotation(geovalues()), " value distribution", sep ="")
plotDensities(geovisualize()$ex, group=geovisualize()$gs, main=title1, legend ="topright")
})

output$expressiondensityplot <- renderPlot({
expressiondensityplotreactive()
})

# umapreactive <- eventReactive(c(input$groupassignment,input$groupassignmentgeo),ignoreInit = T,{
# if (input$groupassignment | input$groupassignmentgeo)
# ex <- na.omit(geovisualize()$ex) # eliminate rows with NAs
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
# par(mar=c(3,3,2,6), xpd=TRUE)
# plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=geovisualize()$gs, pch=20, cex=1.5)
# legend("topleft", inset=c(-0.15,0), legend=levels(geovisualize()$gs), pch=20,
# col=1:nlevels(geovisualize()$gs), title="Group", pt.cex=1.5)
  # # point labels without overlaps
# pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# })

# output$umapplot <- renderPlot({
# umapreactive()
# })

filteredExprSRP <- reactive({
unfilteredExpr <- edgeR::cpm(SRPvalues()$y, log=T)
keep = filterByExpr(SRPvalues()$y, group=sampleTableSRP()$group) ### define the groups
SRPvalues2 <- SRPvalues()$y[keep,]

SRPvalues2 <- calcNormFactors(SRPvalues2)
# Plot the density of filtered gene expression for all samples within groups
filteredExpr <- edgeR::cpm(SRPvalues2, log=T)
list(filteredExpr=filteredExpr,SRPvalues2=SRPvalues2)
})

output$boxplot2 <- renderPlot({
if (input$groupassignment2)
boxplot(filteredExprSRP()$filteredExpr, las=2, col=as.numeric(factor(sampleTableSRP()$sampleTable1$condition)), main="")
})

output$expressiondensityplot2 <- renderPlot({
if(input$groupassignment2)
plotDensities(filteredExprSRP()$filteredExpr, group=sampleTableSRP()$sampleTable1$condition)
})


output$MDSplot2 <- renderPlot({
if(input$groupassignment2)
limma::plotMDS(filteredExprSRP()$filteredExpr, labels=sampleTableSRP()$sampleTable1$condition, col=as.numeric(factor(sampleTableSRP()$sampleTable1$condition)))
})

#### Library size of uploaded data	  
 output$orgLibsizeNormfactor <- renderTable({
        tab <- datareactive()$samples[,-1]
        #colnames(tab) <- c("Library sizes", "Normalization factors")
        tab
  })
  
  
output$boxplotunfiltered <- renderPlot({ 
   Group= datareactive()$samples$group
   unfilteredExpr <- edgeR::cpm(datareactive(), log=T)
   #myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
   par(mar=c(6,5,2,2))
   boxplot(unfilteredExpr, las=2, col=as.numeric(Group), main="")
   
    
  })
 
  output$plotdensityunfiltered <- renderPlot({ 
  meta <- myValues()$meta
    count <- myValues()$count
   unfilteredExpr <- edgeR::cpm(datareactive(), log=T)
   myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
   plotDensities(unfilteredExpr, group=factor(meta[,2]), col=myPalette[1:5], main="")   
  })
  
  
  ### Filter by EdgeR (uploaded data)
  
  FilterReactive <- reactive({    
    if (input$Filter){
	withProgress(message = 'Filter:', value = 0,{
	incProgress(.1, detail = paste("FIltering data"))
    keep = filterByExpr(datareactive(),min.count = 10, min.total.count = 15)
	filteredcounts <- datareactive()[keep,]
	filteredcounts <- calcNormFactors(filteredcounts)
	filteredcounts
    })
  }})
  
  ####Plot filtered data (uploaded data)
  
  output$sampleBoxplot <- renderPlot({ 
  Group = datareactive()$samples$group
  if (input$Filter){
   filteredExpr <- edgeR::cpm(FilterReactive(), log=T)
   #myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
   par(mar=c(6,5,2,2))
   boxplot(filteredExpr, las=2, col=as.numeric(Group), main="")
   }})
  
  
  output$plotdensities <- renderPlot({ 
  meta <- myValues()$meta
  
  if (input$Filter){
   filteredExpr <- edgeR::cpm(FilterReactive(),log=T)
   myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
   plotDensities(filteredExpr, group=factor(meta[,2]), col=myPalette[1:5], main="")
   }})
  
  output$mdsplot <- renderPlot({ 
  Group= datareactive()$samples$group
  if (input$Filter){
   filteredExpr <- edgeR::cpm(FilterReactive(),log=T)
   limma::plotMDS(filteredExpr, labels=Group, col=as.numeric(Group))
  }})
  
 
### Clustering of samples in GEO 

clusterreactive <- eventReactive(input$cluster,{
#words <- unlist(strsplit(input$samplerows, ","))
if (input$samplerows2=="")
{
dt1=sampleGeo()[,c('geo_accession',input$samplerows)]
}else{dt1=sampleGeo()[,c('geo_accession',input$samplerows, input$samplerows2)]}
row.names(dt1) <- dt1$geo_accession
dt1$geo_accession=NULL
dt1
})

 
  ### show the available groups for contrast matrix
 
# i <- reactiveVal()

# observeEvent(input$upload, { 
    # i(paste("The available group levels are: " ,
	# paste(as.character(levels(datareactive()$samples$group)), collapse=", "), sep=""))
  # })

# observeEvent(input$groupassignment, { 
    # i(updateTextInput(session, 'Group1', value=input$factor1) ,
	# updateTextInput(session, 'Group2', value=input$factor2))
  # }) 
# observeEvent(input$groupassignmentgeo, { 
    # i(updateTextInput(session, 'Group1', value=input$factor1geo) ,
	# updateTextInput(session, 'Group2', value=input$factor2geo))
  # }) 
 
# observeEvent(input$uploadSRP, { 
    # i(updateTextInput(session, 'Group1', value=input$factor1_SRP) ,
	# updateTextInput(session, 'Group2', value=input$factor2_SRP))
  # })

# output$grouplevels <- renderUI({
    # i()
  # })  
  output$grouplevels <- renderUI(
    
	if (input$upload){
	paste("The available group levels are: " ,
	paste(as.character(levels(datareactive()$samples$group)), collapse=", "), sep="")
	}else if (input$uploadGEO){
	if (input$groupassignment){
	updateTextInput(session, 'Group1', value=input$factor1)
	updateTextInput(session, 'Group2', value=input$factor2)
	}else if (input$groupassignmentgeo){
	updateTextInput(session, 'Group1', value=input$factor1geo)
	updateTextInput(session, 'Group2', value=input$factor2geo)
	}
	}else if (input$uploadSRP){
	updateTextInput(session, 'Group1', value=input$factor1_SRP)
	updateTextInput(session, 'Group2', value=input$factor2_SRP)
	}
	
  )
  

 
  ####Making contrast matrix
  
 contrastmatrix <- reactive({      
    if (input$degAnalysis)
	withProgress(message = 'Submit:', value = 0,{
	incProgress(.1, detail = paste("DEG analysis"))
	if (input$upload){
    Group=datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)  
    colnames(design) <- levels(Group)
    con2vscon1 = paste(as.character((input$Group2)), as.character((input$Group1)), sep="-")                
                              # levels = colnames(design))
	contr <- makeContrasts(contrasts=con2vscon1, levels = colnames(design))
	}else if (input$uploadGEO) {
	
	design <- model.matrix(~group + 0, geovisualize()$geovalues2)
	colnames(design) <- levels(geovisualize()$gs)
	   
	con2vscon1 = paste(geovisualize()$groups[2], geovisualize()$groups[1], sep="-")
	contr <- makeContrasts(contrasts=con2vscon1, levels = (design))
	}else if (input$uploadSRP){
	group= factor(sampleTableSRP()$group)
	design <- model.matrix(~0+group, data = sampleTableSRP()$sampleTable1)
	colnames(design) <- levels(sampleTableSRP()$group)
	con2vscon1 = paste(sampleTableSRP()$groups[2], sampleTableSRP()$groups[1], sep="-") 
	contr <- makeContrasts(contrasts=con2vscon1, levels = (design))
	}
	contr 
	})
	})
	
	
####Run Voom
	
 voomrun <- reactive ({
 if (input$upload){
 Group=datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)  
    #rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group) 
  	v <- voom(FilterReactive(),design=design, plot=F)
	}else if (input$uploadSRP){
	group= factor(sampleTableSRP()$group)
	print(group)
	design <- model.matrix(~0+group, data = sampleTableSRP()$sampleTable1)
	colnames(design) <- levels(sampleTableSRP()$group)
	v <- voom(filteredExprSRP()$SRPvalues2,design=design, plot=F)
	}
	v 
})

  
#### fit the linear model
 
 fitres <-reactive ({
 if (input$upload){
 Group=datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)  
    #rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group) 
	v1 <- voomrun()
    fit <- lmFit(v1, design)
	}else if (input$uploadSRP){
	group= factor(sampleTableSRP()$group)
	print(group)
	design <- model.matrix(~0+group, data = sampleTableSRP()$sampleTable1)
	colnames(design) <- levels(sampleTableSRP()$group)
	fit <- lmFit(voomrun(), design)
	}
	fit
	})

 ### apply your contrasts and ebayes method
 
voomreactive <- eventReactive(input$degAnalysis,{
if (input$uploadGEO)
{

design <- model.matrix(~group + 0, geovisualize()$geovalues2)
colnames(design) <- levels(geovisualize()$gs)
fit <- lmFit(geovisualize()$geovalues2, design) 
cfit <- contrasts.fit(fit, contrasts=contrastmatrix())
efit <- eBayes(cfit)
}else if (input$upload){
	Group=datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)  
    colnames(design) <- levels(Group) 
	cfit <- contrasts.fit(fitres(), contrasts=contrastmatrix())
    efit <- eBayes(cfit)
	}else if (input$uploadSRP){
	group= factor(sampleTableSRP()$group)
	print(group)
	design <- model.matrix(~0+group, data = sampleTableSRP()$sampleTable1)
	colnames(design) <- levels(sampleTableSRP()$group)
	print (design)
	cfit <- contrasts.fit(fitres(), contrasts=contrastmatrix())
    efit <- eBayes(cfit)
	}
	(efit)
  }) 
  
  
 ###and plot 
 output$voommeanvariance <- renderPlot({
    if (input$degAnalysis){
	      plotSA(voomreactive(), main="Final model: Mean-variance trend")
      }})
 

##### TopTable of DEGs Ebayes method
  
 output$reslimma <- renderDT({
        res_con2vcon1_limma <- topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf)
		if (input$uploadGEO){
		res <- subset(res_con2vcon1_limma, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol",'Gene.ID',"Gene.title"))
		res=res[!(res$Gene.symbol==""),]
		}else if (input$uploadSRP){
		res <- subset(res_con2vcon1_limma, select=c("gene_id","symbol","adj.P.Val","P.Value","t","B","logFC"))
				}else{res = res_con2vcon1_limma[,c('ENSEMBL', 'ENTREZID',  'SYMBOL', 'logFC','P.Value', 'adj.P.Val','t')]}
		res$logFC <- round(res$logFC, digits = 4)
		#res$P.Value <- round(res$P.Value)
		res$adj.P.Val <- round(res$adj.P.Val, digits = 5)
		res
	 
  }, filter = 'top',
    rownames = FALSE,
# options = list(
  # autoWidth = TRUE,
  # columnDefs = list(list(width = '100px')))	
  )
  
 #### summary of ebayes degs
 
 output$summarydegs <- renderTable({
 if (input$pvalfilter){
  dt1 <- as.data.frame(summary(decideTests(voomreactive(),p.value=input$pvalcutoff)))
  }else{dt1 <- as.data.frame(summary(decideTests(voomreactive())))}
  colnames(dt1) <- c('Regulation',"group","Frequency")
  dt1
  })
  
  #### voom download
  output$voomDownload <- downloadHandler(
    filename = function() {
      paste("voom-result", ".csv", sep = "")
    },
    content = function(file) {
	
	res_con2vcon1_limma <- topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf)
	if (input$uploadGEO){
		res <- res <- subset(res_con2vcon1_limma, select=c("ID","adj.P.Val","P.Value","t","logFC","Gene.symbol",'Gene.ID',"Gene.title"))
		}else if (input$uploadSRP){
		res <- subset(res_con2vcon1_limma, select=c("gene_id","symbol","adj.P.Val","P.Value","t","B","logFC"))
				}else{
		res = res_con2vcon1_limma[,c('ENSEMBL', 'ENTREZID',  'SYMBOL', 'logFC','P.Value', 'adj.P.Val','t')]}
      write.csv(res, file, row.names = FALSE)
	  #write.csv(voomreactive(), file)
    }
  )
  
  
  ##### Treat your degs 
  
  #### reactive treat
  treatvals <- eventReactive(input$treatanalysis,
  
  { 
  withProgress(message = 'Submit:', value = 0,{
	incProgress(.1, detail = paste("treating dataset"))
  tfit <- treat(voomreactive(),lfc=log2(as.numeric(input$FC)))
  })
  })
  
    
  ### datatable of treat method
   output$restreat <- renderDT({
   		tfit2 <-treatvals()
        res_con2vcon1_treat <- topTreat(tfit2,sort.by="t", coef=1, n=Inf)
		if (input$uploadGEO){
		res1 <- subset(res_con2vcon1_treat, select=c("ID","adj.P.Val","P.Value","t","logFC","Gene.symbol",'Gene.ID',"Gene.title"))
		res1=res1[!(res1$Gene.symbol==""),]
		}else if (input$uploadSRP){
		res1 <- subset(res_con2vcon1_treat, select=c("gene_id","symbol","adj.P.Val","P.Value","t","logFC"))
				}else{res1 = res_con2vcon1_treat[,c('ENSEMBL', 'ENTREZID',  'SYMBOL', 'logFC','P.Value', 'adj.P.Val','t')]}
	    res1$logFC <- round(res1$logFC, digits = 4)
		#res1$P.Value <- round(res1$P.Value)
		res1$adj.P.Val <- round(res1$adj.P.Val, digits = 5)
		res1 
  }, filter = 'top',
    rownames = FALSE 
  )
  
  #### summary of treat analysis
  output$summaryaftertreat <- renderTable({
   tfit2 <-treatvals()
  dt2 <- as.data.frame(summary(decideTests(tfit2,p.value=input$pvalcutoff2)))
  colnames(dt2) <- c('Regulation',"group","Number")
  dt2
  })

#### treat download
output$treatDownload <- downloadHandler(
    filename = function() {
      paste("treat-result", ".csv", sep = "")
    },
    content = function(file) {
	tfit2 <-treatvals()
        res_con2vcon1_treat <- topTreat(tfit2,sort.by="t", coef=1, n=Inf)
		if (input$uploadGEO){
		res1 <- subset(res_con2vcon1_treat, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol",'Gene.ID',"Gene.title"))
		}else if (input$uploadSRP){
		res1 <- subset(res_con2vcon1_treat, select=c("gene_id","symbol","adj.P.Val","P.Value","t","logFC"))
				}else{res1 = res_con2vcon1_treat[,c('ENSEMBL', 'ENTREZID',  'SYMBOL', 'logFC','P.Value', 'adj.P.Val','t')]}
		#res1 = res_con2vcon1_treat[,c('ENSEMBL', 'ENTREZID',  'SYMBOL', 'logFC','P.Value', 'adj.P.Val','t')]
    write.csv(res1, file, row.names = FALSE)  
    }
  )


#### MD plot after treat analysis
   output$MDplot <- renderPlot({
   tfit2 <-treatvals()
   dt2 <- (decideTests(tfit2))
   plotMD(object=tfit2, column=1, status=dt2[,1], main=paste(as.character((input$Group2)), as.character((input$Group1)), sep="V"), 
       xlim=c(-8,13))
	   })
  
 output$MDplotDownload <-downloadHandler(
	filename = function(){
		paste("MDplot", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 tfit2 <-treatvals()
   dt2 <- (decideTests(tfit2))
   plotMD(object=tfit2, column=1, status=dt2[,1], main=paste(as.character((input$Group2)), as.character((input$Group1)), sep="V")) 
 
  dev.off()
 },
 contentType = 'image/png'
 )
 
  ##### Run Topconfects
  
  ###reactive topconfects
  confectvals <- eventReactive(input$runconfect,
  { 
  withProgress(message = 'Submit:', value = 0,{
	incProgress(.1, detail = paste("running topconfects"))
  confects_1 <- limma_confects(voomreactive(), fdr=as.numeric(input$fdr))
  confects_1 
  })
  })
  
  ### datatable of topconfects
  output$restopconfect <- renderDT({
  		confect2 <-confectvals()
		if (input$uploadGEO){
		#conf_res=confect2$table
		conf_res=subset(confect2$table, select=c('rank','ID',"Gene symbol",'Gene ID',"Gene title", 'confect','effect','AveExpr'))
		#, select=c('rank','ID',"Gene.symbol","Gene.title", 'confect','effect','AveExpr'))
		#[,c('ID',"Gene.symbol","Gene.title", 'confect','effect','AveExpr')]
		}else if (input$uploadSRP){
		conf_res=subset(confect2$table, select=c('rank','gene_id',"symbol", 'confect','effect','AveExpr'))
		}else{
		conf_res=confect2$table[, c('rank','ENSEMBL', 'ENTREZID', 'SYMBOL', 'confect','effect','AveExpr')]
		}
		conf_res
		}, filter = 'top',
    rownames = FALSE 
  )
  
  output$topconfectsDownload <- downloadHandler(
    filename = function() {
      paste("topconfect-result", ".csv", sep = "")
    },
    content = function(file) {
	confect2 <-confectvals()
	if (input$uploadGEO){
		conf_res=subset(confect2$table, select=c('rank','ID',"Gene symbol",'Gene ID',"Gene title", 'confect','effect','AveExpr'))
				}else if (input$uploadSRP){
		conf_res=subset(confect2$table, select=c('rank','gene_id',"symbol", 'confect','effect','AveExpr'))
		}else{
		conf_res=confect2$table[, c('rank','ENSEMBL', 'ENTREZID', 'SYMBOL', 'confect','effect','AveExpr')]
		}    
    write.csv(conf_res, file, row.names = FALSE)  
    }
  )
  
 #### summary of topconfects  
  
 ##### MD plot after Topconfects
 
 observeEvent(input$plotlimma_confect,{
  withProgress(message = 'MD plot:', value = 0,{
	 incProgress(.1, detail = paste("ploting"))
	
output$MD_limma_confects <- renderPlot({
isolate({
confect2 <-confectvals()
top_1 <- topTable(voomreactive(), n=Inf)
if (input$uploadGEO){
design <- model.matrix(~group + 0, geovisualize()$geovalues2)
colnames(design) <- levels(geovisualize()$gs)
fit <- lmFit(geovisualize()$geovalues2, design)  
plotMD(fit , legend="bottomleft", status=paste0(
   ifelse(rownames(voomreactive()) %in% rownames(top_1)[1:input$n], "topTable ",""),
   ifelse(rownames(voomreactive()) %in% confect2$table$name[1:input$n], "confects ","")))
   }else{
      plotMD(fitres() , legend="bottomleft", status=paste0(
   ifelse(rownames(voomreactive()) %in% rownames(top_1)[1:input$n], "topTable ",""),
   ifelse(rownames(voomreactive()) %in% confect2$table$name[1:input$n], "confects ","")))}
       })
	   })
	   })
	   }) 


output$MD_limma_confectsDownload <-downloadHandler(
	filename = function(){
		paste("MDplot_topconfects", ".png", sep="")
 },
 content = function(file){
 png(file)
 confect2 <-confectvals()
   top_1 <- topTable(voomreactive(), n=Inf)
   if (input$uploadGEO){
design <- model.matrix(~group + 0, geovisualize()$geovalues2)
colnames(design) <- levels(geovisualize()$gs)
fit <- lmFit(geovisualize()$geovalues2, design)  
plotMD(fit , legend="bottomleft", status=paste0(
   ifelse(rownames(voomreactive()) %in% rownames(top_1)[1:input$n], "topTable ",""),
   ifelse(rownames(voomreactive()) %in% confect2$table$name[1:input$n], "confects ","")))
   }else{
   plotMD(fitres() , legend="bottomleft", status=paste0(
   ifelse(rownames(voomreactive()) %in% rownames(top_1)[1:input$n], "topTable ",""),
   ifelse(rownames(voomreactive()) %in% confect2$table$name[1:input$n], "confects ","")))}
   dev.off()
 },
 contentType = 'image/png'
 )
 

#### Plot the topconfects
output$confectplot <- renderPlot({
 confects_plot(confectvals())
 })
 
output$confectplotDownload <-downloadHandler(
	filename = function(){
		paste("confectplot", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=confects_plot(confectvals())
 print (p)
  dev.off()
 },
 contentType = 'image/png'
 )
 
 
##### compare voom and topconfect
 output$voomvsconfectplot <- renderPlot({
top_1 <- topTable(voomreactive(), n=Inf)
rank_rank_plot(confectvals()$table$name, rownames(top_1), "limma_confects", "topTable")
})

output$voomvsconfectplotDownload <-downloadHandler(
	filename = function(){
		paste("voomvsconfectplot", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 top_1 <- topTable(voomreactive(), n=Inf)
p=rank_rank_plot(confectvals()$table$name, rownames(top_1), "limma_confects", "topTable")
print (p)
  dev.off()
 },
 contentType = 'image/png'
 )

######### Enrichment analysis##########

## Camera Competitive Enrichment Analysis ##
#Run Camera geneset enrichment on a comparison
# Select ENTREZID
 
##### pathway reactive entrezid 
 pathwayIDreactive <- reactive({
    switch(input$pathwaysname,
           "All Gene_Ontology" = Hs.GO.Symbol,
		   "Gene Ontology: Biological Process (Full)" = Hs.GOBP.Symbol,
           "Gene Ontology: Cellular Component (Full)" = Hs.GOCC.full.Symbol,
		   "Gene Ontology: Molecular Function (Full)" = Hs.GOMF.Symbol,
		   "Human Phenotype Ontology" = Hs.HPO.Symbol,
		   "Reactome"=Hs.Reactome.Symbol,
		   "MSigDB Hallmark"=Hs.Hallmark.Symbol,
           "BioCarta" = Hs.Biocarta.Symbol,
		   "KEGG" = Hs.KEGG.Symbol,
			"PID" = Hs.PID.Symbol,
			"WikiPathways" = Hs.WikiPathways.Symbol,
			"MSigDB Chemical and Genetic Perturbations" = Hs.CGP.Symbol,
			"MSigDB Computational Genesets" = Hs.Comp.Symbol,
			"MSigDB Oncogenic Signature Genesets" = Hs.Oncogenic.Symbol,
			"MSigDB Immunologic signature Genesets" = Hs.Immune.Symbol,
			"MSigDB Cell Types" =Hs.CellType.Symbol)
  })
  
 
msigdbgeneset <- reactive({
if (input$uploadGEO)
  {
  
design <- model.matrix(~group + 0, geovisualize()$geovalues2)
colnames(design) <- levels(geovisualize()$gs)
fit <- lmFit(geovisualize()$geovalues2, design)

  #voomrungeo = voom(geovisualize()$ex)
  #print (input$pathwaysname)
  idx <- ids2indices(pathwayIDreactive(),id=fit[!fit[["genes"]][["Gene title"]]=="",]$genes$`Gene symbol`)
    }else if (input$upload){
	
	idx <- ids2indices(pathwayIDreactive(),id=voomrun()$genes$SYMBOL)
	}else if (input$uploadSRP){
	
	idx <- ids2indices(pathwayIDreactive(),id=voomrun()$genes$symbol)
	}
return(idx)
})

 
#### datatable of CAMERA method
camera_out <- reactive({
if (input$pathwayenrichment)

withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("Enrichment analysis"))
	if (input$uploadGEO){
 design <- model.matrix(~group + 0, geovisualize()$geovalues2)
colnames(design) <- levels(geovisualize()$gs)
cam.HsGO.1 <- camera(geovisualize()$ex,msigdbgeneset(),design,contrast=contrastmatrix())
  }else if (input$uploadSRP){
  design <- model.matrix(~0+group, data = sampleTableSRP()$sampleTable1)
	colnames(design) <- levels(sampleTableSRP()$group)
	cam.HsGO.1 <- camera(voomrun(),msigdbgeneset(),design,contrast=contrastmatrix())}
  else{
	Group=datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f) 
#Run the camera method from limma, make sure to define the contrast you want to analyze from you contrast matrix - in this case, "1" is Con2vCon1
cam.HsGO.1 <- camera(voomrun(),msigdbgeneset(),design,contrast=contrastmatrix())}
# cam.HsGO.1$pval <- round(cam.HsGO.1$pval, digits = 3)
# cam.HsGO.1$padj <- round(cam.HsGO.1$padj, digits = 3)
# cam.HsGO.1$NES <- round(cam.HsGO.1$NES, digits = 3)
cam.HsGO.1
})
})

output$enrichmentout <- renderDT ({
camera_out()
}, filter = 'top'
)
 
 
output$enrichmentontologyDownload <- downloadHandler(
    filename = function() {
      paste("enrichmentontology-result", ".csv", sep = "")
    },
	content = function(file) {
    write.csv(camera_out(), file, row.names = TRUE)  
    }
  )

barcodeplotreactive <- eventReactive(input$barcode,{
id2 = msigdbgeneset()
print (head(id2))
id=id2[[input$geneset]]
print (input$geneset)
print (id)
id
})
output$barcodeout <- renderPlot({

barcodeplot(voomreactive()$t[,1], 
           index=barcodeplotreactive(), 
			# # #index2 = msigdbgeneset()$input$geneset2,
			main = input$geneset)
           #main=paste(as.character((input$Group2)), as.character((input$Group1)), sep=" Vs "))
     
})

output$barcodeDownload <- downloadHandler(
    filename = function(){
		paste("barcodeplot", ".png", sep="")
 },
 content = function(file){
 png(file)
 barcodeplot(voomreactive()$t[,1], 
           index=barcodeplotreactive(), 
			# # #index2 = msigdbgeneset()$input$geneset2,
			main = input$geneset)
   dev.off()
 },
 contentType = 'image/png'
 )

  
####Enrichment analysis using fgsea  

#### reactive element of fgsea

## change the genelist to include either ebayes, treat or topconfects
genelistreactive <- reactive({
    switch(input$genelist,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
           "topConfects" = confectvals()$table)
  })
	 
 
fgseaoutreactive <- eventReactive(input$runfgsea,
    { 
    res_con2vcon1 <- genelistreactive()
	if (input$genelist == "eBayes_tvalue" | input$genelist == "TREAT_tvalue"){
		#res = res_con2vcon1[,c('ENSEMBL', 'ENTREZID', 't', 'SYMBOL', 'logFC','P.Value', 'adj.P.Val')]
		res_con2vcon1$logFC <- round(res_con2vcon1$logFC, digits = 4)
		#res_con2vcon1$P.Value <- round(res_con2vcon1$P.Value)
		res_con2vcon1$adj.P.Val <- round(res_con2vcon1$adj.P.Val, digits = 5)
	if (input$uploadGEO)
  {
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
	ranks_limma1 <- ranks_limma[, c('Gene.symbol','t')]
  }else if (input$uploadSRP){
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
	ranks_limma1 <- ranks_limma[, c('symbol','t')]
  }else{
		res_con2vcon1 <- subset(res_con2vcon1, SYMBOL != "" )
	res_con2vcon1 <- subset(res_con2vcon1, !duplicated(SYMBOL))
ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
ranks_limma1 <- ranks_limma[, c('SYMBOL','t')]}}

if (input$genelist == "topConfects")
{
if (input$uploadGEO)
  {
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
	ranks_limma1 <- ranks_limma[, c('Gene symbol','confect')]
  }else if (input$uploadSRP){
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
	ranks_limma1 <- ranks_limma[, c('symbol','confect')]
  }else{
res_con2vcon1 <- subset(res_con2vcon1, SYMBOL != "" )
res_con2vcon1 <- subset(res_con2vcon1, !duplicated(SYMBOL))
# coerce confects with a NA score to 0
res_con2vcon1[is.na(res_con2vcon1)] <- 0
ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
ranks_limma1 <- ranks_limma[, c('SYMBOL','confect')]}
}
return(ranks_limma1)
 })
 

  
 #### change pathways
 pathwayreactive <- reactive({
    switch(input$pathwaylist,
           "All Gene_Ontology" = Hs.GO.Symbol,
		   "Gene Ontology: Biological Process (Full)" = Hs.GOBP.Symbol,
           "Gene Ontology: Cellular Component (Full)" = Hs.GOCC.full.Symbol,
		   "Gene Ontology: Molecular Function (Full)" = Hs.GOMF.Symbol,
		   "Human Phenotype Ontology" = Hs.HPO.Symbol,
		   "Reactome"=Hs.Reactome.Symbol,
		   "MSigDB Hallmark"=Hs.Hallmark.Symbol,
           "BioCarta" = Hs.Biocarta.Symbol,
		   "KEGG" = Hs.KEGG.Symbol,
			"PID" = Hs.PID.Symbol,
			"WikiPathways" = Hs.WikiPathways.Symbol,
			"MSigDB Chemical and Genetic Perturbations" = Hs.CPG.Symbol,
			"MSigDB Computational Genesets" = Hs.Comp.Symbol,
			"MSigDB Oncogenic Signature Genesets" = Hs.Oncogenic.Symbol,
			"MSigDB Immunologic signature Genesets" = Hs.Immune.Symbol,
			"MSigDB Cell Types" =Hs.CellType.Symbol)
  })
  
### output of fgsea
fgsearesult <- reactive ({
  withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running fgsea"))
 ranks <- deframe(fgseaoutreactive())
 # Run fGSEA 
  fgseaRes <- fgseaMultilevel(pathways=pathwayreactive(), stats=ranks, minSize = 15, maxSize = 500, eps=0)
# Make the Results tidy
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(desc(NES))
fgseaResTidy
})
  })
  
output$fgseaout_gobp <- renderDT ({
 fgsearesult()
 },filter = 'top',
    rownames = FALSE)

output$fgseaDownload <- downloadHandler(
    filename = function() {
      paste("fgsea-result", ".csv", sep = "")
    },
	content = function(file) {
	res=fgsearesult()
    write.csv(res, file, row.names = FALSE)  
    }
  )
 
fgseaplotreactive <- reactive ({
fgRes=fgsearesult()
fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                      "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="fGSEA Results") + 
    theme_minimal()
}) 

output$fgseaplot <- renderPlot({
fgseaplotreactive()
})

output$fgseaplotDownload <- downloadHandler(
    filename = "fgseaplot.png",
		
 content = function(file){
 p=fgseaplotreactive()
 png(file,height=780,width=780)
 
 print(p)
   dev.off()
 }
 )
 
 
 
  
#########fgsea TF#############

genelistTFreactive <- reactive({
    switch(input$genelistTF,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
           "topConfects" = confectvals()$table)
  })
	 
 
fgseaoutreactiveTF <- eventReactive(input$runfgseaTF,
    { 
    res_con2vcon1 <- genelistTFreactive()
	if (input$genelistTF == "eBayes_tvalue" | input$genelistTF == "TREAT_tvalue"){
		#res = res_con2vcon1[,c('ENSEMBL', 'ENTREZID', 't', 'SYMBOL', 'logFC','P.Value', 'adj.P.Val')]
		res_con2vcon1$logFC <- round(res_con2vcon1$logFC, digits = 4)
		#res_con2vcon1$P.Value <- round(res_con2vcon1$P.Value)
		res_con2vcon1$adj.P.Val <- round(res_con2vcon1$adj.P.Val, digits = 5)
	if (input$uploadGEO)
  {
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
	ranks_limma1 <- ranks_limma[, c('Gene.symbol','t')]
  }else if (input$uploadSRP){
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
	ranks_limma1 <- ranks_limma[, c('symbol','t')]
  }else{
		res_con2vcon1 <- subset(res_con2vcon1, SYMBOL != "" )
	res_con2vcon1 <- subset(res_con2vcon1, !duplicated(SYMBOL))
ranks_limma <- res_con2vcon1[order(res_con2vcon1$t), ]
ranks_limma1 <- ranks_limma[, c('SYMBOL','t')]}}

if (input$genelistTF == "topConfects")
{
if (input$uploadGEO)
  {
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
	ranks_limma1 <- ranks_limma[, c('Gene symbol','confect')]
  }else if (input$uploadSRP){
  res_con2vcon1=res_con2vcon1
  ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
	ranks_limma1 <- ranks_limma[, c('symbol','confect')]
  }else{
res_con2vcon1 <- subset(res_con2vcon1, SYMBOL != "" )
res_con2vcon1 <- subset(res_con2vcon1, !duplicated(SYMBOL))
# coerce confects with a NA score to 0
res_con2vcon1[is.na(res_con2vcon1)] <- 0
ranks_limma <- res_con2vcon1[order(res_con2vcon1$confect), ]
ranks_limma1 <- ranks_limma[, c('SYMBOL','confect')]}
}
return(ranks_limma1)

  })
 

 #### change pathways
 pathwayreactiveTF <- reactive({
    switch(input$pathwaylistTF,
           "ENCODE-ChEA Consensus (Enrichr)" = Hs.ECC,
           "ChEA 2016 (Enrichr)"= Hs.ChEA2016,
           "ENCODE 2015 (Enrichr)" = Hs.ENCODE,
           "ReMap ChIP-seq 2018 Human" = Hs.ReMap,
           "TRRUST 2019 Human" = Hs.TRRUST,
           "ChEA3 Literature ChIP-Seq"= Hs.Literature,
           "TRANSFAC/JASPAR PWMs (Enrichr)" = Hs.TRANSFACJASPAR,
           "Gene Transcription Regulation Database (GTRD v20.06)" = Hs.GTRD.sel,
           "MSigDB Legacy TF targets" = Hs.TFLegacy.sel,
           "miRTarBase 2017" = Hs.miR,
           "miRNA TargetScan 2017" = Hs.miRTargetScan,
           "miRDB v6.0" = Hs.miRNA.sel)
  })
  
### output of fgsea
fgsearesultTF <- reactive ({
  withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running fgsea"))
 ranks <- deframe(fgseaoutreactiveTF())
 # Run fGSEA 
  fgseaRes <- fgseaMultilevel(pathways=pathwayreactiveTF(), stats=ranks, minSize = 15, maxSize = 500)
# Make the Results tidy
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(desc(NES))
fgseaResTidy
})
  })
  
output$fgseaout_TF <- renderDT ({
 fgsearesultTF()
 },filter = 'top',
    rownames = FALSE)

output$fgseaTFDownload <- downloadHandler(
    filename = function() {
      paste("fgseaTF-result", ".csv", sep = "")
    },
	content = function(file) {
	res=fgsearesultTF()
    write.csv(res, file, row.names = FALSE)  
    }
  )

##### Enrichr TF########

genesreactive <- reactive({
    switch(input$genes,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
		   "topConfects" = confectvals()$table)
  })

enrichrreactive <- eventReactive(input$runenrichr,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running enrichr"))
res_limma <- genesreactive()
if (input$genes == "eBayes_tvalue" | input$genes == "TREAT_tvalue"){
res_limma_down <- res_limma[which(res_limma$logFC < -as.numeric(input$FC3)), ]
res_limma_up <- res_limma[which(res_limma$logFC > input$FC3), ]	
}else if (input$genes == "topConfects")
{
res_limma_down <- res_limma[which(res_limma$confect < -as.numeric(input$FC3)), ]
res_limma_up <- res_limma[which(res_limma$confect > input$FC3), ]	
}
if (input$uploadGEO){
if ("Gene.symbol" %in% colnames(res_limma)){
down_genes <- as.character(res_limma_down$Gene.symbol)
up_genes <- as.character(res_limma_up$Gene.symbol)
}else if ("Gene symbol" %in% colnames(res_limma)){
down_genes <- as.character(res_limma_down$`Gene symbol`)
up_genes <- as.character(res_limma_up$`Gene symbol`)
}

}else if (input$uploadSRP){
down_genes <- as.character(res_limma_down$symbol)
up_genes <- as.character(res_limma_up$symbol)
}else{
down_genes <- as.character(res_limma_down$SYMBOL)
up_genes <- as.character(res_limma_up$SYMBOL)
}

down_enriched <- enrichr(c(down_genes), input$database)
up_enriched <- enrichr(c(up_genes), input$database)

list(down_enriched1 = down_enriched, up_enriched1 =up_enriched)
})})


down_enrich_res <- reactive ({
down_enriched_ECC <- enrichrreactive()$down_enriched1[[input$database]] 
down_enriched_ECC <- down_enriched_ECC[order(down_enriched_ECC$Combined.Score), ]  
#if (!is.null(down_enriched_ECC$Combined.Score)){
down_enriched_ECC_Tidy <- down_enriched_ECC %>%
  as_tibble() %>%
  arrange(desc(Combined.Score))
down_enriched_ECC_Tidy=dplyr::select(down_enriched_ECC_Tidy, -Old.P.value, -Old.Adjusted.P.value)
down_enriched_ECC_Tidy
#}
})

output$down_enrichr <- renderDT({
#if (!length(down_enrich_res()==0)){
down_enrich_res()
  #}
  }, filter = 'top',
    rownames = FALSE)
  
output$downenrichrDownload <- downloadHandler(
    filename = function() {
      paste("downenrichr-result", ".csv", sep = "")
    },
	content = function(file) {
	res= down_enrich_res()
    write.csv(res, file, row.names = FALSE)  
    }
  )

up_enrich_res <- reactive ({
up_enriched_ECC <- enrichrreactive()$up_enriched1[[input$database]] 
up_enriched_ECC <- up_enriched_ECC[order(up_enriched_ECC$Combined.Score), ] 
#if (!is.null(up_enriched_ECC$Combined.Score)){ 
up_enriched_ECC_Tidy <- up_enriched_ECC %>%
  as_tibble() %>%
  arrange(desc(Combined.Score))
up_enriched_ECC_Tidy=dplyr::select(up_enriched_ECC_Tidy, -Old.P.value, -Old.Adjusted.P.value)
up_enriched_ECC_Tidy
#}
})

 
output$up_enrichr <- renderDT({
#if (!length(up_enrich_res()==0)){
up_enrich_res()
  #}
  }, filter = 'top',
    rownames = FALSE)

output$upenrichrDownload <- downloadHandler(
    filename = function() {
      paste("upenrichr-result", ".csv", sep = "")
    },
	content = function(file) {
	res= up_enrich_res()
    write.csv(res, file, row.names = FALSE)  
    })


#####Enrichr TF plots########

upTFplotreactive <- eventReactive(input$filterenrichrTF,{
toplot = up_enrich_res()
ggplot(data= toplot[1:input$num2,], aes(x=.data[['Term']], y=-log10(.data[[input$pval_FDR]]),fill=.data[[input$pval_FDR]]))+
  geom_bar(stat="identity",color="black", width=0.5)+ coord_flip()
})

output$upTFplots <- renderPlot({
upTFplotreactive()
})

output$upTFplotsDownload <- downloadHandler(
    filename = function(){
		paste("upTFplots", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=upTFplotreactive()
 print (p)
   dev.off()
 },
 contentType = 'image/png'
 )
 

downTFplotreactive <- eventReactive(input$filterenrichrTF,{
toplot = down_enrich_res()
ggplot(data= toplot[1:input$num2,], aes(x=.data[['Term']], y=-log10(.data[[input$pval_FDR]]),fill=.data[[input$pval_FDR]]))+
  geom_bar(stat="identity",color="black", width=0.5)+ coord_flip()
})

output$downTFplots <- renderPlot({
downTFplotreactive()
})

output$downTFplotsDownload <- downloadHandler(
    filename = function(){
		paste("downTFplots", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=downTFplotreactive()
 print (p)
   dev.off()
 },
 contentType = 'image/png'
 )
 

updownTFplotsreactive <- eventReactive(input$filterenrichrTF,{
up_enriched_ECC <- enrichrreactive()$up_enriched1[[input$database]] 
down_enriched_ECC <- enrichrreactive()$down_enriched1[[input$database]]  
up_enriched_ECC$type <- "up"
down_enriched_ECC$type <- "down"
up_enriched_ECC <- up_enriched_ECC[order(up_enriched_ECC$Combined.Score), ]  # sort
down_enriched_ECC <- down_enriched_ECC[order(down_enriched_ECC$Combined.Score), ]  
down_enriched_ECC$Combined.Score <- (-1) * down_enriched_ECC$Combined.Score
down_enriched_ECC1=down_enriched_ECC[1:input$num2,]
up_enriched_ECC1=up_enriched_ECC[1:input$num2,]
gos <- rbind(down_enriched_ECC1,up_enriched_ECC1)
ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) +
geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
scale_fill_manual(name="Expression",
labels = c("Down regulated", "Up regulated"),
values = c("down"="#00ba38", "up"="#f8766d")) +
labs(title=paste0("Combined scores from ", input$databaseOntology)) +
coord_flip()
})

output$updownTFplots <- renderPlot({
updownTFplotsreactive()
})
	
output$updownTFplotsDownload <- downloadHandler(
    filename = function(){
		paste("updownTFplots", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=updownTFplotsreactive()
 print(p)
   dev.off()
 },
 contentType = 'image/png'
 )
###############	
	
##############################
##### Enrichr Ontology########

genesreactiveenrichr <- reactive({
    switch(input$genesenrichr,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
		   "topConfects" = confectvals()$table)
  })

enrichrreactiveOntology <- eventReactive(input$runenrichrOntology,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running enrichr"))
res_limma <- genesreactiveenrichr()
if (input$genesenrichr == "eBayes_tvalue" | input$genesenrichr == "TREAT_tvalue" ){
res_limma_down <- res_limma[which(res_limma$logFC < -as.numeric(input$FC4)), ]
res_limma_up <- res_limma[which(res_limma$logFC > input$FC4), ]	
}else if (input$genesenrichr == "topConfects")
{
res_limma_down <- res_limma[which(res_limma$confect < -as.numeric(input$FC4)), ]
res_limma_up <- res_limma[which(res_limma$confect > input$FC4), ]	
}
if (input$uploadGEO){
if ("Gene.symbol" %in% colnames(res_limma)){
down_genes <- as.character(res_limma_down$Gene.symbol)
up_genes <- as.character(res_limma_up$Gene.symbol)
}else if ("Gene symbol" %in% colnames(res_limma)){
down_genes <- as.character(res_limma_down$`Gene symbol`)
up_genes <- as.character(res_limma_up$`Gene symbol`)
}
}else if (input$uploadSRP){
down_genes <- as.character(res_limma_down$symbol)
up_genes <- as.character(res_limma_up$symbol)
}else{
down_genes <- as.character(res_limma_down$SYMBOL)
up_genes <- as.character(res_limma_up$SYMBOL)
}

down_enriched <- enrichr(c(down_genes), input$databaseOntology)
up_enriched <- enrichr(c(up_genes), input$databaseOntology)

list(down_enriched = down_enriched, up_enriched =up_enriched)
})})


down_enrich_resOntology <- reactive ({
down_enriched_ECC <- enrichrreactiveOntology()$down_enriched[[input$databaseOntology]]  
#if (!is.null(down_enriched_ECC$Combined.Score)){

down_enriched_ECC <- down_enriched_ECC[order(down_enriched_ECC$Combined.Score), ] 
down_enriched_ECC_Tidy <- down_enriched_ECC %>%
  as_tibble() %>%
  arrange(desc(Combined.Score))
down_enriched_ECC_Tidy=dplyr::select(down_enriched_ECC_Tidy, -Old.P.value, -Old.Adjusted.P.value)

down_enriched_ECC_Tidy
#}
})

output$down_ontologyenrichr <- renderDT({
#if (!length(down_enrich_resOntology()==0)){
down_enrich_resOntology()
  #}
  }, filter = 'top',
    rownames = FALSE)
  
output$downenrichrontologyDownload <- downloadHandler(
    filename = function() {
      paste("downenrichrOntology-result", ".csv", sep = "")
    },
	content = function(file) {
	res= down_enrich_resOntology()
    write.csv(res, file, row.names = FALSE)  
    }
  )

up_enrich_resOntology <- reactive ({
up_enriched_ECC <- enrichrreactiveOntology()$up_enriched[[input$databaseOntology]] 
#if (!is.null(up_enriched_ECC$Combined.Score)){ 
up_enriched_ECC <- up_enriched_ECC[order(up_enriched_ECC$Combined.Score), ] 
up_enriched_ECC_Tidy <- up_enriched_ECC %>%
  as_tibble() %>%
  arrange(desc(Combined.Score))
up_enriched_ECC_Tidy=dplyr::select(up_enriched_ECC_Tidy, -Old.P.value, -Old.Adjusted.P.value)

up_enriched_ECC_Tidy
#}
})

 
output$up_ontologyenrichr <- renderDT({
#if (!length(up_enrich_resOntology()==0)){
up_enrich_resOntology()
 # }
  }, filter = 'top',
    rownames = FALSE)

output$upenrichrontologyDownload <- downloadHandler(
    filename = function() {
      paste("upenrichrOntology-result", ".csv", sep = "")
    },
	content = function(file) {
	res= up_enrich_resOntology()
    write.csv(res, file, row.names = FALSE)  
    })

#######Enrichr Ontology Plots########
	
upplotreactive <- eventReactive(input$filterenrichrOntology,{
toplot = up_enrich_resOntology()
ggplot(data= toplot[1:input$num,], aes(x=.data[['Term']], y=-log10(.data[[input$pval_FDR2]]),fill=.data[[input$pval_FDR2]]))+
  geom_bar(stat="identity",color="black", width=0.5)+ coord_flip()
})

output$upplots <- renderPlot({
upplotreactive()
})

output$upplotsenrichrOntologyDownload <- downloadHandler(
    filename = function(){
		paste("upplotsenrichrOntology", ".png", sep="")
 },
 content = function(file){
 png(file, height=780,width=780)
 p=upplotreactive()
 print(p)
   dev.off()
 },
 contentType = 'image/png'
 )


downplotreactive <- eventReactive(input$filterenrichrOntology,{
toplot = down_enrich_resOntology()
ggplot(data= toplot[1:input$num,], aes(x=.data[['Term']], y=-log10(.data[[input$pval_FDR2]]),fill=.data[[input$pval_FDR2]]))+
  geom_bar(stat="identity",color="black", width=0.5)+ coord_flip()
})

output$downplots <- renderPlot({
downplotreactive()
})

output$downplotsenrichrOntologyDownload <- downloadHandler(
    filename = function(){
		paste("downplotsenrichrOntology", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=downplotreactive()
 print (p)
   dev.off()
 },
 contentType = 'image/png'
 )


updownplotsreactive <- eventReactive(input$filterenrichrOntology,{
up_enriched_ECC <- enrichrreactiveOntology()$up_enriched[[input$databaseOntology]] 
down_enriched_ECC <- enrichrreactiveOntology()$down_enriched[[input$databaseOntology]]  
up_enriched_ECC$type <- "up"
down_enriched_ECC$type <- "down"
up_enriched_ECC <- up_enriched_ECC[order(up_enriched_ECC$Combined.Score), ]  # sort
down_enriched_ECC <- down_enriched_ECC[order(down_enriched_ECC$Combined.Score), ]  
down_enriched_ECC$Combined.Score <- (-1) * down_enriched_ECC$Combined.Score
down_enriched_ECC1=down_enriched_ECC[1:input$num,]
up_enriched_ECC1=up_enriched_ECC[1:input$num,]
gos <- rbind(down_enriched_ECC1,up_enriched_ECC1)
ggplot(gos, aes(x=Term, y=Combined.Score , label=Combined.Score)) +
geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
scale_fill_manual(name="Expression",
labels = c("Down regulated", "Up regulated"),
values = c("down"="#00ba38", "up"="#f8766d")) +
labs(title=paste0("Combined scores from ", input$databaseOntology)) +
coord_flip()
})

output$updownplots <- renderPlot({
updownplotsreactive()
})

output$updownplotsenrichrOntologyDownload <- downloadHandler(
    filename = function(){
		paste("updownplotsenrichrOntology", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 p=updownplotsreactive()
 print (p)
   dev.off()
 },
 contentType = 'image/png'
 )
  
 ######DoRothEA#####
  
genesreactivedorothea <- reactive({
    switch(input$genesl,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
		   "topConfects" = confectvals()$table
           )
  })
  
regulonreactive <- reactive ({
  switch(input$dorothearegulon,
		"regulon_a" = regulon_a, 
		"regulon_b" = regulon_b,
		"regulon_c" = regulon_c,
		"regulon_d" = regulon_d,
		"regulon_e" = regulon_e,
  )
  })
  
dorothearesult <- eventReactive(input$rundorothea,{
 withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running dorothea"))
  DEsignature <- genesreactivedorothea()
  regulon_user <- regulonreactive()
  viper_regulon = df2regulon(regulon_user)
 
 
if (input$uploadGEO){
if (input$genesl == "eBayes_tvalue" | input$genesl == "TREAT_tvalue"){
DEsignature=DEsignature[!(DEsignature$Gene.symbol==""),]
DEsignature <- subset(DEsignature, ! duplicated(Gene.symbol))
myStatistics <- matrix(DEsignature$logFC, dimnames = list(DEsignature$Gene.symbol, 'logFC') )
myPvalue <- matrix(DEsignature$P.Value, dimnames = list(DEsignature$Gene.symbol, 'P.Value') )
mySignature <- (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
# Reorder
mySignature <- mySignature[order(mySignature, decreasing = T)]
}else if (input$genesl == "topConfects"){
	DEsignature=DEsignature[!(DEsignature$`Gene symbol`==""),]
	DEsignature=DEsignature[which(!is.na(DEsignature$confect)),]
	myStatistics <- matrix(DEsignature$confect, dimnames = list(DEsignature$`Gene symbol`, 'confect') )
	mySignature <- myStatistics[order(myStatistics, decreasing = T),]
	}

}else if (input$uploadSRP){
if (input$genesl == "eBayes_tvalue" | input$genesl == "TREAT_tvalue"){
DEsignature=DEsignature[!(DEsignature$symbol==""),]
DEsignature <- subset(DEsignature, ! duplicated(symbol))
myStatistics <- matrix(DEsignature$logFC, dimnames = list(DEsignature$symbol, 'logFC') )
myPvalue <- matrix(DEsignature$P.Value, dimnames = list(DEsignature$symbol, 'P.Value') )
mySignature <- (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
# Reorder
mySignature <- mySignature[order(mySignature, decreasing = T)]

}else if (input$genesl == "topConfects"){
DEsignature=DEsignature[!(DEsignature$symbol==""),]
DEsignature=DEsignature[which(!is.na(DEsignature$confect)),]
	myStatistics <- matrix(DEsignature$confect, dimnames = list(DEsignature$symbol, 'confect') )
	mySignature <- myStatistics[order(myStatistics, decreasing = T),]
}}else if (input$upload){
if (input$genesl == "eBayes_tvalue" | input$genesl == "TREAT_tvalue"){
# Exclude probes with unknown or duplicated gene symbol
DEsignature <- subset(DEsignature, SYMBOL != "" )
DEsignature <- subset(DEsignature, ! duplicated(SYMBOL))
# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics <- matrix(DEsignature$logFC, dimnames = list(DEsignature$SYMBOL, 'logFC') )
myPvalue <- matrix(DEsignature$P.Value, dimnames = list(DEsignature$SYMBOL, 'P.Value') )
mySignature <- (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
# Reorder
mySignature <- mySignature[order(mySignature, decreasing = T)]
}else if (input$genesl == "topConfects"){
DEsignature=DEsignature[!(DEsignature$SYMBOL==""),]
DEsignature=DEsignature[which(!is.na(DEsignature$confect)),]
	myStatistics <- matrix(DEsignature$confect, dimnames = list(DEsignature$SYMBOL, 'confect') )
	mySignature <- myStatistics[order(myStatistics, decreasing = T),]
}
}
mrs <- msviper(ges = mySignature, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)
mrs 
})})
 
 
dorotheatable <- reactive ({
TF_activities <- data.frame(Regulon = names(dorothearesult()$es$nes),
                            Size = dorothearesult()$es$size[ names(dorothearesult()$es$nes) ], 
                            NES = dorothearesult()$es$nes, 
                            p.value = dorothearesult()$es$p.value, 
                            FDR = p.adjust(dorothearesult()$es$p.value, method = 'fdr'))
TF_activitiesTidy <- TF_activities %>%
  as_tibble() %>%
  arrange(FDR)	
TF_activitiesTidy 
}) 
 
output$dorothea <-renderDT({
  dorotheatable()
  }, filter = 'top',
    rownames = FALSE)


output$dorotheaDownload <- downloadHandler(
    filename = function() {
      paste("dorothea-result", ".csv", sep = "")
    },
	content = function(file) {
	res=dorotheatable()

    write.csv(res, file, row.names = FALSE)  
    }
  )

### dorothea summary	

dorotheasummary <- reactive ({
mrs1 <- ledge(dorothearesult())
# See the results
summary(mrs1)
}) 


output$genesummary <- renderTable({
dorotheasummary()
})

output$Plotdorothea <- renderPlot({
if (input$rundorothea){
plot(dorothearesult(), cex = .7)
}
})

output$dorotheaplotDownload <- downloadHandler(
	filename = function(){
		paste("dorothea-Plot", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
p=plot(dorothearesult(), cex = .7)
print (p)
dev.off()
 },
 contentType = 'image/png'
 )

#####Shadow analysis

shadowreactive <- eventReactive(input$runshadow,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running shadow analysis"))
 mrshadow <- viper::shadow(dorothearesult(), regulators = as.numeric(input$number), verbose = FALSE)
 summary(mrshadow)
})})

output$shadowanalysis <- renderTable({
shadowreactive()$msviper.results
})

output$shadowpairssummary <- renderTable({

shadowreactive()$Shadow.pairs
})

###Cerno analysis

###
##### gene set reactive 
 genesetreactive <- reactive({
    switch(input$cernogene,
            "Gene Ontology Biological Process (MSigDB Filtered)" = GOBP.sel,
           "Gene Ontology Biological Process (Full)" = GOBPfull.sel,
           "Biocarta" = BioCARTA.sel,
           "Gene Ontology Molecular Function (MSigDB Filtered)" = GOMF.sel,
           "Gene Ontology Molecular Function (Full)" = GOMFfull.sel,
           "Gene Ontology Cellular Compartment (MSigDB Filtered)" = GOCC.sel,
           "Gene Ontology Cellular Compartment (Full)" = GOCCfull.sel,
           "Human Phenotype Ontology" = HPO.sel,
           "KEGG" = kegg.sel,
           "Pathway Interaction Database" =PID.sel,
           "Wikipathways" = Wiki.sel,
           "MSigDB Chemical and Genetic Perturbations" = CGP.sel,
           "MSigDB Computational Genesets" = CM.sel,
           "MSigDB Oncogenic Signature Genesets" = Onc.sel,
           "MSigDB Immunologic signature Genesets" = Imm.sel,
           "MSigDB Cell Types" = CellType.sel,
           "Reactome" = Reactome.sel,
           "Hallmark" = Hallmark.sel)
  })

####Ontology
cernoreactive <- eventReactive(input$runcerno,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running cerno"))
	#if (input$cernogene == "Gene_Ontology" | input$cernogene == 'Reactome'| input$cernogene == "Hallmark"){
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	CERNO.res <- tmodLimmaTest(voomreactive(), 
                                voomreactive()$genes$Gene.symbol, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=msig[genesetreactive()])
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	CERNO.res <- tmodLimmaTest(voomreactive(), 
                                voomreactive()$genes$`Gene symbol`, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=msig[genesetreactive()])}
}else if (input$uploadSRP){
CERNO.res <- tmodLimmaTest(voomreactive(), 
                                voomrun()$genes$symbol, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=msig[genesetreactive()])
}else{
CERNO.res <- tmodLimmaTest(voomreactive(), 
                                voomrun()$genes$SYMBOL, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=msig[genesetreactive()])}

CERNO.res
           	})})

output$cernoanalysis <- renderDT({ 
if (input$uploadSRP) {
con2vscon1 = paste(sampleTableSRP()$groups[2], sampleTableSRP()$groups[1], sep="-") 
res=cernoreactive()[[con2vscon1]]
}else if (input$uploadGEO){
 con2vscon1 = paste(geovisualize()$groups[2], geovisualize()$groups[1], sep="-")
res=cernoreactive()[[con2vscon1]]
}else{
con2vscon1 = paste(as.character((input$Group2)), as.character((input$Group1)), sep="-") 
res=cernoreactive()[[con2vscon1]]               
}
(res)
})

output$cernoDownload <- downloadHandler(
    filename = function() {
      paste("cerno-result", ".csv", sep = "")
    },
content = function(file) {
	res=cernoreactive()
  write.csv(res, file, row.names = FALSE)  
    }
  )
  
piereactive <-reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	pie1 <- tmodLimmaDecideTests(voomreactive(), voomreactive()$genes$Gene.symbol, mset=msig[genesetreactive()])
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	pie1 <- tmodLimmaDecideTests(voomreactive(), voomreactive()$genes$`Gene symbol`, mset=msig[genesetreactive()])
}
}else if (input$uploadSRP){
pie1 <- tmodLimmaDecideTests(voomreactive(), voomrun()$genes$symbol, mset=msig[genesetreactive()])

}else{
pie1 <- tmodLimmaDecideTests(voomreactive(), voomrun()$genes$SYMBOL, mset=msig[genesetreactive()])}
pie1

})


output$Panelplot <- renderPlot({

#Make a panel plot of the enrichment results
tmodPanelPlot(cernoreactive(), pie=piereactive(), text.cex=as.numeric(input$textsize), pie.style="rug", clust = "qval")
})


output$PanelplotDownload <- downloadHandler(
	filename = function(){
		paste("Panelplot", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 #Make a panel plot of the enrichment results
 
p=tmodPanelPlot(cernoreactive(), pie=piereactive(), text.cex=as.numeric(input$textsize), pie.style="rug", clust = "qval")
print (p)
dev.off()
 },
 contentType = 'image/png'
 )
 

###Hypergeometric analysis

##### gene set reactive 
genesreactivehyper <- reactive({
    switch(input$geneshyper,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
		   "topConfects" = confectvals()$table
           )
  })

genesetreactivehyper <- reactive({
    switch(input$hypergene,
           "Gene Ontology Biological Process (MSigDB Filtered)" = GOBP.sel,
           "Gene Ontology Biological Process (Full)" = GOBPfull.sel,
           "Biocarta" = BioCARTA.sel,
           "Gene Ontology Molecular Function (MSigDB Filtered)" = GOMF.sel,
           "Gene Ontology Molecular Function (Full)" = GOMFfull.sel,
           "Gene Ontology Cellular Compartment (MSigDB Filtered)" = GOCC.sel,
           "Gene Ontology Cellular Compartment (Full)" = GOCCfull.sel,
           "Human Phenotype Ontology" = HPO.sel,
           "KEGG" = kegg.sel,
           "Pathway Interaction Database" =PID.sel,
           "Wikipathways" = Wiki.sel,
           "MSigDB Chemical and Genetic Perturbations" = CGP.sel,
           "MSigDB Computational Genesets" = CM.sel,
           "MSigDB Oncogenic Signature Genesets" = Onc.sel,
           "MSigDB Immunologic signature Genesets" = Imm.sel,
           "MSigDB Cell Types" = CellType.sel,
           "Reactome" = Reactome.sel,
           "Hallmark" = Hallmark.sel
		   )
  })

hyperreactive <- eventReactive(input$runhyper,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running hypergeometric"))
DEsignature <- genesreactivehyper()

if (input$geneshyper == "eBayes_tvalue" | input$geneshyper == "TREAT_tvalue" ){
DEsignature_down <- DEsignature[DEsignature$adj.P.Val < input$pvalcutoff3 & (DEsignature$logFC) > input$FC2,]
DEsignature_up <- DEsignature[DEsignature$adj.P.Val < input$pvalcutoff3 & (DEsignature$logFC) > -as.numeric(input$FC2),]
}else if (input$geneshyper == "topConfects")
{
DEsignature_down <- DEsignature[which(DEsignature$confect < -as.numeric(input$FC2)), ]
DEsignature_up <- DEsignature[which(DEsignature$confect > input$FC2), ]	
}

if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(DEsignature))
	{
	DEsignature=DEsignature[!(DEsignature$Gene.symbol==""),]
	DEsignature_down = DEsignature_down[!(DEsignature_down$Gene.symbol==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$Gene.symbol==""),]
	down_genes <- as.character(DEsignature_down$Gene.symbol)
	up_genes <- as.character(DEsignature_up$Gene.symbol)
	all_genes <- as.character(DEsignature$Gene.symbol)

	}else if ("Gene symbol" %in% colnames(DEsignature)){
	DEsignature=DEsignature[!(DEsignature$`Gene symbol`==""),]
	DEsignature_down = DEsignature_down[!(DEsignature_down$`Gene symbol`==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$`Gene symbol`==""),]
	down_genes <- as.character(DEsignature_down$`Gene symbol`)
	up_genes <- as.character(DEsignature_up$`Gene symbol`)
	all_genes <- as.character(DEsignature$`Gene symbol`)

		}}else if (input$uploadSRP){
		DEsignature=DEsignature[!(DEsignature$symbol ==""),]
	DEsignature_down = DEsignature_down[!(DEsignature_down$symbol==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$symbol==""),]
	down_genes <- as.character(DEsignature_down$symbol)
	up_genes <- as.character(DEsignature_up$symbol)

	all_genes <- as.character(DEsignature$symbol)
	
		}else{
		DEsignature=DEsignature[!(DEsignature$SYMBOL ==""),]
DEsignature_down = DEsignature_down[!(DEsignature_down$SYMBOL==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$SYMBOL==""),]
	down_genes <- as.character(DEsignature_down$SYMBOL)
	up_genes <- as.character(DEsignature_up$SYMBOL)

all_genes <- as.character(DEsignature$SYMBOL)}

tmodHG_result_down <- tmodHGtest(down_genes, all_genes, mset = msig[genesetreactivehyper()])
tmodHG_result_up <- tmodHGtest(up_genes, all_genes, mset = msig[genesetreactivehyper()])
list(tmodHG_result_up1 = tmodHG_result_up, tmodHG_result_down1=tmodHG_result_down)
})})

output$hyperanalysis <- renderDT({
hyperreactive()$tmodHG_result_up1
})

output$hyperDownload <- downloadHandler(
    filename = function() {
      paste("hypergeometric_Up-result", ".csv", sep = "")
    },
content = function(file) {
	res=hyperreactive()$tmodHG_result_up1
  write.csv(res, file, row.names = FALSE)  
    }
  )

output$hyperanalysisdown <- renderDT({
hyperreactive()$tmodHG_result_down1
})

output$hyperdownDownload <- downloadHandler(
    filename = function() {
      paste("hypergeometric_down-result", ".csv", sep = "")
    },
content = function(file) {
	res=hyperreactive()$tmodHG_result_down1
  write.csv(res, file, row.names = FALSE)  
    }
  ) 

####Cerno TF analysis

ECC_reactive <- reactive({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	ECC.GENES <- as.data.frame(voomreactive()$genes$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	ECC.GENES <- as.data.frame(voomreactive()$genes$`Gene symbol`)
}}else if (input$uploadSRP){
ECC.GENES <- as.data.frame(voomrun()$genes$symbol)
}else{
ECC.GENES <- as.data.frame(voomrun()$genes$SYMBOL)}
	names(ECC.GENES)[1] <- "ID"
	msetECC <- tmod
	msetECC$GENES2MODULES <- NULL
	msetECC$GENES <- ECC.GENES
	msetECC$MODULES2GENES <- ECC.MODULES2GENES
	msetECC$MODULES <- ECC.MODULES
	colnames(msetECC[["MODULES"]])[1]="ID"
	msetECC
	print (msetECC)
})

ReMap_reactive <- reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	ReMap.GENES <- as.data.frame(voomreactive()$genes$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	ReMap.GENES <- as.data.frame(voomreactive()$genes$`Gene symbol`)
}}else if (input$uploadSRP){
ReMap.GENES <- as.data.frame(voomrun()$genes$symbol)
}else{
ReMap.GENES <- as.data.frame(voomrun()$genes$SYMBOL)}

names(ReMap.GENES)[1] <- "ID"
msetReMap <- tmod
msetReMap$GENES2MODULES <- NULL
msetReMap$GENES <- ReMap.GENES
msetReMap$MODULES2GENES <- ReMap.MODULES2GENES
msetReMap$MODULES <- ReMap.MODULES
colnames(msetReMap[["MODULES"]])[1]="ID"
msetReMap
})

TRRUST_reactive <- reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	TRRUST.GENES <- as.data.frame(voomreactive()$genes$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	TRRUST.GENES <- as.data.frame(voomreactive()$genes$`Gene symbol`)
}}else if (input$uploadSRP){
TRRUST.GENES <- as.data.frame(voomrun()$genes$symbol)
}else{
TRRUST.GENES <- as.data.frame(voomrun()$genes$SYMBOL)}
names(TRRUST.GENES)[1] <- "ID"
msetTRRUST <- tmod
msetTRRUST$GENES2MODULES <- NULL
msetTRRUST$GENES <- TRRUST.GENES
msetTRRUST$MODULES2GENES <- TRRUST.MODULES2GENES
msetTRRUST$MODULES <- TRRUST.MODULES
colnames(msetTRRUST[["MODULES"]])[1]="ID"
msetTRRUST
})

TFJ_reactive <- reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	TFJ.GENES <- as.data.frame(voomreactive()$genes$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	TFJ.GENES <- as.data.frame(voomreactive()$genes$`Gene symbol`)
}}else if (input$uploadSRP){
TFJ.GENES <- as.data.frame(voomrun()$genes$symbol)
}else{
TFJ.GENES <- as.data.frame(voomrun()$genes$SYMBOL)}
names(TFJ.GENES)[1] <- "ID"

msetTFJ <- tmod
msetTFJ$GENES2MODULES <- NULL
msetTFJ$GENES <- TFJ.GENES
msetTFJ$MODULES2GENES <- TFJ.MODULES2GENES
msetTFJ$MODULES <- TFJ.MODULES
colnames(msetTFJ[["MODULES"]])[1]="ID"
msetTFJ
})


miRTarBase_reactive <-reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	miRTarBase.GENES <- as.data.frame(voomreactive()$genes$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	miRTarBase.GENES <- as.data.frame(voomreactive()$genes$`Gene symbol`)
}}else if (input$uploadSRP){
miRTarBase.GENES <- as.data.frame(voomrun()$genes$symbol)
}else{
miRTarBase.GENES <- as.data.frame(voomrun()$genes$SYMBOL)}
names(miRTarBase.GENES)[1] <- "ID"

msetmiRTarBase <- tmod
msetmiRTarBase$GENES2MODULES <- NULL
msetmiRTarBase$GENES <- miRTarBase.GENES
msetmiRTarBase$MODULES2GENES <- miRTarBase.MODULES2GENES
msetmiRTarBase$MODULES <- miRTarBase.MODULES
colnames(msetTFJ[["MODULES"]])[1]="ID"
msetmiRTarBase
})

genesetTFcerno <- reactive({
    switch(input$cernoTFgene,
           "ENCODE/ChEA Consensus (Enrichr)" = ECC_reactive(),
           "ReMap ChIP-Seq" = ReMap_reactive(),
		   "TRRUST" = TRRUST_reactive(),
		   "TRANSFAC/JASPAR PWMs (Enrichr)" = TFJ_reactive(),
		   "Gene Transcription Regulation Database (GTRD v20.06)" = msig[GTRD.sel],
		   "miRTarBase 2017" = miRTarBase_reactive(),
		   "miRDB v6.0" = msig[MIRDB.sel]
		   )
  })

cernoTFreactive <- eventReactive(input$runTFcerno,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running TF cerno"))
	#if (input$cernogene == "Gene_Ontology" | input$cernogene == 'Reactome'| input$cernogene == "Hallmark"){

if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	CERNOTF.res <- tmodLimmaTest(voomreactive(), 
                                voomreactive()$genes$Gene.symbol, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=genesetTFcerno())
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	CERNOTF.res <- tmodLimmaTest(voomreactive(), 
                                voomreactive()$genes$`Gene symbol`, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=genesetTFcerno())}
}else if (input$uploadSRP){
print (genesetTFcerno())
print (ECC_reactive())
CERNOTF.res <- tmodLimmaTest(voomreactive(), 
                                voomrun()$genes$symbol, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=genesetTFcerno())
}else{
CERNOTF.res <- tmodLimmaTest(voomreactive(), 
                                voomrun()$genes$SYMBOL, 
                                sort.by = "msd", 
                                tmodFunc = tmodCERNOtest, 
                                mset=genesetTFcerno())}
			
CERNOTF.res
           	})})

output$cernoTFanalysis <- renderDT({    

if (input$uploadSRP) {
con2vscon1 = paste(sampleTableSRP()$groups[2], sampleTableSRP()$groups[1], sep="-") 
res=cernoTFreactive()[[con2vscon1]]
}else if (input$uploadGEO){
 con2vscon1 = paste(geovisualize()$groups[2], geovisualize()$groups[1], sep="-")
res=cernoTFreactive()[[con2vscon1]]
}else{
con2vscon1 = paste(as.character((input$Group2)), as.character((input$Group1)), sep="-") 
res=cernoTFreactive()[[con2vscon1]]               
}

res
})

output$cernoTFDownload <- downloadHandler(
    filename = function() {
      paste("cernoTF-result", ".csv", sep = "")
    },
content = function(file) {
	res=cernoTFreactive()
  write.csv(res, file, row.names = FALSE)  
    }
  )
 

piereactiveTF <-reactive ({
if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(voomreactive()$genes))
	{
	pie1 <- tmodLimmaDecideTests(voomreactive(), voomreactive()$genes$Gene.symbol, mset=genesetTFcerno())
	}else if ("Gene symbol" %in% colnames(voomreactive()$genes)){
	pie1 <- tmodLimmaDecideTests(voomreactive(), voomreactive()$genes$`Gene symbol`, mset=genesetTFcerno())
}
}else if (input$uploadSRP){
pie1 <- tmodLimmaDecideTests(voomreactive(), voomrun()$genes$symbol, mset=genesetTFcerno())

}else{
pie1 <- tmodLimmaDecideTests(voomreactive(), voomrun()$genes$SYMBOL, mset=genesetTFcerno())}
pie1

})
 
output$PanelplotTF <- renderPlot({

#Make a panel plot of the enrichment results
tmodPanelPlot(cernoTFreactive(), pie=piereactiveTF(), text.cex=as.numeric(input$textsize2), pie.style="rug", clust = "qval")
})


output$PanelplotTFDownload <-downloadHandler(
	filename = function(){
		paste("PanelplotTF", ".png", sep="")
 },
 content = function(file){
 png(file,height=780,width=780)
 #Make a panel plot of the enrichment results
 
p=tmodPanelPlot(cernoTFreactive(), pie=piereactiveTF(), text.cex=as.numeric(input$textsize2), pie.style="rug", clust = "qval")
print (p)
dev.off()
 },
 contentType = 'image/png'
 )


### TF Hypergeometric analysis

##### gene set reactive 
genesetTFhyper <- reactive({
    switch(input$hyperTFgene,
           "ENCODE/ChEA Consensus (Enrichr)" = ECC_reactive(),
           "ReMap ChIP-Seq" = ReMap_reactive(),
		   "TRRUST" = TRRUST_reactive(),
		   "TRANSFAC/JASPAR PWMs (Enrichr)" = TFJ_reactive(),
		   "Gene Transcription Regulation Database (GTRD v20.06)" = msig[GTRD.sel],
		   "miRTarBase 2017" = miRTarBase_reactive(),
		   "miRDB v6.0" = msig[MIRDB.sel]
		   )
  })
  
genesreactivehyperTF <- reactive({
    switch(input$geneshyperTF,
           "eBayes_tvalue" = topTable(voomreactive(), adjust.method="BH", sort.by = "t", n = Inf),
           "TREAT_tvalue" = topTreat(treatvals(),sort.by="t", coef=1, n=Inf),
		   "topConfects" = confectvals()$table
           )
  })

hyperTFreactive <- eventReactive(input$runTFhyper,{
withProgress(message = '', value = 0,{
	incProgress(.5, detail = paste("running TF hypergeometric"))
DEsignature <- genesreactivehyperTF()

if (input$geneshyperTF == "eBayes_tvalue" | input$geneshyperTF == "TREAT_tvalue" ){
DEsignature_down <- DEsignature[DEsignature$adj.P.Val < input$pvalcutoff4 & (DEsignature$logFC) > input$FC4,]
DEsignature_up <- DEsignature[DEsignature$adj.P.Val < input$pvalcutoff4 & (DEsignature$logFC) > -as.numeric(input$FC4),]
}else if (input$geneshyperTF == "topConfects")
{
DEsignature_down <- DEsignature[which(DEsignature$confect < -as.numeric(input$FC4)), ]
DEsignature_up <- DEsignature[which(DEsignature$confect > input$FC4), ]	
}


if(input$uploadGEO){
	if ("Gene.symbol" %in% colnames(DEsignature))
	{
	DEsignature=DEsignature[!(DEsignature$Gene.symbol==""),]
	
	DEsignature_down = DEsignature_down[!(DEsignature_down$Gene.symbol==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$Gene.symbol==""),]
	down_genes <- as.character(DEsignature_down$Gene.symbol)
	up_genes <- as.character(DEsignature_up$Gene.symbol)
	all_genes <- as.character(DEsignature$Gene.symbol)
	}else if ("Gene symbol" %in% colnames(DEsignature)){
	DEsignature=DEsignature[!(DEsignature$`Gene symbol`==""),]
	DEsignature_down = DEsignature_down[!(DEsignature_down$`Gene symbol`==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$`Gene symbol`==""),]
	down_genes <- as.character(DEsignature_down$`Gene symbol`)
	up_genes <- as.character(DEsignature_up$`Gene symbol`)
	all_genes <- as.character(DEsignature$`Gene symbol`)
		}}else if (input$uploadSRP){
		DEsignature=DEsignature[!(DEsignature$symbol==""),]
	DEsignature_down = DEsignature_down[!(DEsignature_down$symbol==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$symbol==""),]
	down_genes <- as.character(DEsignature_down$symbol)
	up_genes <- as.character(DEsignature_up$symbol)
all_genes <- as.character(DEsignature$symbol)
	
		}else{
		DEsignature=DEsignature[!(DEsignature$SYMBOL==""),]
	
	DEsignature_down = DEsignature_down[!(DEsignature_down$SYMBOL==""),]
	DEsignature_up = DEsignature_up[!(DEsignature_up$SYMBOL==""),]
	down_genes <- as.character(DEsignature_down$SYMBOL)
	up_genes <- as.character(DEsignature_up$SYMBOL)
	all_genes <- as.character(DEsignature$SYMBOL)}
tmodHG_result_down <- tmodHGtest(down_genes, all_genes, mset = genesetTFhyper())
tmodHG_result_up <- tmodHGtest(up_genes, all_genes, mset = genesetTFhyper())
list(tmodHG_result_up1 = tmodHG_result_up, tmodHG_result_down1=tmodHG_result_down)
})})

output$hyperTFanalysis <- renderDT({
hyperTFreactive()$tmodHG_result_up1
})

output$hyperTFanalysisdown <- renderDT({
hyperTFreactive()$tmodHG_result_down1
})

output$hyperTFDownload <- downloadHandler(
    filename = function() {
      paste("hypergeometric_TFUp-result", ".csv", sep = "")
    },
content = function(file) {
	res=hyperTFreactive()$tmodHG_result_up1
  write.csv(res, file, row.names = FALSE)  
    }
  )


output$hyperTFdownDownload <- downloadHandler(
    filename = function() {
      paste("hypergeometric_TFdown-result", ".csv", sep = "")
    },
content = function(file) {
	res=hyperTFreactive()$tmodHG_result_down1
  write.csv(res, file, row.names = FALSE)  
    }
  ) 
 
  
  }
  