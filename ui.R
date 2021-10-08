

source(paste(getwd(),'global.R',sep="/"))

header <- dashboardHeader(
  title = "DEGenR"
  #titleWidth = 250
)
header$children[[3]]$children[[3]] <- div(tags$img(src='', align="right", height='50px'))

sidebar <- dashboardSidebar( 
		width =250,
		sidebarMenu(id = "sidebarmenu",
		menuItem("DEGenR Introduction", icon = icon("user"),
                       menuSubItem("Introduction", tabName = "intro")),
		menuItem("Data Upload", tabName = "Manual", icon = icon("dashboard"),
			menuItem("Data from recount2", tabName = "", icon = icon("dashboard"),
			  menuSubItem("Select RNA-seq data", tabName="SRPdataset"),
			  menuSubItem("Data Summary", tabName="groupassign_SRP")
			  ),
			  menuItem("Data from GEO", tabName = "", icon = icon("dashboard"),
			  menuSubItem("Select microarray expression data", tabName="GEOdataset"),
			  #menuSubItem("Sample Clustering", tabName="dataCluster"),
			  menuSubItem("Data Summary", tabName="groupassign")
			  ),
              menuItem("Count Data Upload", tabName = "Manual", icon = icon("dashboard"),
			  menuSubItem("File Upload", tabName="dataInputCounts"),
			  menuSubItem("Data Summary", tabName="dataSummary")
			  )             
			  ),
              menuItem("Sample Contrast", tabName="limmavoom", icon = icon("dashboard")),
			  menuItem("Differential Gene Expression", icon = icon("dashboard"),
						menuSubItem("eBayes", tabName="DEGs"),
                       menuSubItem("TREAT", tabName="treatmethod"),
                       menuSubItem("topConfects", tabName="confect")
			  ),
			  menuItem("Ontology Enrichment Analysis", tabName="enrichment", icon = icon("chart-bar"),
						menuSubItem("camera", tabName="camera"),
                       menuSubItem("fGSEA", tabName="fgsea"),
					   menuSubItem("CERNO", tabName="Cerno"),
					   menuItem("Enrichr","", icon =icon("chart-bar"),
					   menuSubItem("Enrichr Ontology", tabName="enrichrontology"),
					   menuSubItem("Enrichr Plots", tabName="enrichrontologyplots")),
					   menuSubItem("Hypergeometric", tabName="Hypergeometric")
			  ),
			  menuItem("TF Enrichment Analysis", tabName="TFenrichment", icon = icon("chart-bar"),
			  menuItem("Enrichr","", icon =icon("chart-bar"),
			  menuSubItem("Enrichr TF analysis", tabName="enrichr"),
			  menuSubItem("Enrichr Plots", tabName="enrichrTFplots")),
			  menuItem("DoRothEA", tabName="", icon =icon("chart-bar"),
              menuSubItem("DoRothEA TF analysis", tabName="dorothEA"),
			  menuSubItem("VIPER Shadow analysis", tabName="shadowana")),
			  menuSubItem("fGSEA", tabName="fgseaTF"),
			  menuSubItem("CERNO", tabName="CernoTF"),
              menuSubItem("Hypergeometric", tabName="HypergeometricTF"))
			  
  ))
         


body <- dashboardBody(

	tags$head( 
    tags$link(rel = "stylesheet", type = "text/css", href = "my_style.css"),
	 tags$script(
        HTML("
  $(document).ready(function(){
    resize();
  })
  function resize(){
    var h = window.innerHeight - $('.navbar').height() - 150; // Get dashboardBody height
    $('#box').height(h); 
  }"
    )
      )
   ),
   tags$h4( 
    tags$link(rel = "stylesheet", type = "text/css", href = "my_style.css"),
	
	    ),
	tags$style(type="text/css",
  ".shiny-output-error { visibility: hidden; }",
  ".shiny-output-error:before { visibility: hidden; }"
),
  tabItems(
  ############UPLOAD your own###########
  tabItem(tabName="intro",
            #img(src="cri.png"),
            h2("Introduction"),
			p("We developed DEGenR, an interactive web interface that provides integrated tools for 
  performing differential gene expression, rank-based ontological geneset and pathway enrichment analysis, 
  and transcription factor regulatory analysis from user-uploaded raw read counts as well as microarray and sequencing datasets available at the NCBI Gene Expression Omnibus 
  (GEO) and Sequencing Read Archive (SRA)."),
  h3("Data Upload", style="padding-left: 1em"),
  h4("Data from recount2"),
  p("The recount2 project has processed over 2000 RNA-seq studies in the Sequencing Read Archive (SRA) 
  and other sources, included GTEx and TCGA datasets, using the RNA-seq alignment program Rail-RNA. 
  Entering a SRP dataset here will download the recount2-processed RNA expression and metadata for the 
  dataset you select and prepare it for analysis with DEGenR. For a full list of datasets, 
  see https://jhubiostatistics.shinyapps.io/recount/ and select the accession for the dataset you would like to analyze.
References: Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, Jaffe AE, Langmead B, Leek JT. Reproducible RNA-seq analysis using recount2. Nature Biotechnology, 2017. doi: 10.1038/nbt.3838.
Nellore A, Collado-Torres L, Jaffe AE, Alquicira-Hernández J, Wilks C, Pritt J, Morton J, Leek JT, Langmead B. Rail-RNA: scalable analysis of RNA-seq splicing and coverage. Bioinformatics, 2017. doi: 10.1093/bioinformatics/btw575."),
h4("Data from GEO"),
p("The Gene Expression Omnibus (GEO) housed at NCBI is repository of gene expression data, including 
numerous human microarray gene expression studies. This step uses the R package GEOquery to download 
the expression data and metadata for a user-selected microarray study. Search for datasets to 
analyze at https://www.ncbi.nlm.nih.gov/geo/. Only those datatsets having matrix file can be analyzed."),
h4("Count Data Upload"),
p("Users are required to upload two files."),
p("1. RNA-seq raw count data."),
p("2. Metadata table"),
p("Please check the example files in .data/public_data folder")

			),
   tabItem(tabName="dataInputCounts",
			fluidRow(
			
              box(title = "Upload Data",
                   solidHeader = T, status = "primary",
                   width = 6,
                        helpText("Upload your RNA-seq raw count data here. Format the count data file for your dataset as a comma separated values (.csv) file, 
						with the first column providing the gene info and all subsequent columns providing the 
						raw counts for each sample. The header for each sample column should correspond to the sample name. To test the DEGenR pipeline with an example dataset (GSE76987, left_colon_processed), click the submit button"
                                 ,style="color:black; padding-right:0em; font-size:16px;"),
                       # popify(placement = "bottom", title = "File-input info",
                        fileInput("countFile","Choose CSV File", 
                                  accept=c( "text/csv",
											"text/comma-separated-values,text/plain",
											".csv")),
                                                              
                        helpText("Upload your metadata table here, formatted as a tab-delimited file. The first column header should be “Sample” with the sample names provided, 
						and the second column header should be “Condition” with the sample conditions provided."
                                 ,style="color:black; padding-right:0em;font-size:16px;"),
                        
                        fileInput("metaTab", 
                                  "Choose tab-delimited file",
                                  accept = c('text/tab-separated-values',                                             
                                             'text/tab-separated-values',
                                             '.txt',                                             
                                             '.tsv')
                        ),
						selectInput(inputId = "geneinfo",
						label = "Choose gene ID information:",
						choices = c("ENSEMBL",
							"ENTREZID", 
							"SYMBOL")),
						
                  actionButton("upload", label = "Submit",icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),
            ##Input information summary under Data Input Panel
            
              box(title = "Input data Summary",
                  solidHeader = T, status = "primary",
                   width = 6, 
                  h4("samples in the datasets", style="font-size:20px"),
                           tableOutput("sampleInfo")
                  )
				  )),
  
  ############# GEO Dataset ##############
  
tabItem(tabName="GEOdataset",
			fluidRow(
			
              box(title = "Enter a GEO dataset ID",
                   solidHeader = T, status = "primary",  width=12,
                   textInput(inputId="geodataset", 
                             label="Enter GSE accession from GEO"),
                   actionButton("uploadGEO", label = "Submit",icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	)),
            ##Input information summary under Data Input Panel
            fluidRow(
              box(title = "Input data Summary",
                  solidHeader = T, status = "primary",  width=12,
                  
                 # width = 6, 
					h4("samples in the datasets", style="font-size:20px"),
                           tableOutput("sampleInfo2") 
                  ),
				  
				  )
				  ),

############# SRP Dataset ##############

tabItem(tabName="SRPdataset",
			fluidRow(
              box(title = "Select RNA-seq data",
                   solidHeader = T, status = "primary",  width=12,
                   textInput(inputId="srpdataset", 
                             label="Enter SRP accession from recount2. For GTEx data, enter SRP012682; for TCGA, enter TCGA"),
                   actionButton("uploadSRP", label = "Submit",icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	)),
            ##Input information summary under Data Input Panel
            fluidRow(
              box(title = "Input data Summary",
                  solidHeader = T, status = "primary",  width=12,
                  
                 # width = 6, 
					h4("samples in the datasets", style="font-size:20px"),
                           tableOutput("sampleInfo3") 
                  )
				  )),

				  
    #########################################
    ## Second tab content for data summarization panel
    tabItem(tabName="dataSummary",
            
            ##First 3 box under Data Summarization panel
            fluidRow(
              box(title = "Data summary without filtering",
                  solidHeader = T, status = "primary",  width=12,                               
                  
                  fluidRow(
                    ##Raw Count Summary box under the Data Summarization panel
                    
					box(title = "Raw Count Summary",
                       
                        status = "primary",
                        width = 4,
                        #fluidRow(
                        #  column(4,
                                 div(tableOutput("orgLibsizeNormfactor"),style = "font-size:70%"),
                                 tags$style("#orgLibsizeNormfactor table {border: 1px solid black; float: center; position:relative;}","#orgLibsizeNormfactor th {border: 1px solid black;}","#orgLibsizeNormfactor td {border: 1px solid black;}")
                          ),
						  box(title = "Density Plot of unfiltered gene expression",
                        
                        status = "primary",
                        width = 4,
                        
                          #column(4,
						  plotOutput("plotdensityunfiltered")
                                 ),
								 box(title = "Box Plot of unfiltered gene expression",
                        
                        status = "primary",
                        width = 4,
						  plotOutput("boxplotunfiltered")
						  )))),
		  
				  actionButton("Filter", label = "Filter via EdgeR", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                    
    
    fluidRow(
      box(title = "After filtering via EdgeR",
                  solidHeader = T, status = "primary",  width=12, 
      ## Sample Normalization box plot under Data Summarization panel
	  fluidRow(
	  
	  box(title = "Plot densities after filteration",
          status = "primary",
          width = 4,
          plotOutput("plotdensities")
      ),
      box(title = "Boxplot after filteration",
          status = "primary",
          width = 4,
          plotOutput("sampleBoxplot")
		  ),
	box(title = "MDS plot",
           status = "primary",
          width = 4,
          plotOutput("mdsplot")
		  )	  
		  )))
		  ),
		
tabItem(tabName='groupassign',
box (title ="Group assignment",
		solidHeader = T, status = "primary",                  
                  width = 12,
				helpText('Enter GEO accsion ids separated by comma for both Groups (no space between commas). And name the group accordingly
				e.g: control, diseased',style="color:black; padding-right:0em;font-size:16px;"),
				splitLayout(
					textInput("factor1", "Baseline Group"),
					textInput("nam1", "Define Group")
					),
				  splitLayout(
				  textInput("factor2", "Comparison Group"),
				  textInput("nam2", "Define Group")
                   ),
				  
				actionButton(inputId="groupassignment", label="Submit", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
				helpText('Alternatively, you can enter column name and assign sample names to groups. And name the group accordingly
				e.g: control, diseased. Use either of the two options.',style="color:black; padding-right:0em;font-size:16px;"),
				
				textInput(inputId="colfactor", 
                          label="Column name of the factor" 
                          ),
				splitLayout(
					textInput("factor1geo", "Baseline Group"),
					textInput("namgeo1", "Define Group")
					),
				  splitLayout(
				  textInput("factor2geo", "Comparison Group"),
				  textInput("namgeo2", "Define Group")
                   ),
				actionButton(inputId="groupassignmentgeo", label="Submit", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
	
	fluidRow(
	  	  box(title = "Box Plot",
          status = "primary",
          width = 6,
          plotOutput("boxplot1")
      ),
      box(title = "Expression Density Plot",
          status = "primary",
          width = 6,
          plotOutput("expressiondensityplot")
		  )	  
		  )
),

tabItem(tabName='groupassign_SRP',
box (title ="Group assignment",
		solidHeader = T, status = "primary",                  
                  width = 12,
				  textInput(inputId="colfactor_SRP", 
                                           label="Column name of the factor" 
                                           ),
				  textInput(inputId="factor1_SRP", 
                                           label="Baseline Group" 
                                           ),
				  textInput(inputId="factor2_SRP", 
                                           label="Comparison Group" 
                                           ),

				actionButton(inputId="groupassignment2", label="Submit", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
	
	fluidRow(
	  	  box(title = "Box Plot",
          status = "primary",
          width = 4,
          plotOutput("boxplot2")
      ),
      box(title = "Expression Density Plot",
          status = "primary",
          width = 4,
          plotOutput("expressiondensityplot2")
		  ),
	box(title = "MDS Plot",
          status = "primary",
          width = 4,
          plotOutput("MDSplot2")
		  )
		  )
),
 tabItem(tabName="limmavoom",            
            ## 1st row with 2 boxes under limma-voom tab panel
            fluidRow(
              box(title = "Sample contrast matrix and mean-variance trend",
                  solidHeader = T, status = "primary",                  
                  width = 12,
                  fluidRow(
                     
              box(title = "Create the group for contrast matrix and DEGs",
                        
                        status = "primary",
                        uiOutput("grouplevels"),
                        #tags$style("#voomGroupLevel{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                        p("Please select any two groups for comparison"
                          , style="font-weight: bold"),
                        fluidRow(
                          column(6,
                                 textInput(inputId="Group1", 
                                           label="Baseline Group" 
                                           )
                          ),
                          column(6,
                                 textInput(inputId="Group2", 
                                           label="Comparison group" 
                                           )
                          ),
						 actionButton(inputId="degAnalysis", label="Submit", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4") 
                        )
                   # )
					),
					column(6,
                     box(title = "Plot of fitted microarray linear model",
                         width = NULL,
                         status = "primary",
                         plotOutput("voommeanvariance")
                     )))))),          

	tabItem(tabName="DEGs",            
           
                  fluidRow(
                       #column(8,              
                     ## Estimated dispersion under limma-voom tab panel
                     box(title = "topTable of differentially expressed genes (eBayes method)",
                         width = 12,
                         solidHeader = T, status = "primary",
                         div(DTOutput("reslimma"),style = "font-size:100%"),
						 downloadButton("voomDownload", 
                                        label = "Download") 
                     ),
					
                     ## DEGs Summary to summarize the up- and down-regulated DEGs
                     box(title = "DEG Summary",
                         width = 12,
                         
                         solidHeader = T, status = "primary",
						 helpText('Default cutoff is P-value=0.05', style="color:black; padding-right:0em; font-size:16px;"),
                        # h4(textOutput("voomTestDGEtitle"), align="center" ),
                         tableOutput("summarydegs"),
						 
                                 textInput(inputId="pvalcutoff", 
                                           label="FDR adjusted p-value", 
                                           value="0.05"),
                         
						actionButton(inputId="pvalfilter", label="Submit", icon("paper-plane"),
						style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                     )
                     )
					 
			),

tabItem(tabName="treatmethod",            
            ## 1st row with 2 boxes under limma-voom tab panel
            fluidRow(
              box(title = "TREAT analysis of differentially expressed genes",
                  solidHeader = T, status = "primary",                  
                  width = 12,
				   box(title = "topTable of differentially expressed genes (TREAT method)",
                         width = NULL,
                         status = "primary",
                         div(DTOutput("restreat"), style = "font-size:100%"),
						 downloadButton("treatDownload", 
                                        label = "Download",
                                        class = NULL)
                     ),
					 column(6,
					  box(title = "TREAT analysis",
                         width = NULL,
                         status = "primary",
						 h4("TREAT analysis tests differential gene expression relative to a fold-change threshold",style="font-size:20px"),
						 helpText('Default cutoff is Adj.P-value=0.05', style="color:black; padding-right:0em; font-size:16px;"),
                        
						 textInput(inputId="FC", 
                                           label=HTML("log<sub>2</sub>FC cutoff"), 
                                           value="1.5"),
						textInput(inputId="pvalcutoff2", 
                                           label="FDR adjusted p-value cutoff", 
                                           value="0.05"),
                         
						actionButton(inputId="treatanalysis", label="Treat", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
						 tableOutput("summaryaftertreat"),
                         
                     )),
					 column(6,
					 box(title = "MD-plot after TREAT analysis",
                         width = NULL,
                         status = "primary",
						 #tableOutput("summaryaftertreat"),
						 downloadButton("MDplotDownload", 
                                        label = "Download"),
                         plotOutput("MDplot")
						 
                     ))
					 ))
					 ),

						 
	tabItem(tabName="confect",            
            # # 1st row with 2 boxes under limma-voom tab panel
            fluidRow(
              box(title = "topConfects analysis of differentially expressed genes",
                  solidHeader = T, status = "primary",                  
                  width = 12,
				box(title = "top differentially expressed genes (topConfects method)",
                         width = NULL,
                         status = "primary",
						 h4("topConfects builds on the TREAT method to rank genes by confident effect size (based on the Confidence Interval) at a fixed FDR",style="font-size:20px"),
						 textInput(inputId="fdr", 
                                           label="FDR adjusted p-value", 
                                           value="0.05"),
						actionButton(inputId="runconfect", label="Run topconfects", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                         div(DTOutput("restopconfect"), style = "font-size:100%"),
						 downloadButton("topconfectsDownload", 
                                        label = "Download")						 
))),
				fluidRow(
				column(4,
					 box(title = "Plot topConfects results",
                         width = NULL,
                         status = "primary",
						 #tableOutput("summaryaftertreat"),
                         plotOutput("confectplot"),
						 downloadButton("confectplotDownload", 
                                        label = "Download")
                     )),
					  column(4,
					 box(title = "Compare between eBayes and topConfects DEGs",
                         width = NULL,
                         status = "primary",
						 #tableOutput("summaryaftertreat"),
                         plotOutput("voomvsconfectplot"),
						 downloadButton("voomvsconfectplotDownload", 
                                        label = "Download")
                     )),
					 column(4,
					 box(title = "MD plot of the top n number of genes ranked by eBayes or topConfects",
                         width = NULL,
                         status = "primary",
						 #tableOutput("summaryaftertreat"),
						 textInput(inputId="n", 
                                           label="No. of genes", 
                                           value="500"),
						actionButton(inputId="plotlimma_confect", label="MD plot", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                         plotOutput("MD_limma_confects"),
						 downloadButton("MD_limma_confectsDownload", 
                                        label = "Download")
                     ))
)
),

tabItem(tabName="camera",
	fluidRow(                           
             box(title = "Competitive Geneset Test Accounting for Inter-gene Correlation (camera method)",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 #helpText('Choose Enriched gene sets',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "pathwaysname",
					label = "Select geneset database:",
					choices = c(
							"All Gene_Ontology",
							"Gene Ontology: Biological Process (Full)",
							"Gene Ontology: Cellular Component (Full)",
							"Gene Ontology: Molecular Function (Full)",
							"Human Phenotype Ontology",
							"Reactome",
							"MSigDB Hallmark",
							"BioCarta",
							"KEGG",
							"PID",
							"WikiPathways",
							"MSigDB Chemical and Genetic Perturbations ",
							"MSigDB Computational Genesets",
							"MSigDB Oncogenic Signature Genesets",
							"MSigDB Immunologic signature Genesets",
							"MSigDB Cell Types"
							)),
							actionButton(inputId="pathwayenrichment", label="Enrichment (CAMERA)", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                         div(DTOutput("enrichmentout"), style = "font-size:90%"),
						 downloadButton("enrichmentontologyDownload", 
                                        label = "Download")
                     )),
					 fluidRow(
					 box(title = "Barcode plot",width = NULL,
                         solidHeader = T, status = "primary",
						 h4("Enter Gene-set for plotting Barcode Plot"),
						 textInput(inputId="geneset", 
                                           label="Gene Set",
											value= ""),
						actionButton(inputId="barcode", label="Barcode Plot", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),										   
                         plotOutput("barcodeout"),
						 downloadButton("barcodeDownload", 
                                        label = "Download")
                     )
					 )
					 ), 
					 
tabItem(tabName="fgsea",
	fluidRow(
               box(title = "fast Gene Set Enrichment Analysis (fGSEA)",
                         width = NULL,
                         solidHeader = T, status = "primary",  
					 
					 #column(4,
					  box(title = "Gene ranking and geneset database selection",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 
						helpText('Make sure to run TREAT and/or topConfects in order to use as a ranking method',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "genelist",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue", 
							"topConfects")),
							
					helpText('Choose Enriched gene sets',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "pathwaylist",
                  label = "Choose Enriched gene sets:",
                  choices = c("All Gene_Ontology",
							"Gene Ontology: Biological Process (Full)",
							"Gene Ontology: Cellular Component (Full)",
							"Gene Ontology: Molecular Function (Full)",
							"Human Phenotype Ontology",
							"Reactome",
							"MSigDB Hallmark",
							"BioCarta",
							"KEGG",
							"PID",
							"WikiPathways",
							"MSigDB Chemical and Genetic Perturbations",
							"MSigDB Computational Genesets",
							"MSigDB Oncogenic Signature Genesets",
							"MSigDB Immunologic signature Genesets",
							"MSigDB Cell Types")),		
						actionButton(inputId="runfgsea", label="Run fgsea", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
	                     ),
					box(title = "fGSEA Results",
                         width = NULL,
                         solidHeader = T, status = "primary",
					  div(DTOutput("fgseaout_gobp"), style = "font-size:90%"),
						 downloadButton("fgseaDownload", 
                                        label = "Download")),
					box(title = "fGSEA Plot",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                     plotOutput('fgseaplot'),
					 downloadButton("fgseaplotDownload", 
                                        label = "Download"))
					 ))
					 ),
					 
	tabItem(tabName="fgseaTF",
	fluidRow(
               box(title = "Fast Gene Set TF Enrichment Analysis (fgsea)",
                         width = NULL,
                         solidHeader = T, status = "primary",  
					 
					 #column(4,
					  box(title = "Filter criteria",
                         width = NULL,
                         solidHeader = T, status = "primary",
						helpText('Choose from the gene ranking matrix, be careful to run Treat and/or Topconfects if you choose either of these',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "genelistTF",
                  label = "Choose a gene ranking matrix:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue", 
							"topConfects")),
							
					helpText('Choose Enriched gene sets',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "pathwaylistTF",
                  label = "Choose Enriched gene sets:",
                  choices = c("ENCODE-ChEA Consensus (Enrichr)",
							"ChEA 2016 (Enrichr)",
							"ENCODE 2015 (Enrichr)",
							"ReMap ChIP-seq 2018 Human",
							"TRRUST 2019 Human",
							"ChEA3 Literature ChIP-Seq",
							"TRANSFAC/JASPAR PWMs (Enrichr)",
							"Gene Transcription Regulation Database (GTRD v20.06)",
							"MSigDB Legacy TF targets",
							"miRTarBase 2017",
							"miRNA TargetScan 2017",
							"miRDB v6.0")),		
						actionButton(inputId="runfgseaTF", label="Run fgsea", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
	                     ),
					 div(DTOutput("fgseaout_TF"), style = "font-size:90%"),
						 downloadButton("fgseaTFDownload", 
                                        label = "Download"),
						plotOutput('fgseaplotTF')
					 ))
					 ),
	
					 
	tabItem(tabName="enrichrontology",
	fluidRow(
				box(title = "Ontology enrichment using Enrichr",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    box(title = "Filter criteria",
                         width = NULL,
                         solidHeader = T, status = "primary",
			helpText('Make sure to run TREAT and/or topConfects in order to use as a ranking method',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "genesenrichr",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue",
							"topConfects")),		 
			selectInput(inputId = "databaseOntology",
                  label = "Select geneset database:",
                  choices = c("GO_Biological_Process_2018",
                  "GO_Biological_Process_2017b", 
                  "GO_Molecular_Function_2018",
                  "GO_Molecular_Function_2017b", 
                  "GO_Cellular_Component_2018",
                  "GO_Cellular_Component_2017b",
                  "MSigDB_Hallmark_2020",
                  "Reactome_2016",
                  "BioCarta_2016",
                  "KEGG_2019_Human",
                  "Panther_2016",
                  "WikiPathways_2019_Human",
                  "BioPlanet_2019")),		 
				  textInput(inputId="FC4", 
                                           label=HTML("log<sub>2</sub>FC or confect (topConfects) cutoff"), 
                                           value="0.5"),
					textInput(inputId="pvalenrichonto",
					label="P value cutoff", 
                                           value="0.05"
					),					   
						actionButton(inputId="runenrichrOntology", label="Run enrichr Ontology", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 
					 box(title = "Ontology enrichment on up-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
                         div(DTOutput("up_ontologyenrichr"), style = "font-size:90%; width: 90%"),
						 downloadButton("upenrichrontologyDownload", 
                                        label = "Download")
                     ),
					
					  box(title = "Ontology Enrichment analysis on down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 						
						div(DTOutput("down_ontologyenrichr"), style = "font-size:90%; width: 90%"),
						downloadButton("downenrichrontologyDownload", 
                                        label = "Download")
										)
                      ))                
					 ),
	
	tabItem(tabName="enrichrontologyplots",
	fluidRow(
		box(title = "Filter criteria",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 textInput(inputId="num", 
                                           label="Number of terms to plot", 
                                           value="20"),
						selectInput(inputId = "pval_FDR2",
                  label = "Select P-val or FDR:",
                  choices = c("P.value", 
							"Adjusted.P.value")),
					actionButton(inputId="filterenrichrOntology", label="Plotting", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),
	),
	fluidRow(column(6,
				box(title = "Bar plot of up-regulated terms",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    plotOutput('upplots'),
					downloadButton("upplotsenrichrOntologyDownload", 
                                        label = "Download")
					)),
			column(6,
				box(title = "Bar plot of down-regulated terms",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    plotOutput('downplots'),
					downloadButton("downplotsenrichrOntologyDownload", 
                                        label = "Download")
					))),		
					 
			fluidRow(
			box(title = "Combined scores from up-regulated and down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 helpText("Combined plot will have double the number of terms used in filter crieria"),
                         plotOutput('updownplots'),
						 downloadButton("updownplotsenrichrOntologyDownload", 
                                        label = "Download")
                     ))),
					
					  

	tabItem(tabName="enrichr",
	fluidRow(
				box(title = "TF enrichment using Enrichr",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    box(title = "Filter criteria",
                         width = NULL,
                         solidHeader = T, status = "primary",
			helpText('Make sure to run TREAT and/or topConfects in order to use as a ranking method',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "genes",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue",
							"topConfects")),		 
			selectInput(inputId = "database",
                  label = "Select geneset database:",
                  choices = c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
            "ENCODE_TF_ChIP-seq_2015", 
            "ChEA_2016",
            "TRANSFAC_and_JASPAR_PWMs", 
            "TargetScan_microRNA",
            "ARCHS4_TFs_Coexp",
            "TRRUST_Transcription_Factors_2019",
            "TargetScan_microRNA_2017",
            "miRTarBase_2017"))		 
					 ,	

						textInput(inputId="FC3", 
                                           label=HTML("log<sub>2</sub>FC or confect (topConfects) cutoff"), 
                                           value="0.5"),
                                           textInput(inputId="pvalenrichTF",
					label="P value cutoff", 
                                           value="0.05"
					),
										   
						actionButton(inputId="runenrichr", label="Run enrichr", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 
					 box(title = "TF enrichment on up-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
                         div(DTOutput("up_enrichr"), style = "font-size:90%; width: 90%"),
						 downloadButton("upenrichrDownload", 
                                        label = "Download")
                     ),
					
					  box(title = "TF Enrichment analysis on down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 						
						div(DTOutput("down_enrichr"), style = "font-size:90%; width: 90%"),
						downloadButton("downenrichrDownload", 
                                        label = "Download")
										)
                      ))                
					#,uiOutput("Next_stepdorothEA", align="center")  
					 ),

tabItem(tabName="enrichrTFplots",
	fluidRow(
		box(title = "Filter criteria",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 textInput(inputId="num2", 
                                           label="Number of terms to plot", 
                                           value="20"),
						selectInput(inputId = "pval_FDR",
                  label = "Select P-val or FDR:",
                  choices = c("P.value", 
							"Adjusted.P.value")),
					actionButton(inputId="filterenrichrTF", label="Plotting", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),
	
	),
	fluidRow(column(6,
				box(title = "Bar plot of up-regulated terms",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    plotOutput('upTFplots'),
					downloadButton("upTFplotsDownload", 
                                        label = "Download")
					)),
			column(6,
				box(title = "Bar plot of down-regulated terms",
                         width = NULL,
                         solidHeader = T, status = "primary",  
                    plotOutput('downTFplots'),
					downloadButton("downTFplotsDownload", 
                                        label = "Download")
					
					))),		
					 
			fluidRow(
			box(title = "Combined scores from up-regulated and down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 helpText("Combined plot will have double the number of terms used in filter crieria"),
                         plotOutput('updownTFplots'),
						 downloadButton("updownTFplotsDownload", 
                                        label = "Download")
                     ))),

tabItem(tabName="dorothEA",
	fluidRow(
                  box(title = "TF activity analysis using viper algorithm and DoRothEA regulons",
                         width = NULL,
                         solidHeader = T, status = "primary",
					selectInput(inputId = "dorothearegulon",
                  label = "Select the DoRothEA regulon:",
                  choices = c("regulon_a", 
							"regulon_b",
							"regulon_c",
							"regulon_d",
							"regulon_e")),	 
					
				helpText('Make sure to run TREAT and/or topConfects in order to use as a ranking method',style="color:black; padding-right:0em; font-size:16px;"),
						selectInput(inputId = "genesl",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue",
							"topConfects")),
 						actionButton(inputId="rundorothea", label="Run DoRothEA", icon("paper-plane"), 
		style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
				
					 box(title = "DoRothEA TF Activity Analysis",
                         width = NULL,
                         solidHeader = T, status = "primary",
                         div(DTOutput("dorothea"), style = "font-size:90%"),
						 downloadButton("dorotheaDownload", 
                                        label = "Download")
                     )),
					fluidRow(
					  box(title = "Genes contributing most to these TF activity",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 					
                         tableOutput("genesummary")
						 ),
						box(title = "A graphics representation of the results (msVIPER plot)",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 					
                         plotOutput("Plotdorothea"),
						 downloadButton("dorotheaplotDownload", 
                                        label = "Download")
						 )
					 )),
tabItem(tabName="shadowana",
	fluidRow(
                  box(title = "Shadow analysis",
                         width = NULL,
                         solidHeader = T, status = "primary",
					
				helpText('A regulator may appear to be significantly activated because it may share its regulon 
of its with a activated TF (shadow effect). To account for this shadow analysis is performed which can list shadow pairs.',style="color:black; padding-right:0em; font-size:16px;"),
						textInput(inputId="number", 
                                           label="Number of top regulators", 
                                           value="25"),
 						actionButton(inputId="runshadow", label="Run shadow analysis", icon("paper-plane"), 
		style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
				
					 box(title = "Result of Shadow Analysis",
                         width = NULL,
                         solidHeader = T, status = "primary",
                         tableOutput("shadowanalysis"),
						 downloadButton("shadowDownload", 
                                        label = "Download")
                     ),
					
					  box(title = "Shadow pairs",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 					
                         tableOutput("shadowpairssummary")
						 )
						
					 )),
					 
	tabItem(tabName="Cerno",
	fluidRow(
                  box(title = "Coincident Extreme Ranks in Numerical Observations (CERNO)",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 
					selectInput(inputId = "cernogene",
                  label = "Select geneset database:",
                  choices = c("Hallmark", 
							"Reactome", 
							"Gene Ontology Biological Process (MSigDB Filtered)",
							"Gene Ontology Biological Process (Full)",
							"Biocarta",
							"Gene Ontology Molecular Function (MSigDB Filtered)",
							"Gene Ontology Molecular Function (Full)",
							"Gene Ontology Cellular Compartment (MSigDB Filtered)",
							"Gene Ontology Cellular Compartment (Full)",
							"Human Phenotype Ontology",
							"KEGG",
							"Pathway Interaction Database",
							"Wikipathways",
							"MSigDB Chemical and Genetic Perturbations",
							"MSigDB Computational Genesets",
							"MSigDB Oncogenic Signature Genesets",
							"MSigDB Immunologic signature Genesets",
							"MSigDB Cell Types"
							)),
					textInput(inputId="textsize", 
                    label="Text size for the Panel Plot", 
                              value="0.5"),
						actionButton(inputId="runcerno", label="Run cerno", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 	
					box(title = "Results from CERNO analysis",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("cernoanalysis"),
						 div(DTOutput("cernoanalysis"), style = "font-size:90%"),
							downloadButton("cernoDownload", 
                                        label = "Download")						 
					 ),
					 box(title = "Panel Plot",
                         width = NULL, status = "primary",
						downloadButton("PanelplotDownload", 
                                        label = "Download"),
                         plotOutput("Panelplot")
						 )
					 )),
					 
	tabItem(tabName="Hypergeometric",
	fluidRow(
                  box(title = "Hypergeometric Enrichment",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 
				helpText("Make sure to run TREAT and/or topConfects in order to use as a ranking method"),
				selectInput(inputId = "geneshyper",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue",
							"topConfects")),	 
					selectInput(inputId = "hypergene",
                  label = "Select geneset database:",
                  choices = c("Hallmark", 
							"Reactome", 
							"Gene Ontology Biological Process (MSigDB Filtered)",
							"Gene Ontology Biological Process (Full)",
							"Biocarta",
							"Gene Ontology Molecular Function (MSigDB Filtered)",
							"Gene Ontology Molecular Function (Full)",
							"Gene Ontology Cellular Compartment (MSigDB Filtered)",
							"Gene Ontology Cellular Compartment (Full)",
							"Human Phenotype Ontology",
							"KEGG",
							"Pathway Interaction Database",
							"Wikipathways",
							"MSigDB Chemical and Genetic Perturbations",
							"MSigDB Computational Genesets",
							"MSigDB Oncogenic Signature Genesets",
							"MSigDB Immunologic signature Genesets",
							"MSigDB Cell Types")),	

					textInput(inputId="pvalcutoff3", 
                                           label="FDR adjusted p-value cutoff", 
                                           value="0.05"),
					textInput(inputId="FC2", 
                                           label=HTML("log<sub>2</sub>FC or confect (topConfects) cutoff"), 
                                           value="0.5"),
										   
						actionButton(inputId="runhyper", label="Run Hypergeometric", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 	
					box(title = "Hypergeometric test of up-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("hyperanalysis"),
						 div(DTOutput("hyperanalysis"), style = "font-size:90%"),
							downloadButton("hyperDownload", 
                                        label = "Download")						 
					 ),
					 box(title = "Hypergeometric test of down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("hyperanalysisdown"),
						 div(DTOutput("hyperanalysisdown"), style = "font-size:90%"),
							downloadButton("hyperdownDownload", 
                                        label = "Download")						 
					 )
					 )),
					 
	tabItem(tabName="CernoTF",
	fluidRow(
                  box(title = "TF-gene target enrichment using Coincident Extreme Ranks in Numerical Observations (CERNO)",
                         width = NULL,
                         solidHeader = T, status = "primary",
						 
					selectInput(inputId = "cernoTFgene",
                  label = "Select geneset database:",
                  choices = c("ENCODE/ChEA Consensus (Enrichr)",
							"ReMap ChIP-Seq",
							"TRRUST" ,
							"TRANSFAC/JASPAR PWMs (Enrichr)",
							"Gene Transcription Regulation Database (GTRD v20.06)",
							"miRTarBase 2017" ,
							"miRDB v6.0")),
					textInput(inputId="textsize2", 
                    label="Text size for the Panel Plot", 
                              value="0.5"),							
						actionButton(inputId="runTFcerno", label="Run cerno", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 	
					box(title = "Results from CERNO TF analysis",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("cernoTFanalysis"),
						 div(DTOutput("cernoTFanalysis"), style = "font-size:90%"),
							downloadButton("cernoTFDownload", 
                                        label = "Download")						 
					 ),
					 box(title = "Panel Plot",
                         width = NULL, status = "primary",
						downloadButton("PanelplotTFDownload", 
                                        label = "Download"),
                         plotOutput("PanelplotTF")
						 )
					 )),
	tabItem(tabName="HypergeometricTF",
	fluidRow(
                  box(title = "Hypergeometric TF Enrichment",
                         width = NULL,
                         solidHeader = T, status = "primary",
					helpText("Make sure to run TREAT and/or topConfects in order to use as a ranking method"),
				selectInput(inputId = "geneshyperTF",
                  label = "Select a gene ranking method:",
                  choices = c("eBayes_tvalue", 
							"TREAT_tvalue",
							"topConfects")),	 
					selectInput(inputId = "hyperTFgene",
                  label = "Select geneset database:",
                  choices = c("ENCODE/ChEA Consensus (Enrichr)",
							"ReMap ChIP-Seq",
							"TRRUST" ,
							"TRANSFAC/JASPAR PWMs (Enrichr)",
							"Gene Transcription Regulation Database (GTRD v20.06)",
							"miRTarBase 2017" ,
							"miRDB v6.0")),	

					textInput(inputId="pvalcutoff4", 
                                           label="FDR adjusted p-value cutoff", 
                                           value="0.05"),
					textInput(inputId="FC4", 
                                           label=HTML("log<sub>2</sub>FC or confect (topConfects) cutoff"), 
                                           value="0.5"),
										   
						actionButton(inputId="runTFhyper", label="Run Hypergeometric", icon("paper-plane"), 
    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
	),			 	
					box(title = "Hypergeometric TF enrichment of up-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("hyperTFanalysis"),
						 div(DTOutput("hyperTFanalysis"), style = "font-size:90%"),
							downloadButton("hyperTFDownload", 
                                        label = "Download")						 
					 ),
					 box(title = "Hypergeometric TF enrichment of down-regulated genes",
                         width = NULL,
                         solidHeader = T, status = "primary",	
                         #tableOutput("hyperTFanalysisdown"),
						 div(DTOutput("hyperTFanalysisdown"), style = "font-size:90%"),
							downloadButton("hyperTFdownDownload", 
                                        label = "Download")						 
					 )
					 ))			 
            
    )
    )
	
	ui <- dashboardPage(header, sidebar, body)
	
  























