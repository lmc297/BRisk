

library(shiny)
library(readr)
library(stringr)
library(fitdistrplus)
library(mc2d)
library(ggplot2)


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel(HTML("<i>BRisk")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        fileInput("file1","Choose a BTyper final results file to analyze",
                  multiple=FALSE,accept=c("_final_results.txt")),
        helpText("Final results files are text files that BTyper
                 creates. They have the extension '_final_results.txt'."),
        selectInput("foodmatrix",
                  label="Select a food matrix",
                  choices=c(
                  "Select a food matrix"=0,
                  "Milk, pasteurized fluid"=1, 
                  "Rice, cooked white"=2,
                  "Soup, vegetable"=3),selected=0),
      fileInput("file1m","Upload a bacterial count data file",
                multiple=FALSE,accept=c(".txt")),
      helpText("A count data file is a text (.txt) file that contains B. cereus group isolate counts in log(CFU/g) or log(CFU/mL),
               with one count per line.")),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("hist"),
         htmlOutput("riskText"),
         tableOutput("riskTable")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   get.infile<-function(){
     infile <- input$file1
     if (is.null(infile)){
       return(NULL)}
     else{
       rfile<-read_file(infile$datapath)}
     lines<-strsplit(rfile,"\n")
     pgenes<-lapply(lines,function(x) which(x%in%"panC Clade Name\tClosest Strain\tPercent (%) Identity\tPercent (%) Coverage"))[[1]]
     pts<-c()
     for (i in 1:length(lines)){
       pgroup<-strsplit(lines[[i]][pgenes[[i]][1]+1],split="\t")[[1]][1:2]
       punit<-paste(pgroup[1],pgroup[2],sep="_")
       if (!grepl("Predicted",punit)){
         pts<-append(pts,punit)}
     }
     pts<-gsub("_.*","\n",pts)
     pts<-gsub("\\*","",pts)
     clade<-pts[which(as.character(pts)!="NA\nNA")]
     clade<-gsub("[\r\n]", "", clade)
     anthrax<-c("\tcya\\|","\tlef\\|","\tpagA\\|")
     emetic<-c("\tcesA\\|","\tcesB\\|","\tcesC\\|","\tcesD\\|")
     dia1<-c("\tcytK1\\|","\tcytK2\\|")
     dia2<-c("\thblA\\|","\thblB\\|","\thblC\\|","\thblD\\|")
     vgenes<-rownames(sapply(lines,function(x) which(vapply(x,function(r) grepl("\\|",r)&&!(grepl("rpoB\\|",r))&&!(grepl("\\(",r)),FUN.VALUE = 1)==1)))
     agenes<-vgenes[grep(paste(anthrax,collapse="|"), vgenes)]
     egenes<-vgenes[grep(paste(emetic,collapse="|"), vgenes)]
     d1genes<-vgenes[grep(paste(dia1,collapse="|"),vgenes)]
     d2genes<-vgenes[grep(paste(dia2,collapse="|"),vgenes)]
     print(clade)
     if(length(agenes)>0){
       disease<-'<p style="color:darkred;">high-risk isolate (anthrax)</p>'}
     else if(length(egenes)==4 && clade=="clade3"){
       disease<-'<p style="color:darkred;">high-risk isolate (emetic disease)</p>'}
     else if(length(d1genes)>1||length(d2genes)==4){
       disease<-'<p style="color:darkred;">high-risk isolate (diarrheal disease)</p>'}
     else{
       disease<-'<p style="color:mediumseagreen;">low-risk isolate (no known disease outcome)</p>'}
     return(list(clade,disease))
    }#end get.infile
   
   get.counts<-function(){
     infile<-input$file1m
     if (is.null(infile)){
       return(NULL)}
     if (is.null(input$file1)){
       return(NULL)}
     if (input$foodmatrix==0){
       return(NULL)}
     else{
       mtable<-read.table(infile$datapath,header=FALSE,sep = "\t")
       validate(need(ncol(mtable)==1,"Your count data file needs to have exactly 1 column with B. cereus group counts in log(CFU/g) or log(CFU/mL)."))
       # fit lognormal dist'n per USDA-FSIS 2003 Listeria in deli meats, page 49: https://www.fsis.usda.gov/OPPDE/rdad/FRPubs/97-013F/ListeriaReport.pdf
       if(length(mtable$V1)>1){
         dist0<-fitdist(mtable$V1,distr = "lnorm",method="mle")
         sample0<-rlnorm(n = 1,meanlog = dist0$estimate[1],sdlog = dist0$estimate[2])}
       else{
         sample0<-rlnorm(n = 1,meanlog = mtable$V1,sdlog = 1)}
       return(sample0)
      }
   }# end get.counts
   
   get.growth<-function(clade){
     # using t-distribution for parameters: A Toolbox for Nonlinear Regression in R: The Package nlstools (page 8)
     if (input$foodmatrix==0){
       return(NULL)}
     else{
     food<-input$foodmatrix}
     if(clade=="clade3"&&food==1){
       mumax<-1.161991
       log10nmax<-8.62642
     }#clade3
     else if(clade=="clade3"&&food==2){
       mumax<-2
       log10nmax<-9
     }#clade3
     else if(clade=="clade3"&&food==3){
       mumax<-0.5
       log10nmax<-6
     }#clade3
     # clade4
     else if(clade=="clade4"&&food==1){
       mumax<-1.5
       log10nmax<-5
     }#clade4
     else if(clade=="clade4"&&food==2){
       mumax<-1.9
       log10nmax<-9.5
     }#clade4
     else if(clade=="clade4"&&food==3){
       mumax<-0.3
       log10nmax<-7
     }#clade4
     ######## only adding this for demo version
     else{
        mumax<-0.01
        log10nmax<-1}
     return(list(mumax,log10nmax))
   }#end get.growth
   
   # log10N
   log10N_func = function(t,mumax,LOG10N0,LOG10Nmax){
     ans <- LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) * log(10)/mumax)) *     mumax * t/log(10) + (t > ((LOG10Nmax - LOG10N0) * log(10)/mumax)) *     (LOG10Nmax - LOG10N0)
     return(ans)
   }
   
   
   output$hist<-renderPlot({
     if(input$foodmatrix==0){
       return(NULL)}
     if(is.null(input$file1)){
       return(NULL)}
     if(is.null(input$file1m)){
       return(NULL) }
     virulence<-get.infile()
     cfuvec<-c()
     for (i in 1:100){
     counts<-get.counts()
     print("counts")
     print(counts)
     mumax<-get.growth(virulence[[1]])[[1]]
     log10nmax<-get.growth(virulence[[1]])[[2]]
     if(input$foodmatrix==1){
       # usda table III-12: https://www.fda.gov/downloads/Food/FoodScienceResearch/UCM197330.pdf
       retailtemp<-runif(1,1,5)
       retailstorage<-runif(1,1,3)
       EGRratio<-(5.18/(retailtemp+1.18))^2
       mumaxNEW<-mumax/EGRratio
       logCFU1<-log10N_func(t = retailstorage,mumax = mumaxNEW,LOG10N0 = counts,LOG10Nmax = log10nmax)
       # USDA table III-7: https://www.fda.gov/downloads/Food/FoodScienceResearch/UCM197330.pdf
       bin1<-runif(n = 1,min = 0,max = 32)
       bin2<-runif(n = 1,min = 33,max = 35)
       bin3<-runif(n = 1,min = 36,max = 38)
       bin4<-runif(n = 1,min = 39,max = 41)
       bin5<-runif(n = 1,min = 42,max = 44)
       bin6<-runif(n = 1,min = 45,max = 47)
       bin7<-runif(n = 1,min = 48,max = 50)
       bin8<-runif(n = 1,min = 51,max = 53)
       bin9<-runif(n = 1,min = 54,max = 56)
       bin10<-runif(n = 1,min = 57,max = 59)
       bin11<-runif(n = 1,min = 60,max = 63)
       mybin<-sample(x = c(paste("bin",c(1:11),sep="")),size = 1,prob = c(0.09,0.1,0.25,0.29,0.18,0.05,0.03,0.004,0.005,0.004,0.001),replace = TRUE)
       ftemp<-eval(parse(text=mybin))
       consumertemp<-(ftemp-32)/1.8
       # USDA table III-5: https://www.fda.gov/downloads/Food/FoodScienceResearch/UCM197330.pdf
       consumerstorage<-rpert(n = 1,min = 0.5,mode = c(3,5),max = c(10,15))
       EGRratio2<-(5.18/(consumertemp+1.18))^2
       mumaxNEW2<-mumax/EGRratio2
       logCFU2<-log10N_func(t = consumerstorage,mumax = mumaxNEW2,LOG10N0 = logCFU1,LOG10Nmax = log10nmax)
       print("logCFU2")
       print(logCFU2)
     # USDA table III-3: https://www.fda.gov/downloads/Food/FoodScienceResearch/UCM197330.pdf
       serving.size<-sample(x = c(rep(x = 244,50),rep(245,25),rep(488,20),rep(732,5)),size = 1,replace = TRUE)
       print("serving size")
       print(serving.size)}
     else if (input$foodmatrix==2 ){
       # "expert" opinion: min 15.555
       ftemp<-rpert(n = 1,min = 60,mode = c(135,140),max = 176)
       consumertemp<-(ftemp-32)/1.8
       EGRratio<-(136.18/(consumertemp+1.18))^2
       mumaxNEW<-mumax/EGRratio
       consumerstorage<-rpert(n = 1,min = 0,mode = c(0.000694444,0.1666666),max = 1)
       logCFU2<-log10N_func(t = consumerstorage,mumax = mumaxNEW,LOG10N0 = counts,LOG10Nmax = log10nmax)
       # "expert" opinion: 1/8 c min, 1/4-2 c mode, 4 c max
       serving.size<-rpert(n = 1,min = 30,mode = c(65,500),max = 1000)}
     else if (input$foodmatrix==3){
       ftemp<-rpert(n = 1,min = 90,mode = c(135,140),max = 212)
       consumertemp<-(ftemp-32)/1.8
       EGRratio<-(136.18/(consumertemp+1.18))^2
       mumaxNEW<-mumax/EGRratio
       consumerstorage<-rpert(n = 1,min = 0,mode = c(0.000694444,0.1666666),max = 1)
       logCFU2<-log10N_func(t = consumerstorage,mumax = mumaxNEW,LOG10N0 = counts,LOG10Nmax = log10nmax)
       # "expert" opinion: 1/8 c min, 1/4 c-4 c, max=1 gallon
       serving.size<-rpert(n = 1,min = 29,mode = c(59,946.3),max = 3785.41)}
       realCFU<-10^logCFU2
       CFUsim<-realCFU*serving.size
       cfuvec<-c(cfuvec,CFUsim)
     }# end 1 simulation
     #return(virulence[[2]])
     cfu.df<-data.frame(cfuvec)
     cfu.df$color<-ifelse(test = cfu.df$cfuvec>=100000,yes = "High",no = 
                            ifelse(cfu.df$cfuvec>=100,yes = "Medium",no = "Low"))
     finalhist<-ggplot(data = cfu.df,aes(cfuvec))
     finalhist<-finalhist+
       geom_histogram(data = cfu.df,aes(fill=color),binwidth = 0.1)+
       scale_fill_manual("Risk",values = c("High"="red3","Medium"="darkorange1","Low"="springgreen3"))+
       #geom_histogram(data=subset(cfu.df,color=="High"), fill="red3") +
       #geom_histogram(data=subset(cfu.df,color=="Medium"), fill="darkorange1") +
       #geom_histogram(data=subset(cfu.df,color=="Low"), fill="springgreen3") +
       scale_x_log10()+
       scale_y_continuous()+
       xlab(label = "Colony Forming Units (CFU)")+
       ylab(label = "Servings")+
       ggtitle(label = "Histogram of Colony Forming Units (CFU) per Serving")
     histtable<-ggplot_build(finalhist)$data
     print(histtable)
     return(finalhist)} #end hist
    )
   
   output$riskText<-renderUI({
     if(is.null(input$file1)){
       return(NULL)}
     else{
     infile<-get.infile()
     clade<-infile[[1]]
     disease<-infile[[2]]
     finalstring<-HTML(paste("<h2> Your isolate is a <b>",disease,"</b></h2>",sep = ""))
     return(finalstring)}
   })# end riskText
   
  
   

}

# Run the application 
shinyApp(ui = ui, server = server)

