#!/usr/bin/env Rscript
# date: 1/9/2023 by Hongjian Jin @ St Jude Children's Research Hospital
library(optparse)
option_list <- list(
  make_option(c("-d", "--dataDir"), type="character",default=".", 
    help="character. fastq_screen output directory ")
  ,make_option(c("-o", "--output"), type="character", default=NA, 
    help="character. prefix for output filename")
)

parser <- OptionParser(usage="%prog -d <datadir> -o <outPrefix> ",
                        description = "Generate summary table for fastq_screen outputs",
                        option_list = option_list,
                        epilogue =paste( "Examples:",
            " summarize_fastq_screen.R -d fastq_screen -o YANG2_287210_CUTTAG",
            "\n ", sep="\n\n"))
args<-NA
result<-tryCatch({
    args <- parse_args(parser, positional_arguments = TRUE)  # TRUE = c(0, Inf), FALSE,1,  c(1,2)
    }, warning=function(w){
        message(w)
        cat("\n\n")
    }, error=function(e){
        message(e)
        cat("\n\n")
 })

opt <- args$options
if (any(is.na(args))|| any(is.na(opt)) ) {
    print_help(parser)
    quit("no")
}

#########################
#########################

dataDir<- normalizePath(opt$dataDir) 
prefix<- opt$output

filelist <- dir(dataDir, pattern = ".*_screen.txt$",full.names=TRUE,recursive=FALSE)   

res <- NULL
rn <- NULL
for(inFile in filelist) {
    cat("\n", basename(inFile))
    dat <- read.table(inFile, sep="\t",header=TRUE,fill=TRUE,stringsAsFactors = FALSE, quote="",row.names=1 ,check.names=F,comment.char = "" ,skip=1)
    used <- dat[c("Human","Mouse","Yeast","Ecoli","Vectors"),1:11 ]
    used$mapped<-  rowSums(used[,c(5,7,9,11)])
    tmp <- data.frame(t(used$mapped))
    res <- rbind(res,tmp )
    rn <- c(rn, gsub("_screen.txt","", basename(inFile)))    
}

colnames(res) <- c("Human","Mouse","Yeast","Ecoli","Vectors")
rownames(res) <- rn

outFile <- paste0(prefix,"_fastq_screen_table.txt")
write.table(cbind(ID=rownames(res),res), outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
cat("\n[",basename(outFile),"[saved]")

cat("\n\nCheers!\n\n")
garbage<- gc()
q("no")

######################################################
#update log
######################################################
# 1/10/2023, first version