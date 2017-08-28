BedToBigWig <- R6Class(
    classname = "BedToBigWig",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, bedInput = NULL, bwOutput = NULL, toWig = FALSE, editable=FALSE){
            super$initialize("BedToBigWig",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            
            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }
            if(toWig){
                sfx<-".wig"
            }else{
                sfx<-".bw"
            }
            
            if(is.null(bwOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["bwOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),sfx))
                }
            }else{
                private$paramlist[["bwOutput"]] <- bwOutput;
            }
            
            private$paramlist[["toWig"]] <- toWig
            print(private$paramlist)
            private$paramValidation()
        }
    ),
    
    private = list(
        processing = function(){
            genome <- Seqinfo(genome = "hg19")
            bedranges <- import(private$paramlist[["bedInput"]], genome = genome)
            cov <- coverage(bedranges)
            ans <- GRanges(cov)
            ans <- subset(ans, score > 0)
            if(private$paramlist[["toWig"]]){
                export.wig(ans,private$paramlist[["bwOutput"]])
            }else{
                export.bw(ans,private$paramlist[["bwOutput"]])
            }
            
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop(paste("bedInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileCreatable(private$paramlist[["bwOutput"]]);
        }
    )
    
    
)




atacBedToBigWig <- function(atacProc, bedInput = NULL, bwOutput = NULL, toWig = FALSE){
    atacproc <- BedToBigWig$new(atacProc = atacProc, bedInput = bedInput, bwOutput = bwOutput, toWig = toWig)
    atacproc$process()
    return(atacproc)
}
