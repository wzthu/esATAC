PeakQC <-R6Class(
    classname = "PeakQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),bedInput = NULL,editable=FALSE){
            super$initialize("PeakQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            qcbedInput <- qcbedInput[1]
            if(qcbedInput == "DHS"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("DHS");
            }else if(qcbedInput == "blacklist"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("blacklist");
            }else{
                private$paramlist[["qcbedInput"]]<-qcbedInput;
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report.txt"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }
            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            genome <- Seqinfo(genome = .obtainConfigure("genome"))

            inputbed <- import(private$paramlist[["bedInput"]], genome = genome)


            qcbedInput<-import(private$paramlist[["qcbedInput"]], genome = genome)



            qcval=list();

            qcval[["totalInput"]]<-length(inputbed)
            qcval[["qcbedInput"]]<-length(subsetByOverlaps(inputbed, qcbedInput,ignore.strand = TRUE))
            qcval[["qcbedRate"]]<-qcval[["qcbedInput"]]/qcval[["totalInput"]]

            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["qcbedInput"]])){
                stop("qcbedInput is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["qcbedInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","totalInput","qcbedInput","qcbedRate"))
        }
    )


)

#' @name atacPeakQC
#' @aliases atacPeakQC
#' @aliases peakQC
#' @title Quality control for peak overlap 
#' @description 
#' These functions are used to generate fregment distribution plot.  
#' The fourier transform of fregment distribution will be calculated.
#' Strength distribution around period at 10.4bp and 180bp 
#' will be shown in another two plots.
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}, 
#' \code{\link{atacBedUtils}}.
#' @param reportOutput \code{Character} scalar. 
#' The report file path. 
#' @param qcbedInput \code{Character} scalar. 
#' It can be "DHS","blacklist" or 
#' Other quality control BED file input path.
#' @param bedInput \code{Character} scalar. 
#' BED file input path for quality control.
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{peakQC} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}
#' @importFrom Rcpp  evalCpp
#' @importFrom R6  R6Class
#' @importFrom igraph  graph
#' @importFrom igraph vertex.attributes
#' @importFrom igraph `vertex.attributes<-`
#' @importFrom igraph are.connected
#' @importFrom Rcpp sourceCpp
#' @importFrom rJava .jpackage
#' @importFrom rJava .jnew
#' @importFrom rJava .jcall
#' @importFrom rtracklayer import
#' @importFrom rtracklayer export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_ribbon
#' @importFrom DiagrammeR render_graph
#' @importFrom DiagrammeR to_igraph
#' @importFrom DiagrammeR select_nodes
#' @importFrom DiagrammeR trav_in
#' @importFrom DiagrammeR trav_out
#' @importFrom DiagrammeR set_node_attrs_ws
#' @importFrom DiagrammeR clear_selection
#' @importFrom DiagrammeR create_node_df
#' @importFrom DiagrammeR create_edge_df
#' @importFrom DiagrammeR create_graph
#' @importFrom DiagrammeR set_global_graph_attrs
#' @importFrom DiagrammeR export_graph
#' @importFrom magrittr `%>%`
#' @importFrom digest digest
#' @importFrom BSgenome getBSgenome
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomeInfoDb seqnames
#' @importFrom AnnotationDbi saveDb
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures makeTxDbFromUCSC
#' @importFrom R.utils isGzipped
#' @importFrom R.utils gunzip
#' @importFrom R.utils isBzipped
#' @importFrom R.utils bunzip2
#' @importFrom GenomicRanges coverage
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics subset
#' @importFrom rmarkdown render
#' @importFrom knitr knit
#' @importFrom markdown markdownToHTML
#' @useDynLib ATACFlow
#' @rdname atacPeakQC
#' @export 
atacPeakQC<-function(atacProc, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"), bedInput = NULL){
    atacproc<-PeakQC$new(atacProc, reportOutput=reportOutput,qcbedInput = qcbedInput,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacPeakQC
#' @export
peakQC<-function(bedInput, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed")){
    atacproc<-PeakQC$new(atacProc=NULL, reportOutput=reportOutput,qcbedInput = qcbedInput,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    invisible(atacproc)
}
