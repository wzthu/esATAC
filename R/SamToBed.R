SamToBed <- R6::R6Class(
  classname = "SamToBed",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= NULL,
                          reportPrefix =NULL,samInput = NULL, bedOutput = NULL, sortBed = TRUE, uniqueBed = TRUE,
                          minFregLen = 0,maxFregLen = 100, saveExtLen = FALSE, editable=FALSE){
      super$initialize("SamToBed",editable,list(arg1=atacProc))
      merge=merge[1]
      if(!is.null(atacProc)){
          private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
          regexProcName<-sprintf("(SAM|Sam|sam|%s)",atacProc$getProcName())
          if(merge=="auto"){
              if(private$singleEnd){
                  merge=FALSE
              }else{
                  merge=TRUE
              }

          }else if(merge=="yes"){
              merge=TRUE
          }else if(merge=="no"){
              merge=FALSE
          }else{
              stop(paste0("Invalid value of merge: ",merge))
          }

      }else{
          regexProcName<-"(SAM|Sam|sam)"
          if(!editable){
              if(merge=="auto"){
                  if(is.null(samInput)||!file.exists(samInput)){
                      stop(paste0("samInput file does not exist! ",samInput))
                  }
                  asamfile <- file(samInput, "r")
                  lines<-readLines(asamfile,n=1000)
                  close(asamfile)
                  lines<-lines[!grepl("^@",lines)]
                  merge=FALSE
                  for(i in 2:length(lines)){
                      code1=strsplit(lines[i-1],"\t")[[1]][1]
                      code2=strsplit(lines[i],"\t")[[1]][1]
                      if(code1==code2){
                          merge=TRUE
                      }
                  }
              }else if(merge=="yes"){
                    merge=TRUE
              }else if(merge=="no"){
                    merge=FALSE
              }else{
                  stop(paste0("Invalid value of merge: ",merge))
              }
          }else{
              merge=TRUE
          }


      }

      if(!is.null(samInput)){
          private$paramlist[["samInput"]] <- samInput;
      }
      if(is.null(bedOutput)){
          if(!is.null(private$paramlist[["samInput"]])){
              prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
              private$paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".bed"))
          }
      }else{
          private$paramlist[["bedOutput"]] <- bedOutput;
      }
      if(is.null(reportPrefix)){
          if(!is.null(private$paramlist[["samInput"]])){
              prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
              private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
          }
      }else{
          private$paramlist[["reportPrefix"]] <- reportPrefix;
      }


      private$paramlist[["merge"]] <- merge;
      private$paramlist[["posOffset"]] <- posOffset;
      private$paramlist[["negOffset"]] <- negOffset;
      private$paramlist[["filterList"]] <- chrFilterList;
      private$paramlist[["sortBed"]] <- sortBed
      private$paramlist[["uniqueBed"]] <- uniqueBed
      private$paramlist[["minFregLen"]] <- minFregLen
      private$paramlist[["maxFregLen"]] <- maxFregLen
      private$paramlist[["saveExtLen"]] <- saveExtLen


      private$paramValidation()
    }
  ),

    private = list(
        processing = function(){
            sink(private$paramlist[["reportPrefix"]])####------------------------------------------------------------std::cout>>Rcpp::Rout
            if(private$paramlist[["merge"]]){
                .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = private$paramlist[["bedOutput"]],
                                                                       posOffset = private$paramlist[["posOffset"]], negOffset = private$paramlist[["negOffset"]],
                                                                       sortBed = private$paramlist[["sortBed"]],uniqueBed = private$paramlist[["uniqueBed"]],
                                                                       filterList = private$paramlist[["filterList"]],minFregLen = private$paramlist[["minFregLen"]],
                                                                       maxFregLen = private$paramlist[["maxFregLen"]],saveExtLen = private$paramlist[["saveExtLen"]] )
            }else{
                .sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = private$paramlist[["bedOutput"]],
                                                                 posOffset = private$paramlist[["posOffset"]], negOffset = private$paramlist[["negOffset"]],
                                                                 sortBed = private$paramlist[["sortBed"]],uniqueBed = private$paramlist[["uniqueBed"]],
                                                                 filterList = private$paramlist[["filterList"]])
            }
            sink()

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["samInput"]])){
                stop(paste("samInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkFileCreatable(private$paramlist[["bedOutput"]]);
        }
  )


)



#' convert sam to bed
#' @param ATAC_obj obj returned from ATAC_mapping
#' @param samfile sam file dir
#' @param bedfile bed file dir
#' @param readlen reads length
#' @export
atacSam2Bed <- function(atacProc, reportPrefix =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                        samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                        saveExtLen = FALSE,uniqueBed = TRUE){
    atacproc <- SamToBed$new(atacProc=atacProc, merge=merge, posOffset=posOffset, negOffset=negOffset, chrFilterList=chrFilterList,
                             samInput=samInput, bedOutput=bedOutput, sortBed=sortBed, uniqueBed=uniqueBed, minFregLen=minFregLen,
                             maxFregLen=maxFregLen, saveExtLen=saveExtLen,reportPrefix=reportPrefix)
    atacproc$process()
    return(atacproc)
}
