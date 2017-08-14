
atacPipe <- function(fastqInput1,fastqInput2=NULL,tmpdir=NULL,saveTmpFiles=FALSE,renamer=TRUE,removeAdapter=TRUE){
  atacInputFile(fastqInput1 = fastqInput1, fastqInput2 = fastqInput2)
  gph<-GraphMng$new()
  vlist<-gph$getProcList();
  if(!renamer){
    vlist$renamer<-NULL
  }
  if(!removeAdapter){
    vlist$removeAdapter<-NULL
  }
  return(gph$getSubGraphTopo(as.numeric(vlist),atacInputFile,list(Renamer=Renamer,RemoveAdapter=RemoveAdapter)))
}

atacPrintMap <-function(atacProc=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE){
    if(is.null(atacProc)){
        GraphMng$new()$printMap()
    }else if(class(atacProc)=="character"){
        GraphMng$new()$printMap(atacProc,preProc,nextProc,curProc)
    }else{
        atacProc$printMap(preProc,nextProc,curProc)
    }
}

atacPrintNextList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getNextProcs(atacProc)
    }else{
        GraphMng$new()$getNextProcs(atacProc$getProcName())
    }
}

atacPrintPrevList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getPrevProcs(atacProc)
    }else{
        GraphMng$new()$getPrevProcs(atacProc$getProcName())
    }
}














#' Mapping reads to the reference using Rbowtie, if output file do not be specified, the output will be named mapping_result.sam
#' @param seq_file A full path of the fa file(containing fa file). For single end, using a list; for paired end, using a list(length = 2).
#' @param ref_file Character scalar. The path to the bowtie index and prefix to align against, in the form </path/to/index>/<prefix>.
#' @param out_file path and output file name, 'E:\\RATAC_test\\output\\output.sam'
#' @param seq_type sequence type, "single", "paired"
#' @export
atacMappingBt <- function(atacProc = NULL, fileInput = NULL, Reference = NULL, fileOutput = NULL, In_type = NULL){
  tmp <- BowtieMapping$new(atacProc, fileInput, Reference, fileOutput, In_type)
  tmp$processing()
  return(tmp)
}




#' Quality control using Quasr::qQCreport
#' @param input_file a c()
#' @export
atacQCReport <- function(atacProc = NULL, input_file = NULL, output_file = NULL){
  tmp <- QCreporter$new(atacProc, input_file, output_file)
  tmp$processing()
  return(tmp)
}



