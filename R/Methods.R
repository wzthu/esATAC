
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
        return(atacProc)
    }
}

atacPrintNextList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getNextProcs(atacProc)
    }else{
        GraphMng$new()$getNextProcs(atacProc$getProcName())
        return(atacProc)
    }
}

atacPrintPrevList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getPrevProcs(atacProc)
    }else{
        GraphMng$new()$getPrevProcs(atacProc$getProcName())
        return(atacProc)
    }
}


#' Quality control using Quasr::qQCreport
#' @param input_file a c()
#' @export
atacQCReport <- function(atacProc = NULL, input_file = NULL, output_file = NULL){
  tmp <- QCreporter$new(atacProc, input_file, output_file)
  tmp$processing()
  return(tmp)
}



