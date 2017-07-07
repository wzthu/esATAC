fseq <- function(argvs){
  argvs<-as.character(argvs)
  .jcall("edu/duke/igsp/gkde/Main","V","main",argvs)
}
