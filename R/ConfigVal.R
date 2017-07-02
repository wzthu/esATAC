.ConfigClass<-R6Class(classname = ".ConfigClass",
  public = list(
    getAllConfigure = function(){
      print(private$configList)
    },
    getConfigure=function(item){
      return(private$configList[[item]]);
    },
    setConfigure=function(item,val){
      private$configList[[item]]<-val;
    }
    
  ),
  private = list(
    configList=list(threads=1,tmpdir=NULL)
  )
)

.configObj<-.ConfigClass$new()

getAllConfigure<-function(){
  .configObj$getAllConfigure();
}

getConfigure <- function(item){
  .configObj$getConfigure(item);
}

setConfigure<- function(item,val){
  .configObj$setConfigure(item,val);
}



