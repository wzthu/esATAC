GraphMng <- R6Class(
  classname = "GraphMng",

  public = list(
    initialize = function(){
      edges1<-c(
       "SeqFile","Renamer",
        "SeqFile","FastQC",
        "SeqFile","RemoveAdapter",
        "Renamer","RemoveAdapter",
        "Renamer","FastQC",
       "Renamer", "QCreporter",
       "SeqFile", "QCreporter",
       "SamToBam", "Rsortbam",
       "SamToBam", "QCreporter",
       "RemoveAdapter","BowtieMapping",
       "BowtieMapping", "SamToBed",
       "BowtieMapping", "SamToBam",
       "RemoveAdapter","Bowtie2Mapping",
       "Bowtie2Mapping","SamToBed",
       "Bowtie2Mapping", "SamToBam"
      )
      private$graphDep1<-graph(edges = edges1)
      private$vtx1<-vertex.attributes(private$graphDep1)
    },
    checkRelation1 = function(resultProcName,procName){
      return(are.connected(private$graphDep1,resultProcName,procName))
    },
    printGraph = function(){
      #plot(private$graphDep1, layout=layout.reingold.tilford)
      plot(private$graphDep1)
    },
    getNextProcs1 = function(procName){
      v <- neighbors(graph = private$graphDep1,v = procName, mode = "out");
      nextProc<-private$vtx1[v];
      print("This object is a valid input for:")
      print(nextProc);
      return(nextProc);
    },
    getPrevProcs1 = function(procName){
      v <- neighbors(graph = private$graphDep1,v = procName, mode = "in");
      preProc<-private$vtx1[v];
      print("This object accepts the objects below for initialization:")
      print(preProc);
      return(preProc);
    },
    getProcList = function(){
      vtx<-as.character(vertex.attributes(private$graphDep1)$name)
      vsize<-length(vtx)
      idx<-list()
      for(i in 1:vsize){
        idx[[vtx[i]]]<-i
      }
      return(idx)
    },
    getSubGraphTopo = function(idx,inputAtacPorc,procList){
      sbgh<-induced_subgraph(private$graphDep1,as.numeric(idx))
      vtx<-as.character(vertex.attributes(sbgh)$name)
      vorder<-as.numeric(topo_sort(sbgh,"out"))
      nameList<-vtx[vorder];
      vsize<-length(nameList)
      atacProcList<-list()
      i<-1
      atacProc<-procList[[nameList[i]]]$new(inputAtacPorc)
      atacProc$process()
      atacProcList[[nameList[i]]]<-atacProc
      for(i in 2:vsize){
        v <- neighbors(graph = sbgp,v = nameList[i], mode = "in")
        if(length(v)>=1){
          atacProc<-procList[[nameList[i]]]$new(atacProcList[vtx[v[0]]])
        }else if(length(v)>=2){
          atacProc<-procList[[nameList[i]]]$new(atacProcList[vtx[v[0]]],atacProcList[vtx[v[1]]])
        }
        atacProc$process()
        atacProcList[[nameList[i]]]<-atacProc
      }
      return(atacProcList)
    }


  ),
  private = list(
    graphDep1=NULL,
    vtx1=NULL
    #graphDep2,
    #graphDep3,
  )

)

