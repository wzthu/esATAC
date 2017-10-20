#' @importFrom igraph neighbors 
setClass(Class = "GraphMng",
         slots = list(
             private = "list"
         ))
setMethod(f = "initialize",
          signature = "GraphMng", 
          definition = function(.Object,...){
              .Object@private <- list(
                  graph=NULL,
                  graph1=NULL,
                  graph2=NULL,
                  graphDep=NULL,
                  graphDep1=NULL,
                  graphDep2=NULL,
                  vtx1=NULL,
                  vtx2=NULL,
                  allProcName=NULL)
              edges1<-c(
                  "UnzipAndMerge", "Renamer",
                  "UnzipAndMerge", "FastQC",
                  "UnzipAndMerge", "RemoveAdapter",
                  "Renamer", "RemoveAdapter",
                  "Renamer", "FastQC",
                  "SamToBam", "Rsortbam",
                  "SamToBam", "FastQC",
                  "SamToBam", "BamToBed",
                  "Rsortbam", "BamToBed",
                  "RemoveAdapter", "Bowtie2Mapping",
                  "Bowtie2Mapping", "SamToBed",
                  "Bowtie2Mapping", "SamToBam",
                  "Bowtie2Mapping", "LibComplexQC",
                  "SamToBed", "PeakCallingFseq",
                  "PeakCallingFseq", "RPeakComp",
                  "SamToBed", "FregLenDistr",
                  "SamToBed", "TSSQC",
                  "SamToBed", "FRiPQC",
                  "SamToBed", "BedToBigWig",
                  "BamToBed", "BedUtils",
                  "BedUtils", "BedToBigWig",
                  "BedUtils", "BedUtils",
                  "BedUtils", "PeakCallingFseq",
                  "BedUtils", "FregLenDistr",
                  "BedUtils", "TSSQC",
                  "BedUtils", "FRiPQC",
                  "BedUtils", "CutSitePre",
                  "SamToBed", "BedUtils",
                  "PeakCallingFseq", "PeakQC",
                  "PeakCallingFseq", "RMotifScan",
                  "PeakCallingFseq", "RPeakAnno",
                  "PeakCallingFseq", "RSNPs",
                  "SamToBed", "CutSitePre",
                  "CutSitePre", "CutSiteCountR",
                  "RMotifScan", "CutSiteCountR",
                  "RPeakAnno", "RGo",
                  "RMotifScan", "RSNPs",
                  "RPeakComp", "RMotifScan",
                  "RPeakComp", "RPeakAnno",
                  "RPeakComp", "RSNPs"
                  
              )
              #edges1<-sapply(edges1,function(x) x$classname)
              
              
              nodes1<-unique(edges1)
              edges2<-c(
                  "PeakCallingFseq","FRiPQC",
                  "PeakCallingFseq","RPeakComp"
              )
              #edges2<-sapply(edges2,function(x) x$classname)
              nodes2<-unique(edges2)
              
              allNodes<-unique(edges1)
              .Object@private$allProcName<-allNodes
              #create first parameter graph
              gph1<-getGraph(.Object,edges1)
              
              #create second parameter graph
              gph2<-getGraph(.Object,edges2)
              
              
              #create merged graph
              edges<-c(edges1,edges2)
              # for(i in seq(1,length(edges2),2)){
              #     for(j in seq(1,length(edges),2)){
              #         if(!(edges[j]==edges1[i]&&edges[j+1]==edges1[i+1])){
              #             edges<-c(edges,edges2[i],edges2[i+1])
              #         }
              #     }
              # }
              
              gph<-getGraph(.Object,edges)
              
              
              
              
              .Object@private$graph<-gph
              .Object@private$graph1<-gph1
              .Object@private$graph2<-gph2
              .Object@private$graphDep<-to_igraph(gph)
              vertex.attributes(.Object@private$graphDep)$name<-vertex.attributes(.Object@private$graphDep)$label
              .Object@private$graphDep1<-to_igraph(gph1)
              vertex.attributes(.Object@private$graphDep1)$name<-vertex.attributes(.Object@private$graphDep1)$label
              .Object@private$graphDep2<-to_igraph(gph2)
              vertex.attributes(.Object@private$graphDep2)$name<-vertex.attributes(.Object@private$graphDep2)$label
              .Object@private$vtx1<-vertex.attributes(.Object@private$graphDep1)
              .Object@private$vtx2<-vertex.attributes(.Object@private$graphDep2)
              .Object
          })

setGeneric(name = "checkRelation1",
           def = function(graphMngObj,resultProcName,procName,...){
               standardGeneric("checkRelation1")
           })
setMethod(f = "checkRelation1",
          signature = "GraphMng",
          definition = function(graphMngObj,resultProcName,procName,...){
              return(are.connected(graphMngObj@private$graphDep1,resultProcName,procName))
          })


setGeneric(name = "checkRelation2",
           def = function(graphMngObj,...){
               standardGeneric("checkRelation2")
           })
setMethod(f = "checkRelation2",
          signature = "GraphMng",
          definition = function(graphMngObj,...){
              return(are.connected(graphMngObj@private$graphDep2,resultProcName,procName))
          })

setGeneric(name = "printMap",
           def = function(graphMngObj,...,procName=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
               standardGeneric("printMap")
           })
setMethod(f = "printMap",
          signature = "GraphMng",
          definition = function(graphMngObj,...,procName=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
              if(is.null(procName)){
                  if(display){
                      graphMngObj@private$graph%>%render_graph()
                  }else{
                      graphMngObj@private$graph%>%export_graph(file_name = file.path(.obtainConfigure("tmpdir"),"currentMap.pdf"),file_type="pdf")
                  }
                  
              }else{
                  tempMap<-graphMngObj@private$graph
                  if(preProc){
                      tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                          trav_in() %>%set_node_attrs_ws("fillcolor", "Gold") %>%clear_selection()
                  }
                  if(nextProc){
                      tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                          trav_out() %>%set_node_attrs_ws("fillcolor", "SpringGreen") %>%clear_selection()
                  }
                  if(curProc){
                      tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                          set_node_attrs_ws("fillcolor", "Red") %>%clear_selection()
                  }
                  if(display){
                      render_graph(tempMap)
                  }else{
                      tempMap%>%export_graph(file_name = file.path(.obtainConfigure("tmpdir"),"currentMap.pdf"),file_type="pdf")
                      
                  }
                  
              }
          })

setGeneric(name = "getNextProcs",
           def = function(graphMngObj,procName,...){
               standardGeneric("getNextProcs")
           })
setMethod(f = "getNextProcs",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              getNextProcs1(graphMngObj,procName)
              getNextProcs2(graphMngObj,procName)
          })

setGeneric(name = "getNextProcs1",
           def = function(graphMngObj,procName,...){
               standardGeneric("getNextProcs1")
           })
setMethod(f = "getNextProcs1",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              v <- neighbors(graph = graphMngObj@private$graphDep1,v = procName, mode = "out");
              nextProc<-graphMngObj@private$vtx1[v];
              print("This object is a valid first parameter for:")
              print(nextProc);
              return(nextProc);
          })

setGeneric(name = "getNextProcs2",
           def = function(graphMngObj,procName,...){
               standardGeneric("getNextProcs2")
           })
setMethod(f = "getNextProcs2",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              v <- neighbors(graph = graphMngObj@private$graphDep2,v = procName, mode = "out");
              nextProc<-graphMngObj@private$vtx2[v];
              print("This object is a valid second parameter for:")
              print(nextProc);
              return(nextProc);
          })


setGeneric(name = "getPrevProcs",
           def = function(graphMngObj,procName,...){
               standardGeneric("getPrevProcs")
           })
setMethod(f = "getPrevProcs",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              getPrevProcs1(graphMngObj,procName)
              getPrevProcs2(graphMngObj,procName)
          })

setGeneric(name = "getPrevProcs1",
           def = function(graphMngObj,procName,...){
               standardGeneric("getPrevProcs1")
           })
setMethod(f = "getPrevProcs1",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              v <- neighbors(graph = graphMngObj@private$graphDep1,v = procName, mode = "in");
              preProc<-graphMngObj@private$vtx1[v];
              print("This object accepts the objects below as first parameter for initialization:")
              print(preProc);
              return(preProc);
          })

setGeneric(name = "getPrevProcs2",
           def = function(graphMngObj,procName,...){
               standardGeneric("getPrevProcs2")
           })
setMethod(f = "getPrevProcs2",
          signature = "GraphMng",
          definition = function(graphMngObj,procName,...){
              v <- neighbors(graph = graphMngObj@private$graphDep2,v = procName, mode = "in");
              preProc<-graphMngObj@private$vtx2[v];
              print("This object accepts the objects below as second parameter for initialization:")
              print(preProc);
              return(preProc);
          })

setGeneric(name = "getProcList",
           def = function(graphMngObj,...){
               standardGeneric("getProcList")
           })
setMethod(f = "getProcList",
          signature = "GraphMng",
          definition = function(graphMngObj,...){
              vtx<-as.character(vertex.attributes(graphMngObj@private$graphDep1)$name)
              vsize<-length(vtx)
              idx<-list()
              for(i in 1:vsize){
                  idx[[vtx[i]]]<-i
              }
              return(idx)
          })

setGeneric(name = "getGraph",
           def = function(graphMngObj,edges,...){
               standardGeneric("getGraph")
           })
setMethod(f = "getGraph",
          signature = "GraphMng",
          definition = function(graphMngObj,edges,...){
              if(is.null(graphMngObj@private$allProcName)){
                  stop("all ProcName can not be NULL")
              }
              nodes <- unique(edges)
              idx=as.vector(sapply(nodes,function(x) which(x==graphMngObj@private$allProcName)))
              
              ndf<- create_node_df(
                  idx = idx,
                  n = length(nodes),
                  label = nodes,
                  #shape = "circle",
                  #fixedsize = TRUE,
                  fontname = "Helvetica",
                  fillcolor = "DeepSkyBlue",
                  color = "White",
                  fontcolor = "Black",
                  style = "filled"
                  
                  
              )
              
              startidx<-as.vector(sapply(edges[seq(1,length(edges),2)],function(x) which(x==nodes)))
              endidx<-as.vector(sapply(edges[seq(2,length(edges),2)],function(x) which(x==nodes)))
              edf<-create_edge_df(
                  from = startidx,
                  to = endidx,
                  color = "DeepSkyBlue"
              )
              
              gph<-create_graph(nodes_df = ndf,edges_df = edf)%>%
                  set_global_graph_attrs(attr_type = "graph",attr = "layout",value = "dot")
              #set_global_graph_attrs(attr_type = "node",attr = "fontname",value = "Helvetica")%>%
              #set_global_graph_attrs(attr_type = "node",attr = "shape",value = "circle")%>%
              #set_global_graph_attrs(attr_type = "node",attr = "fixedsize",value = "true")
              return(gph)
          })

GraphMng <- R6Class(
  classname = "GraphMng",

  public = list(
    initialize = function(){
        edges1<-c(
            "UnzipAndMerge", "Renamer",
            "UnzipAndMerge", "FastQC",
            "UnzipAndMerge", "RemoveAdapter",
            "Renamer", "RemoveAdapter",
            "Renamer", "FastQC",
            "SamToBam", "Rsortbam",
            "SamToBam", "FastQC",
            "SamToBam", "BamToBed",
            "Rsortbam", "BamToBed",
            "RemoveAdapter", "Bowtie2Mapping",
            "Bowtie2Mapping", "SamToBed",
            "Bowtie2Mapping", "SamToBam",
            "Bowtie2Mapping", "LibComplexQC",
            "SamToBed", "PeakCallingFseq",
            "PeakCallingFseq", "RPeakComp",
            "SamToBed", "FregLenDistr",
            "SamToBed", "TSSQC",
            "SamToBed", "FRiPQC",
            "SamToBed", "BedToBigWig",
            "BamToBed", "BedUtils",
            "BedUtils", "BedToBigWig",
            "BedUtils", "BedUtils",
            "BedUtils", "PeakCallingFseq",
            "BedUtils", "FregLenDistr",
            "BedUtils", "TSSQC",
            "BedUtils", "FRiPQC",
            "BedUtils", "CutSitePre",
            "SamToBed", "BedUtils",
            "PeakCallingFseq", "PeakQC",
            "PeakCallingFseq", "RMotifScan",
            "PeakCallingFseq", "RPeakAnno",
            "PeakCallingFseq", "RSNPs",
            "SamToBed", "CutSitePre",
            "CutSitePre", "CutSiteCountR",
            "RMotifScan", "CutSiteCountR",
            "RPeakAnno", "RGo",
            "RMotifScan", "RSNPs",
            "RPeakComp", "RMotifScan",
            "RPeakComp", "RPeakAnno",
            "RPeakComp", "RSNPs"

        )
        #edges1<-sapply(edges1,function(x) x$classname)


        nodes1<-unique(edges1)
        edges2<-c(
            "PeakCallingFseq","FRiPQC",
            "PeakCallingFseq","RPeakComp"
        )
        #edges2<-sapply(edges2,function(x) x$classname)
        nodes2<-unique(edges2)

        allNodes<-unique(edges1)
        private$allProcName<-allNodes
        #create first parameter graph
        gph1<-private$getGraph(edges1)

        #create second parameter graph
        gph2<-private$getGraph(edges2)


        #create merged graph
        edges<-c(edges1,edges2)
        # for(i in seq(1,length(edges2),2)){
        #     for(j in seq(1,length(edges),2)){
        #         if(!(edges[j]==edges1[i]&&edges[j+1]==edges1[i+1])){
        #             edges<-c(edges,edges2[i],edges2[i+1])
        #         }
        #     }
        # }

        gph<-private$getGraph(edges)




        private$graph<-gph
        private$graph1<-gph1
        private$graph2<-gph2
        private$graphDep<-to_igraph(gph)
        vertex.attributes(private$graphDep)$name<-vertex.attributes(private$graphDep)$label
        private$graphDep1<-to_igraph(gph1)
        vertex.attributes(private$graphDep1)$name<-vertex.attributes(private$graphDep1)$label
        private$graphDep2<-to_igraph(gph2)
        vertex.attributes(private$graphDep2)$name<-vertex.attributes(private$graphDep2)$label
        private$vtx1<-vertex.attributes(private$graphDep1)
        private$vtx2<-vertex.attributes(private$graphDep2)
    },

    checkRelation1 = function(resultProcName,procName){
      return(are.connected(private$graphDep1,resultProcName,procName))
    },
    checkRelation2 = function(resultProcName,procName){
        return(are.connected(private$graphDep2,resultProcName,procName))
    },

    printMap = function(procName=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
      #plot(private$graphDep1, layout=layout.reingold.tilford)

        if(is.null(procName)){
            if(display){
                private$graph%>%render_graph()
            }else{
                private$graph%>%export_graph(file_name = file.path(.obtainConfigure("tmpdir"),"currentMap.pdf"),file_type="pdf")
            }

        }else{
            tempMap<-private$graph
            if(preProc){
                tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                    trav_in() %>%set_node_attrs_ws("fillcolor", "Gold") %>%clear_selection()
            }
            if(nextProc){
                tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                    trav_out() %>%set_node_attrs_ws("fillcolor", "SpringGreen") %>%clear_selection()
            }
            if(curProc){
                tempMap<-tempMap%>%select_nodes(conditions = label==procName)%>%
                    set_node_attrs_ws("fillcolor", "Red") %>%clear_selection()
            }
            if(display){
                render_graph(tempMap)
            }else{
                tempMap%>%export_graph(file_name = file.path(.obtainConfigure("tmpdir"),"currentMap.pdf"),file_type="pdf")

            }

        }


    },
    getNextProcs = function(procName){
        self$getNextProcs1(procName)
        self$getNextProcs2(procName)
    },
    getNextProcs1 = function(procName){
      v <- neighbors(graph = private$graphDep1,v = procName, mode = "out");
      nextProc<-private$vtx1[v];
      print("This object is a valid first parameter for:")
      print(nextProc);
      return(nextProc);
    },
    getNextProcs2 = function(procName){
        v <- neighbors(graph = private$graphDep2,v = procName, mode = "out");
        nextProc<-private$vtx2[v];
        print("This object is a valid second parameter for:")
        print(nextProc);
        return(nextProc);
    },
    getPrevProcs = function(procName){
        self$getPrevProcs1(procName)
        self$getPrevProcs2(procName)
    },
    getPrevProcs1 = function(procName){
      v <- neighbors(graph = private$graphDep1,v = procName, mode = "in");
      preProc<-private$vtx1[v];
      print("This object accepts the objects below as first parameter for initialization:")
      print(preProc);
      return(preProc);
    },
    getPrevProcs2 = function(procName){
        v <- neighbors(graph = private$graphDep2,v = procName, mode = "in");
        preProc<-private$vtx2[v];
        print("This object accepts the objects below as second parameter for initialization:")
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
    graph=NULL,
    graph1=NULL,
    graph2=NULL,
    graphDep=NULL,
    graphDep1=NULL,
    graphDep2=NULL,
    vtx1=NULL,
    vtx2=NULL,
    allProcName=NULL,
    getGraph = function(edges){
        if(is.null(private$allProcName)){
            stop("all ProcName can not be NULL")
        }
        nodes <- unique(edges)
        idx=as.vector(sapply(nodes,function(x) which(x==private$allProcName)))

        ndf<- create_node_df(
            idx = idx,
            n = length(nodes),
            label = nodes,
            #shape = "circle",
            #fixedsize = TRUE,
            fontname = "Helvetica",
            fillcolor = "DeepSkyBlue",
            color = "White",
            fontcolor = "Black",
            style = "filled"


        )

        startidx<-as.vector(sapply(edges[seq(1,length(edges),2)],function(x) which(x==nodes)))
        endidx<-as.vector(sapply(edges[seq(2,length(edges),2)],function(x) which(x==nodes)))
        edf<-create_edge_df(
            from = startidx,
            to = endidx,
            color = "DeepSkyBlue"
        )

        gph<-create_graph(nodes_df = ndf,edges_df = edf)%>%
            set_global_graph_attrs(attr_type = "graph",attr = "layout",value = "dot")
            #set_global_graph_attrs(attr_type = "node",attr = "fontname",value = "Helvetica")%>%
            #set_global_graph_attrs(attr_type = "node",attr = "shape",value = "circle")%>%
            #set_global_graph_attrs(attr_type = "node",attr = "fixedsize",value = "true")
        return(gph)
    }

  )

)

graphMng<-new("GraphMng")


