ChrDivi <- R6::R6Class(
  classname = "ChrDivi",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, ReadsIfile = NULL, ReadsOpath = NULL, name = NULL,
                          sort = NULL, sort_col = NULL, editable = FALSE){
      super$initialize("ChrDivi", editable, list(arg1 = atacProc))

      # necessary parameters
      private$paramlist[["ReadsIfile"]] <- ReadsIfile
      private$paramlist[["sort"]] <- sort
      private$paramlist[["sort_col"]] <- sort_col
      # unnecessary parameters
      if(is.null(ReadsOpath)){
        private$paramlist[["ReadsOpath"]] <- paste0(dirname(private$paramlist[["ReadsIfile"]]),
                                                    "/", collapse = "")
      }else{
        private$paramlist[["ReadsOpath"]] <- ReadsOpath
        # add "/" in the output path
        if(substr(private$paramlist[["ReadsOpath"]], nchar(private$paramlist[["ReadsOpath"]]),
                  nchar(private$paramlist[["ReadsOpath"]])) != "/"){
          private$paramlist[["ReadsOpath"]] <- paste0(private$paramlist[["ReadsOpath"]],
                                                      "/", collapse = "")
        }
      }
      if(is.null(name)){
        private$paramlist[["name"]] <- "Output"
      }else{
        private$paramlist[["name"]] <- name
      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), #public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["ReadsIfile"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["ReadsOpath"]]))
      .chr_separate_call(ReadsIfile = private$paramlist[["ReadsIfile"]], ReadsOpath = private$paramlist[["ReadsOpath"]],
                         Name = private$paramlist[["name"]])
      if(!is.null(private$paramlist[["sort"]]) && (private$paramlist[["sort"]])){
        chr <- as.list(c(1:22, "X", "Y"))
        for(i in seq(1:24)){
          sfile <- paste(private$paramlist[["ReadsOpath"]], private$paramlist[["name"]],
                         "_chr", chr[[i]], ".bed", sep = "", collapse = "")
          data <- read.table(sfile)
          data <- data[order(data[, private$paramlist[["sort_col"]]]), ]
          write.table(data, file = sfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      }
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["ReadsIfile"]])){
        stop("Parameter ReadsIfile is required!")
      }
      if(!is.null(private$paramlist[["sort"]]) && (private$paramlist[["sort"]])){
        if(is.null(private$paramlist[["sort_col"]])){
          stop("Please specify which column you want sort!")
        }
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["ReadsIfile"]])
      private$checkPathExist(private$paramlist[["ReadsOpath"]])
    } # checkAllPath end

  ) # private end

) # class end


#' separate genome information by chromatin name.
#' @param atacProc Not using now, we will use it in the future.
#' @param Ifile Input bed file path, the first column is chromatin name, must be sorted by chromatin.
#' @param Opath The output path, an empty folder would be great.
#' If not specified, the output file will be writen in ReadsIfile's dirctionary.
#' @param prefix the prefix of the output name, format:prefix_chr*.bed, default:output.
#' @param sort TRUE or FALSE, sort every output file by a column.
#' @param sort_col Which column you want to sort, if sort = TRUE, this parameter must be specified.
#' @export
atacChrDivi <- function(atacProc = NULL, Ifile = NULL, Opath = NULL,
                        prefix = NULL, sort = NULL, sort_col = null){
  tmp <- ChrDivi$new(atacProc, Ifile, Opath, prefix, sort, sort_col)
  tmp$process()
  return(tmp)
}
