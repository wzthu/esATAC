RSNPs <- R6::R6Class(
  classname = "RSNPs",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, snp.regions.file = NULL, bio.features.loc = NULL,
                          built.in.biofeatures = NULL, par.threads = NULL,
                          verbose = NULL, method.p = NULL, search.window = NULL,
                          output = NULL, editable = FALSE){
      super$initialize("RSNPs",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["snp.regions.file"]] <- snp.regions.file
      private$paramlist[["bio.features.loc"]] <- bio.features.loc
      private$paramlist[["built.in.biofeatures"]] <- built.in.biofeatures
      private$paramlist[["par.threads"]] <- par.threads
      private$paramlist[["verbose"]] <- verbose
      private$paramlist[["method.p"]] <- method.p
      private$paramlist[["search.window"]] <- search.window
      private$paramlist[["output"]] <- output
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # piblic end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("SNP source:%s", private$paramlist[["snp.regions.file"]]))
      private$writeLog(sprintf("bio feature source:%s", private$paramlist[["bio.features.loc"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["output"]]))
      tmp <- FunciSNP::getFSNPs(snp.regions.file = private$paramlist[["snp.regions.file"]],
                                bio.features.loc = private$paramlist[["bio.features.loc"]],
                                built.in.biofeatures = private$paramlist[["built.in.biofeatures"]],
                                par.threads = private$paramlist[["par.threads"]],
                                verbose = private$paramlist[["verbose"]],
                                method.p = private$paramlist[["method.p"]],
                                search.window = private$paramlist[["search.window"]])
      saveRDS(tmp, private$paramlist[["output"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["snp.regions.file"]])){
        stop("Parameter snp.regions.file is required!")
      }
      if(is.null(private$paramlist[["bio.features.loc"]])){
        stop("Parameter bio.features.loc is required!")
      }
      if(is.null(private$paramlist[["output"]])){
        stop("Parameter output is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["snp.regions.file"]])
      private$checkFileExist(private$paramlist[["bio.features.loc"]])
      private$checkPathExist(private$paramlist[["ReadsOpath"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Using FunciSNP to do SNP analysis.
#'
#' @param snp.regions.file path: Location of the regions file.
#' @param bio.features.loc path: Location of the biological features folder.
#'
SNPana <- function(atacProc = NULL, snp.regions.file = NULL, bio.features.loc = NULL,
                   built.in.biofeatures = TRUE,
                   par.threads = parallel::detectCores()/2,
                   verbose = par.threads < 2, method.p = "BH",
                   search.window = 200000, output = NULL){
  tmp <- RSNPs$new(atacProc, snp.regions.file, bio.features.loc,
                   built.in.biofeatures, par.threads,
                   verbose, method.p, search.window, output)
  tmp$process()
  return(tmp)
}
