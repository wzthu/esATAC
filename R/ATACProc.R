
#' @importFrom pipeFrame Step
#' @importFrom pipeFrame addEdges
#' @importFrom pipeFrame checkAndInstallBSgenome
#' @importFrom pipeFrame checkAndInstallGenomeFa
#' @importFrom pipeFrame checkFileCreatable
#' @importFrom pipeFrame checkFileExist
#' @importFrom pipeFrame checkPathExist
#' @importFrom pipeFrame getBasenamePrefix
#' @importFrom pipeFrame getGenome
#' @importFrom pipeFrame getJobDir
#' @importFrom pipeFrame getJobName
#' @importFrom pipeFrame getNextSteps
#' @importFrom pipeFrame getPathPrefix
#' @importFrom pipeFrame getPrevSteps
#' @importFrom pipeFrame getRef
#' @importFrom pipeFrame getRefDir
#' @importFrom pipeFrame getRefFiles
#' @importFrom pipeFrame getRefRc
#' @importFrom pipeFrame getThreads
#' @importFrom pipeFrame getTmpDir
#' @importFrom pipeFrame getValidGenome
#' @importFrom pipeFrame initPipeFrame
#' @importFrom pipeFrame processing
#' @importFrom pipeFrame genReport
#' @importFrom pipeFrame runWithFinishCheck
#' @importFrom pipeFrame setGenome
#' @importFrom pipeFrame setJobName
#' @importFrom pipeFrame setRefDir
#' @importFrom pipeFrame setThreads
#' @importFrom pipeFrame setTmpDir
#' @importFrom pipeFrame addFileSuffix
#' @importFrom pipeFrame loadConfig
#' @importFrom pipeFrame saveConfig
#' @importFrom pipeFrame input "input<-" output "output<-" "param" "param<-" 
#' @importFrom pipeFrame "property" "property<-" "report" "report<-"
#' @importFrom pipeFrame regAttachedStep
#' @importFrom pipeFrame checkAndInstallTxDb
#' @importFrom pipeFrame checkAndInstallOrgDb
#' @importFrom pipeFrame ignoreCheck
#' @importFrom pipeFrame setPipeName
#' @importMethodsFrom pipeFrame checkAllPath
#' @importMethodsFrom pipeFrame checkRequireParam
#' @importMethodsFrom pipeFrame clearStepCache
#' @importMethodsFrom pipeFrame getAutoPath
#' @importMethodsFrom pipeFrame getParam
#' @importMethodsFrom pipeFrame getParamItems
#' @importMethodsFrom pipeFrame getParamMD5Path
#' @importMethodsFrom pipeFrame stepID
#' @importMethodsFrom pipeFrame stepName
#' @importMethodsFrom pipeFrame getStepWorkDir
#' @importMethodsFrom pipeFrame init
#' @importMethodsFrom pipeFrame isReady
#' @importMethodsFrom pipeFrame writeLog
#' 
#' 


#' @name ATACProc-class
#' @rdname  ATACProc-class
#' @title Base class of this package
#' @description This class is inherit from \code{Step} in pipeFrame package,
#' no more method is extended or override. Please see \code{Step} class for detail.
#' @export
setClass(Class = "ATACProc",
         contains = "Step"
)
