library(rmarkdown)

#' Render report
#'
#' @param script full path to R script to render
#' @return saves a pdf/html
#' @export
#' main()

main <- function(script){
    render(script)
}

############
## Inputs
############

args <- commandArgs(trailingOnly = TRUE)
print(args)

main(args[1])
