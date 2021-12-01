#' Check for working HOMER installation
#'
#' Used to check that HOMER is installed and working on package attach,
#' and reads/sets the \code{options('homer_path')} if user has not.
#'
#' @return logical; \code{TRUE} for success, \code{FALSE} otherwise
#' @export
check_homer <- function() {
    cmd <- "findMotifsGenome.pl"
    output <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))

    ## Success - found in PATH
    if (output == 0) {
        message("HOMER installation found.")

        ## Parse homer_base path and assign via options
      # remove -a parameter as it breaks on many systems
        loc <- system('type findMotifsGenome.pl', intern = TRUE)
        tmp <- stringr::str_split(loc, ' ')[[1]][3]
        hp <- dirname(dirname(tmp)) # expectation is cmd is nested two levels deep
        options('homer_path' = hp)
        return(invisible(0))
    }

    ## Not a default success, but user has set `options('homer_path')` (success or fail)
    if (output != 0 && !is.null(options('homer_path')$homer_path)) {
        homer_base <- options('homer_path')
        cmd2 <- paste0(homer_base, '/bin/findMotifsGenome.pl')
        output2 <- suppressWarnings(system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE))
        
        if (output2 == 0) {
            message("HOMER installation found via options('homer_path')")
            message("Amending R's $PATH via `Sys.setenv` to add all HOMER utils ..")

            homer_bin <- paste0(homer_base$homer_path, '/bin')
            path.env <- Sys.getenv('PATH')
            Sys.setenv('PATH' = paste0(homer_bin, ':', path.env))

            ## Recheck if amending worked
            output2 <- suppressWarnings(system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE))
            output3 <- suppressWarnings(system('findKnownMotifs.pl', ignore.stdout = TRUE, ignore.stderr = TRUE))

            if (output2 == 0 & output3 == 0) {
                message("Successfully added HOMER utils.")
                return(invisible(0))
            } else {
                message("Failed to load HOMER utils. Check your HOMER install, and failing that, file a Github issue.")
                return(invisible(127))
            }
        } else {
            message("HOMER installation not found.\n")
            message("Check `options('homer_path')` and that HOMER is properly installed. Be sure to set the path to point to the base directory of your HOMER install.\n")
            return(invisible(127))
        }
    }

    ## Failure - not findable by default and not set `options('homer_path)`
    if (output != 0 & is.null(options('homer_path')$homer_path)) {
        message("HOMER installation not found.\n")
        message(paste("Make sure HOMER is properly installed on your system.",
                      "See http://homer.ucsd.edu/homer/index.html",
                      "",
                      "Following that, the first and safest fix is to add the following",
                      "line in your ~/.Rprofile (or run it interactively upon starting R):",
                      "",
                      "options('homer_path' = '/path/to/homer')",
                      "",
                      "If `check_homer()` still doesn't work, then HOMER may not be in your $PATH.",
                      "(e.g. your system cannot find the HOMER utils)",
                      "Add HOMER to your PATH by adding the following line to your ~/.bashrc:",
                      "",
                      "PATH=$PATH:/your/path/to/homer/bin/",
                      "",
                      "Then `source ~/.bashrc` to add HOMER to your path.",
                      "",
                      "If still no good, file an issue on Github:",
                      "https://github.com/robertamezquita/marge",
                      sep = '\n'))
        return(invisible(127))
    }
}

#' Get the base installation directory of HOMER
#'
#' Grabs the base directory of HOMER installation from \code{options('homer_path')};
#' the correct base path should be set with package attachment, or alert user
#' to check/set the path manually. This function is kept to keep code from breaking.
#'
#' @return character vector with base directory of HOMER installation
#' @export
get_homer_bin <- function() {
    homer_base <- options('homer_path')$homer_path
    homer_bin <- paste0(homer_base, '/bin/')
    return(homer_bin)
}


#' List installed/available HOMER packages
#'
#' Runs the configureHomer.pl script to get list of installed and
#' available packages, parsing the output into a tidy tibble.
#'
#' @param installed_only only show installed packages (default: TRUE)
#' 
#' @return \code{tibble} with status of installed (+) or not (-),
#'   the package, version, and description.
#' @export
list_homer_packages <- function(installed_only = TRUE) {
    homer_base <- dirname(get_homer_bin())

    homer_conf <- paste0(homer_base, '/configureHomer.pl')

    if (!file.exists(homer_conf)) {
        stop(paste("The configureHomer.pl script was not found at the",
                   "base of the HOMER directory:",
                   "",
                   homer_base,
                   "",
                   "Place the configureHomer.pl script in the above directory.",
                   sep = '\n'))
    }

    cmd <- paste('perl', homer_conf, '-list')
    output <- system(cmd, intern = TRUE, ignore.stderr = TRUE) %>%
        stringr::str_replace_all("\t", " ") %>%
        .[stringr::str_detect(., "v[0-9]")]
    
    ## Parse each line further to generate tibble of packages
    .parse_out <- function(x) { x %>% t() %>% tibble::as_tibble() }
    parsed <- stringr::str_split(output, " ", n = 4) %>%
        purrr::map(.parse_out) %>%
        dplyr::bind_rows() %>%
        dplyr::rename(status = V1,
                       package = V2,
                       version = V3,
                       description = V4)

    ## Return either only installed or all packages
    if (installed_only) {
        dplyr::filter(parsed, status == '+') %>% return()
    } else {
        return(parsed)
    }
}
