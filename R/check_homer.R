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
        tmp <- stringr::str_split(Sys.which(cmd), '/')
        h <- map(tmp, stringr::str_detect, 'homer') %>% unlist() %>% which()
        hp <- paste0(tmp[[1]][1:h] %>% paste0(collapse = '/'), '/')
        options('homer_path' = hp)

        return(invisible(0))
    }

    ## Not a default success, but user has set `options('homer_path')` (success or fail)
    if (output != 0 & !is.null(options('homer_path')$homer_path)) {
        homer_base <- options('homer_path')
        cmd2 <- paste0(homer_base, '/bin/findMotifsGenome.pl')
        output2 <- suppressWarnings(system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE))

        if (output2 == 0) {
            message("HOMER installation found via options('homer_path')")
            return(invisible(0))
        } else {
            message("HOMER installation not found\n")
            message("Check `options('homer_path')` and that HOMER is properly installed. Be sure to set the path to point to the base directory of your HOMER install.\n")
            return(invisible(127))
        }
    }

    ## Failure - not findable by default and not set `options('homer_path)`
    if (output != 0 & is.null(options('homer_path')$homer_path)) {
        message("HOMER installation not found - check $PATH and/or set `options('homer_path')`\n")
        message("\n")
        message(paste("If HOMER is not in $PATH (run `Sys.getenv()[['PATH']]`",
                      "try adding to PATH via `~/.bashrc` or equivalent ",
                      "(see HOMER website for details). If that ",
                      "does not work, provide the path by setting in R: \n",
                      "`options('homer_path' = '/your/full/path/to/homer/')` \n",
                      "You can add this line to your `~/.Rprofile` so it is ",
                      "loaded automagically every session.\n", sep = '\n'))
        message("\n")
        message("Note that conda installs are not well supported at this time.\n")
        message("\n")
        message(paste0("Lastly, try reinstalling HOMER from source.\n",
                       "See the HOMER install instructions at:\n",
                       "http://homer.ucsd.edu/homer/introduction/install.html"))
        return(invisible(127))
    }
}

#' Get the base installation directory of HOMER
#'
#' Grabs the base directory of HOMER installation from \code{options('homer_path')};
#' the correct base path should be set with package attachment, or alert user
#' to check/set the path manually. This function is kept to keep code from breaking.
#'
#' @param bin appends \code{bin/} to the homer base path (default: TRUE)
#' 
#' @return character vector with base directory of HOMER installation
#' @export
get_homer_base <- function(bin = TRUE) {
    homer_base <- options('homer_path')$homer_path
    if (bin == TRUE) {
        homer_base <- paste0(homer_base, 'bin/')
    }
    return(homer_base)
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
    homer_base <- get_homer_base(bin = FALSE)
    cmd <- "perl"
    args <- paste0(homer_base, "/configureHomer.pl ",
                   "-list")

    ## Parse output to remove header
    output <- system2(cmd, args, stdout = TRUE, stderr = TRUE) %>%
        stringr::str_replace_all("\t", " ") %>%
        .[stringr::str_detect(., "v[0-9]")]

    ## Parse each line further to generate tibble of packages
    .parse_out <- function(x) { x %>% t() %>% tibble::as_tibble() }
    parsed <- stringr::str_split(output, " ", n = 4) %>%
        purrr::map(.parse_out) %>%
        dplyr::bind_rows() %>%
        dplyr::rename_("status" = "V1",
                       "package" = "V2",
                       "version" = "V3",
                       "description" = "V4")

    ## Return either only installed or all packages
    if (installed_only == TRUE) {
        dplyr::filter_(parsed, "status == '+'") %>%
            return
    } else {
        return(parsed)
    }
}
