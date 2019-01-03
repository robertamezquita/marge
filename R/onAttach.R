.onAttach <- function(libname, pkgname) {
    ## Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        ## Grab version
        v <- utils::packageVersion("marge")

        ## Check if its a dev version (9999)
        if (!is.na(v[1, 4])) {
            if (v[1, 4] == 9999) {
                dev <- 'dev'
            } else {
                dev <- 'bugfix'
            }
        } else {
            dev <- 'stable'
        }

        ## Startup messages and check for HOMER install
        packageStartupMessage("marge ", v, paste0(" (", dev, " version)"))
        packageStartupMessage("Checking for HOMER installation...")
        value <- check_homer() ## does all the error checking and messaging
        ## packageStartupMessage("Read the docs at https://marge.aerobatic.io")
    }
}

