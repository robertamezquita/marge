rmarkdown::render("README.Rmd", "md_document", "README.md")
devtools::document()
pkgdown::build_site()
devtools::build(path = ".")

