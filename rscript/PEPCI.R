## ---------------------------------------------------------------------------------------------------------------------
## Install dependencies

install.dependencies <- function(pkgs, repos) {
	pkgs <- pkgs[!(pkgs %in% rownames(installed.packages()))]
	if (length(pkgs) != 0) {
		if (identical(repos, "bioconductor")) {
			biocLite(pkgs, ask = FALSE)
		} else {
			install.packages(pkgs, repos = repos)
		}
	}
}
source("http://bioconductor.org/biocLite.R")

DEPENDENCIES <- list(
	"http://cran.rstudio.com" = c("RPMM", "ggplot2", "grid", "nlme", "quadprog"),
	"bioconductor" = c("illuminaio", "knitr"))
repos <- mapply(install.dependencies, pkgs = DEPENDENCIES, repos = as.list(names(DEPENDENCIES)))
rm(install.dependencies, DEPENDENCIES, repos)

## ---------------------------------------------------------------------------------------------------------------------
## Install PEPCI

fname <- c("Windows" = "zip", "Darwin" = "tgz", "Linux" = "tar.gz")
if (Sys.info()["sysname"] %in% names(fname)) {
	fname <- unname(fname[Sys.info()["sysname"]])
	package.type <- ifelse(fname == "tar.gz", "source", "binary")
	fname <- paste0("PEPCI_0.1.0.", fname)
	fname.server <- paste0("http://computational-epigenomics.com/downloads/PEPCI/", fname)
	fname.local <- file.path(tempdir(), fname)
	if (file.exists(fname.local)) {
		if (!file.remove(fname.local)) {
			stop(paste("Could not remove existing file", fname.local))
		}
	}
	if (download.file(fname.server, fname.local, quiet = TRUE) != 0) {
		stop(paste("Could not download file", fname.server))
	}
	install.packages(fname.local, repos = NULL, type = package.type)
	rm(fname, package.type, fname.server, fname.local)
} else {
	stop("Unsupported operating system")
}
