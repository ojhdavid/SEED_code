#!/usr/bin/env Rscript

commit <- "2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0"
repo <- paste0("YingMa0107/CARD@", commit)
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Install the environment from visium/environment-card.yml first.")
}

# Install MuSiC, required by CARD at runtime.
remotes::install_github("xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c",
                        upgrade = "never", dependencies = NA)
if (!requireNamespace("MuSiC", quietly = TRUE) || as.character(utils::packageVersion("MuSiC")) != "1.0.0") {
  stop("MuSiC 1.0.0 could not be loaded.")
}
remotes::install_github(repo, upgrade = "never", dependencies = NA)

if (!requireNamespace("CARD", quietly = TRUE)) {
  stop("CARD installation completed without a loadable CARD namespace.")
}

cat("CARD version:", as.character(utils::packageVersion("CARD")), "\n")
desc <- utils::packageDescription("CARD")
cat("RemoteSha:", ifelse(is.null(desc$RemoteSha), "not recorded", desc$RemoteSha), "\n")
