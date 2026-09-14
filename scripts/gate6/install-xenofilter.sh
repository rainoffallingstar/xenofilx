#!/usr/bin/env bash
# Install the pinned XenofilteR release into an `enva` analysis environment.
#
# XenofilteR imports GenomicAlignments, which loads GenomeInfoDb, which requires the
# GenomeInfoDbData annotation package at load time. The conda spec
# `bioconductor-genomeinfodbdata` resolves but does not place an R package named
# `GenomeInfoDbData` into the environment R library (linux-64, r-base 4.3, verified
# 2026-09-14), so `R CMD INSTALL` fails with:
#
#   package or namespace load failed for 'GenomeInfoDb' ...
#   there is no package called 'GenomeInfoDbData'
#
# The data package is therefore installed from the matching Bioconductor release and its
# version is asserted, so a silent Bioconductor drift fails the job instead of changing the
# analysis environment unnoticed.
#
# usage: install-xenofilter.sh <xenofilter-source-directory> <environment-name>
set -euo pipefail

xenofilter_source_directory="${1:?usage: install-xenofilter.sh <xenofilter-source-directory> <environment-name>}"
environment_name="${2:?usage: install-xenofilter.sh <xenofilter-source-directory> <environment-name>}"
enva_binary="${ENVA_BINARY:-.ci/enva/target/release/enva}"
bioconductor_release="${BIOCONDUCTOR_RELEASE:-3.18}"
expected_genomeinfodbdata_version="${GENOMEINFODBDATA_VERSION:-1.2.11}"

if [[ ! -d "${xenofilter_source_directory}" ]]; then
  echo "XenofilteR source directory not found: ${xenofilter_source_directory}" >&2
  exit 1
fi

install_script="$(mktemp --suffix=.R)"
trap 'rm -f "${install_script}"' EXIT

cat > "${install_script}" <<RSCRIPT
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(
  "GenomeInfoDbData",
  version = "${bioconductor_release}",
  ask = FALSE,
  update = FALSE
)
genomeinfodbdata_version <- as.character(utils::packageVersion("GenomeInfoDbData"))
cat("GenomeInfoDbData version:", genomeinfodbdata_version, "\n")
if (genomeinfodbdata_version != "${expected_genomeinfodbdata_version}") {
  stop(paste0(
    "expected GenomeInfoDbData ${expected_genomeinfodbdata_version}",
    " but installed ", genomeinfodbdata_version
  ))
}
if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
  stop("GenomeInfoDb is not loadable after installing GenomeInfoDbData")
}
cat("GenomeInfoDb is loadable\n")
RSCRIPT

"${enva_binary}" run --name "${environment_name}" --command "Rscript '${install_script}'"

"${enva_binary}" run --name "${environment_name}" \
  --command "R CMD INSTALL '${xenofilter_source_directory}'"

"${enva_binary}" run --name "${environment_name}" \
  --command "Rscript -e 'library(XenofilteR); cat(\"XenofilteR\", as.character(packageVersion(\"XenofilteR\")), \"loads\n\")'"
