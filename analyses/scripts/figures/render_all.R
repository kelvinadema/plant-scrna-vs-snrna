.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(.script_file)), "_setup.R"))
scripts <- list.files(dirname(normalizePath(.script_file)), pattern = "^(fig|table).*\\.R$", full.names = TRUE)
for (script in scripts) {
  message(basename(script))
  source(script, local = new.env(parent = globalenv()))
}
