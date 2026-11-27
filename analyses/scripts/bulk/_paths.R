script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
chapter_dir <- normalizePath(file.path(dirname(script_file), "../../.."))
bulk_out <- file.path(chapter_dir, "output/bulk")
