

writeLines(getwd())
writeLines(dir())

# Replace printf with scs_printf
file <- "src/scs/src/aa.c"
src <- readLines(file)
src <- gsub("^(\\s*)printf\\(", "\\1scs_printf(", src)
changed_src <- c("Changed 'aa.c'", grep("printf", src, value = TRUE))
# writeLines(changed_src)
writeLines(src, file)


# Replace sprintf with snprintf
file <- "src/scs/src/cones.c"
src <- paste(readLines(file), collapse = "\n")
m <- gregexpr("\\ssprintf\\(.*?;", src)
smatches <- regmatches(src, m)[[1]]
patched <- gsub("sprintf(", "snprintf(", smatches, fixed = TRUE)
patched <- strsplit(patched, ",", fixed = "TRUE")
for (i in seq_along(smatches)) {
  x <- patched[[i]]
  if (i == 1) {
    replacement <- paste(c(x[1], " 512", tail(x, -1)), collapse = ",")
  } else {
    replacement <- paste(c(x[1], " 512 - strlen(tmp)", tail(x, -1)), collapse = ",")
  }
  src <- gsub(smatches[i], replacement, src, fixed = TRUE)
}
writeLines(src, file)



