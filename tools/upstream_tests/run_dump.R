## Run upstream SCS 3.3.1's test functions inside this R process.
## Usage: Rscript run_dump.R <upstream tree root> <dump_tests.so> <out dir>
args <- commandArgs(trailingOnly = TRUE)
tree <- args[1]; so <- args[2]; out <- args[3]
setwd(tree)          ## tests open test/problems/<data> relative to the tree root
dyn.load(so)
res <- .Call("dump_tests_run", out)
cat(sprintf("\nR process: %d upstream tests run, %d assertion failures\n", res[1], res[2]))
quit(status = if (res[2] > 0) 2L else 0L)
