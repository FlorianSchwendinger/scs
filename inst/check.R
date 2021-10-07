library(scs)
library(Matrix)
rds_files  <- c(
  "degenerate1.RDS",
  "degenerate2.RDS",
  "degenerate3.RDS",
  "degenerate4.RDS",
  "small_lp.RDS",
  "rob_gauss_cov_est1.RDS",
  "rob_gauss_cov_est2.RDS",
  "hs21_tiny_qp1.RDS",
  "hs21_tiny_qp2.RDS",
  "qafiro_tiny_qp1.RDS",
  "qafiro_tiny_qp2.RDS",
  "qafiro_tiny_qp3.RDS",
  "infeasible_tiny_qp.RDS",
  "unbounded_tiny_qp.RDS",
  "random_prob.RDS"
)

info_vars_to_check  <- c(
  "iter",
  "status",
  "statusVal",
  "scaleUpdates",
  "pobj",
  "dobj",
  ## "resPri",
  ## "resDual",
  ## "gap",
  "resInfeas",
  "resUnbddA",
  "resUnbddP",
  ## "setupTime",
  ## "solveTime",
  "scale",
  ## "compSlack",
  "rejectedAccelSteps",
  "acceptedAccelSteps"
  ## "linsysTime",
  ## "coneTime",
  ## "accelTime"
)

sol_vars_to_check  <- c("x", "y", "s")

fname  <- rds_files[10]
d  <- readRDS(file.path("../tests/testthat/problem_results", fname))
sol <- scs(A = d$data$A, P = d$data$P, b = d$data$b, obj = d$data$obj, initial = d$sol, cone = d$cone, control = d$settings)

