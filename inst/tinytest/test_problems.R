## "Saved Problems"

rds_files  <- c(
  "degenerate1.RDS",
  "degenerate2.RDS",
  "degenerate3.RDS",
  "degenerate4.RDS",
  "hs21_tiny_qp1.RDS",
  "hs21_tiny_qp2.RDS",
  "infeasible_tiny_qp.RDS",
  "max_ent.RDS",
  "qafiro_tiny_qp1.RDS",
  "qafiro_tiny_qp2.RDS",
  "qafiro_tiny_qp3.RDS",
  "random_prob.RDS",
  "rob_gauss_cov_est.RDS",
  "small_lp.RDS",
  "small_qp1.RDS",
  "unbounded_tiny_qp.RDS"
)

info_vars_to_check  <- c(
  "iter",
  "status",
  "status_val",
  "scale_updates",
  "pobj",
  "dobj",
  ## "res_pri",
  ## "res_dual",
  ## "gap",
  ## "res_infeas",
  "res_unbdd_a",
  "res_unbdd_p",
  ## "setup_time",
  ## "solve_time",
  "scale",
  ## "comp_slack",
  "rejected_accel_steps",
  "accepted_accel_steps"
  ## "lin_sys_time",
  ## "cone_time",
  ## "accel_time"
)

sol_vars_to_check  <- c("x", "y", "s")

for (fname in rds_files) {
  d <- readRDS(file.path("problem_results", fname))
  sol <- scs(A = d$data$A, P = d$data$P, b = d$data$b, obj = d$data$obj, cone = d$cone, control = d$settings)
  expect_equal(sol$info[info_vars_to_check], d$sol$info[info_vars_to_check], tolerance = d$settings$eps_abs)
  expect_equal(sol[sol_vars_to_check], d$sol[sol_vars_to_check], tolerance = d$settings$eps_abs)
}
