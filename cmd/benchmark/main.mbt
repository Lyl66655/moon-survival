///|
fn main {
  let fixture = @src.benchmark_fixture(5000)
  let start_km = @env.now()
  let km = @src.kaplan_meier(fixture.data)
  let end_km = @env.now()
  let start_na = @env.now()
  let na = @src.nelson_aalen(fixture.data)
  let end_na = @env.now()
  let cox_sample = @src.simulate_cox(
    500,
    [0.55, -0.35, 0.2],
    0.06,
    0.02,
    20260818U,
  )
  let start_cox = @env.now()
  let fit = @src.cox_fit(cox_sample.data, max_iterations=30)
  let end_cox = @env.now()
  let diagnostics = fixture.diagnostics(30.0)
  let start_parametric = @env.now()
  let parametric = @src.fit_parametric_exponential(fixture.data)
  let end_parametric = @env.now()
  let start_ipcw = @env.now()
  let ipcw = @src.ipcw_weights(fixture.data)
  let ipcw_steps = @src.ipcw_kaplan_meier(ipcw)
  let end_ipcw = @env.now()
  let start_group = @env.now()
  let groups = @src.survival_group_km(
    [fixture.data, fixture.data],
    ["replicated_a", "replicated_b"],
    30.0,
  )
  let end_group = @env.now()
  println("MoonBit survival benchmark fixture")
  println("sample_size=\{fixture.data.length()}")
  println("events=\{fixture.data.event_count()}")
  println("censoring_fraction=\{diagnostics.censoring_fraction}")
  println(
    "kaplan_meier_steps=\{km.steps.length()},elapsed_ms=\{end_km - start_km}",
  )
  println(
    "nelson_aalen_steps=\{na.steps.length()},elapsed_ms=\{end_na - start_na}",
  )
  println(
    "cox_n=\{cox_sample.data.length()},iterations=\{fit.iterations},converged=\{fit.converged},elapsed_ms=\{end_cox - start_cox}",
  )
  println("cox_coefficient_count=\{fit.coefficients.length()}")
  println("rmst_30=\{diagnostics.observed_rmst}")
  println(
    "parametric_family=\{@src.parametric_family_name(parametric.family)},elapsed_ms=\{end_parametric - start_parametric}",
  )
  println(
    "ipcw_rows=\{ipcw.length()},ipcw_steps=\{ipcw_steps.length()},effective_n=\{ipcw.effective_n()},elapsed_ms=\{end_ipcw - start_ipcw}",
  )
  println(
    "group_curves=\{groups.length()},group_survival_30=\{groups[0].survival_at_horizon},elapsed_ms=\{end_group - start_group}",
  )
}
