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
}
