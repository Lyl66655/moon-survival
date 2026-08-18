///|
fn main {
  let data = @src.SurvivalData::new([2.0, 4.0, 5.0, 7.0, 8.0, 10.0], [
    1, 0, 1, 1, 0, 1,
  ])
  let km = @src.kaplan_meier(data)
  let na = @src.nelson_aalen(data)
  let reliability = @src.reliability_report(data, 10.0)
  println("MoonBit survival analysis example")
  println("n=\{data.length()}, events=\{data.event_count()}")
  println("S(6)=\{km.survival_at(6.0)}")
  println("H(6)=\{na.cumulative_hazard_at(6.0)}")
  println("RMST(10)=\{reliability.empirical_rmst}")
}
