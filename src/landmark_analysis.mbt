///|
/// Landmark survival estimate for subjects alive and uncensored at a time.
pub(all) struct LandmarkEstimate {
  landmark : Double
  eligible : Int
  subsequent_events : Int
  horizon : Double
  survival : Double
  restricted_mean : Double
} derive(Eq, Debug)

///|
/// Rebase follow-up at a landmark and estimate conditional survival.
pub fn landmark_estimate(
  data : SurvivalData,
  landmark : Double,
  horizon : Double,
) -> LandmarkEstimate {
  let times = []
  let events = []
  for observation in data.to_array() {
    if observation.time > landmark {
      times.push(observation.time - landmark)
      events.push(observation.event)
    }
  }
  let rebased = SurvivalData::new(times, events)
  let km = kaplan_meier(rebased)
  {
    landmark,
    eligible: rebased.length(),
    subsequent_events: rebased.event_count(),
    horizon,
    survival: km.predict(horizon),
    restricted_mean: km.restricted_mean(horizon),
  }
}

///|
/// Difference in restricted mean survival between two fitted curves.
pub fn restricted_mean_loss(
  baseline : SurvivalData,
  follow_up : SurvivalData,
  horizon : Double,
) -> Double {
  kaplan_meier(baseline).restricted_mean(horizon) -
  kaplan_meier(follow_up).restricted_mean(horizon)
}

///|
/// A discrete hazard row generated from a KM event table.
pub(all) struct DiscreteHazardRow {
  time : Double
  interval_width : Double
  n_at_risk : Int
  events : Int
  hazard : Double
  survival_after : Double
} derive(Eq, Debug)

///|
/// Convert event-table steps into an application-friendly hazard table.
pub fn discrete_hazard_table(data : SurvivalData) -> Array[DiscreteHazardRow] {
  let rows = event_table(data)
  let result = []
  let mut previous = 0.0
  let mut survival = 1.0
  for row in rows {
    let hazard = if row.n_at_risk == 0 {
      0.0
    } else {
      row.n_events.to_double() / row.n_at_risk.to_double()
    }
    survival *= 1.0 - hazard
    result.push({
      time: row.time,
      interval_width: row.time - previous,
      n_at_risk: row.n_at_risk,
      events: row.n_events,
      hazard,
      survival_after: survival,
    })
    previous = row.time
  }
  result
}

///|
/// Convert a landmark result to one CSV row.
pub fn LandmarkEstimate::to_csv(self : LandmarkEstimate) -> String {
  "landmark,eligible,subsequent_events,horizon,survival,restricted_mean\n" +
  "\{self.landmark},\{self.eligible},\{self.subsequent_events},\{self.horizon},\{self.survival},\{self.restricted_mean}\n"
}
