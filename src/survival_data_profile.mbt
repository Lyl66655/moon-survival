///|
/// A compact data profile for a survival dataset.
pub(all) struct SurvivalDataProfile {
  rows : Int
  events : Int
  censored : Int
  minimum_time : Double
  maximum_time : Double
  mean_time : Double
  event_fraction : Double
  profile_row : Bool
} derive(Eq, Debug)

///|
/// Create a data profile.
pub fn survival_data_profile(data : SurvivalData) -> SurvivalDataProfile {
  let times : Array[Double] = []
  for observation in data.to_array() {
    times.push(observation.time)
  }
  {
    rows: data.length(),
    events: data.event_count(),
    censored: data.censored_count(),
    minimum_time: survival_ext_min(times),
    maximum_time: survival_ext_max(times),
    mean_time: survival_ext_mean(times),
    event_fraction: if data.length() == 0 {
      0.0
    } else {
      data.event_count().to_double() / data.length().to_double()
    },
    profile_row: true,
  }
}

///|
/// Return the event-time profile.
pub fn survival_event_time_profile(data : SurvivalData) -> SurvivalDataProfile {
  survival_data_profile(
    SurvivalData::new(
      data.event_times(),
      Array::make(data.event_times().length(), 1),
    ),
  )
}

///|
/// Return the censoring-time profile.
pub fn survival_censoring_time_profile(
  data : SurvivalData,
) -> SurvivalDataProfile {
  survival_data_profile(
    SurvivalData::new(
      data.censoring_times(),
      Array::make(data.censoring_times().length(), 0),
    ),
  )
}

///|
/// Compare two profiles.
pub fn survival_profile_difference(
  first : SurvivalDataProfile,
  second : SurvivalDataProfile,
) -> SurvivalDataProfile {
  {
    rows: first.rows - second.rows,
    events: first.events - second.events,
    censored: first.censored - second.censored,
    minimum_time: first.minimum_time - second.minimum_time,
    maximum_time: first.maximum_time - second.maximum_time,
    mean_time: first.mean_time - second.mean_time,
    event_fraction: first.event_fraction - second.event_fraction,
    profile_row: true,
  }
}

///|
/// Export a data profile.
pub fn survival_data_profile_to_csv(profile : SurvivalDataProfile) -> String {
  "metric,value\nrows," +
  profile.rows.to_string() +
  "\nevents," +
  profile.events.to_string() +
  "\ncensored," +
  profile.censored.to_string() +
  "\nminimum_time," +
  profile.minimum_time.to_string() +
  "\nmaximum_time," +
  profile.maximum_time.to_string() +
  "\nmean_time," +
  profile.mean_time.to_string() +
  "\nevent_fraction," +
  profile.event_fraction.to_string() +
  "\n"
}

///|
/// Validate profile ranges.
pub fn survival_data_profile_is_valid(profile : SurvivalDataProfile) -> Bool {
  profile.rows >= profile.events &&
  profile.events >= 0 &&
  profile.censored >= 0 &&
  profile.minimum_time >= 0.0 &&
  profile.maximum_time >= profile.minimum_time &&
  profile.event_fraction >= 0.0 &&
  profile.event_fraction <= 1.0
}
