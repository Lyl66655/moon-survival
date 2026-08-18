///|
/// Result of a two-group log-rank test.
pub(all) struct LogRankResult {
  statistic : Double
  degrees_of_freedom : Int
  p_value : Double
  observed : Array[Double]
  expected : Array[Double]
  variance : Double
  events : Int
} derive(Eq, Debug)

///|
/// Compare two right-censored groups using the Mantel-Haenszel log-rank
/// statistic.  Ties are handled with the standard finite population
/// correction.
pub fn logrank_test(
  group_a : SurvivalData,
  group_b : SurvivalData,
) -> LogRankResult {
  let groups = [group_a, group_b]
  logrank_test_many(groups)
}

///|
/// General k-group log-rank test.  The returned statistic uses k-1 degrees of
/// freedom and is calculated from the first k-1 group event residuals.
pub fn logrank_test_many(groups : Array[SurvivalData]) -> LogRankResult {
  if groups.length() < 2 {
    panic()
  }
  let times = []
  for group in groups {
    for observation in group.to_array() {
      if observation.event == 1 {
        times.push(observation.time)
      }
    }
  }
  times.sort()
  let unique_times = []
  for time in times {
    if unique_times.length() == 0 ||
      unique_times[unique_times.length() - 1] != time {
      unique_times.push(time)
    }
  }
  let group_count = groups.length()
  let observed = zeros(group_count)
  let expected = zeros(group_count)
  let mut total_events = 0
  let mut statistic = 0.0
  let mut variance = 0.0
  for time in unique_times {
    let risk = []
    let events = []
    let mut total_risk = 0
    let mut total_at_time = 0
    for group in groups {
      let mut group_risk = 0
      let mut group_events = 0
      for observation in group.to_array() {
        if observation.time >= time {
          group_risk += 1
        }
        if observation.time == time && observation.event == 1 {
          group_events += 1
        }
      }
      risk.push(group_risk)
      events.push(group_events)
      total_risk += group_risk
      total_at_time += group_events
    }
    if total_at_time > 0 && total_risk > 0 {
      total_events += total_at_time
      for i in 0..<group_count {
        observed[i] += events[i].to_double()
        expected[i] += total_at_time.to_double() *
          risk[i].to_double() /
          total_risk.to_double()
      }
      if total_risk > 1 {
        let correction = (total_risk - total_at_time).to_double() /
          (total_risk - 1).to_double()
        for i in 0..<group_count {
          let proportion = risk[i].to_double() / total_risk.to_double()
          variance += total_at_time.to_double() *
            correction *
            proportion *
            (1.0 - proportion)
        }
      }
    }
  }
  for i in 0..<(group_count - 1) {
    let residual = observed[i] - expected[i]
    statistic += residual * residual
  }
  if variance > survival_epsilon {
    statistic /= variance
  } else {
    statistic = 0.0
  }
  let degrees = group_count - 1
  let p_value = chi_square_survival(statistic, degrees)
  {
    statistic,
    degrees_of_freedom: degrees,
    p_value,
    observed,
    expected,
    variance,
    events: total_events,
  }
}

///|
/// A one-degree-of-freedom chi-square survival approximation.
pub fn chi_square_survival(
  statistic : Double,
  degrees_of_freedom : Int,
) -> Double {
  if statistic <= 0.0 {
    1.0
  } else if degrees_of_freedom == 1 {
    2.0 * (1.0 - normal_cdf(double_sqrt(statistic)))
  } else {
    double_exp(-0.5 * statistic)
  }
}

///|
/// Compute a standardized log-rank z-score for a two-group comparison.
pub fn LogRankResult::z_score(self : LogRankResult) -> Double {
  if self.observed.length() < 2 || self.variance <= survival_epsilon {
    0.0
  } else {
    (self.observed[0] - self.expected[0]) / double_sqrt(self.variance)
  }
}

///|
/// Return a directionally interpretable hazard contrast.  Values above one
/// indicate more observed events in group A than expected under the null.
pub fn LogRankResult::event_ratio(self : LogRankResult) -> Double {
  if self.expected.length() == 0 || self.expected[0] <= survival_epsilon {
    0.0
  } else {
    self.observed[0] / self.expected[0]
  }
}

///|
/// A pairwise result matrix for exploratory multi-arm screening.
pub fn pairwise_logrank(groups : Array[SurvivalData]) -> Array[Array[Double]] {
  let result = []
  for i in 0..<groups.length() {
    let row = []
    for j in 0..<groups.length() {
      if i == j {
        row.push(1.0)
      } else {
        row.push(logrank_test(groups[i], groups[j]).p_value)
      }
    }
    result.push(row)
  }
  result
}

///|
/// Estimate a restricted mean survival difference between two groups.
pub fn restricted_mean_difference(
  group_a : SurvivalData,
  group_b : SurvivalData,
  horizon : Double,
) -> Double {
  kaplan_meier(group_a).restricted_mean(horizon) -
  kaplan_meier(group_b).restricted_mean(horizon)
}
