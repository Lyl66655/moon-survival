///|
/// A group count and event summary.
pub(all) struct SurvivalGroupCount {
  group : Int
  sample_size : Int
  events : Int
  censored : Int
  event_fraction : Double
  mean_time : Double
  group_count_row : Bool
} derive(Eq, Debug)

///|
/// A group-level survival contrast.
pub(all) struct SurvivalGroupContrast {
  first_group : Int
  second_group : Int
  horizon : Double
  first_survival : Double
  second_survival : Double
  survival_difference : Double
  rmst_difference : Double
  logrank_p_value : Double
} derive(Eq, Debug)

///|
/// Count groups in a sorted survival sample.
pub fn survival_group_counts(
  data : SurvivalData,
  groups : Array[Int],
) -> Array[SurvivalGroupCount] {
  let result : Array[SurvivalGroupCount] = []
  if groups.length() != data.length() {
    return result
  }
  let unique : Array[Int] = []
  for group in groups {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  unique.sort()
  for group in unique {
    let mut sample_size = 0
    let mut events = 0
    let mut time_total = 0.0
    for i in 0..<groups.length() {
      if groups[i] == group {
        sample_size += 1
        events += data.get(i).event
        time_total += data.get(i).time
      }
    }
    result.push({
      group,
      sample_size,
      events,
      censored: sample_size - events,
      event_fraction: if sample_size == 0 {
        0.0
      } else {
        events.to_double() / sample_size.to_double()
      },
      mean_time: if sample_size == 0 {
        0.0
      } else {
        time_total / sample_size.to_double()
      },
      group_count_row: true,
    })
  }
  result
}

///|
/// Extract one group as core data.
pub fn survival_group_extract(
  data : SurvivalData,
  groups : Array[Int],
  target : Int,
) -> SurvivalData {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let count = if groups.length() < data.length() {
    groups.length()
  } else {
    data.length()
  }
  for i in 0..<count {
    if groups[i] == target {
      times.push(data.get(i).time)
      events.push(data.get(i).event)
    }
  }
  SurvivalData::new(times, events)
}

///|
/// Calculate one pairwise group contrast.
pub fn survival_group_contrast(
  data : SurvivalData,
  groups : Array[Int],
  first_group : Int,
  second_group : Int,
  horizon : Double,
) -> SurvivalGroupContrast {
  let first = survival_group_extract(data, groups, first_group)
  let second = survival_group_extract(data, groups, second_group)
  let first_fit = kaplan_meier(first)
  let second_fit = kaplan_meier(second)
  let logrank = logrank_test(first, second)
  {
    first_group,
    second_group,
    horizon,
    first_survival: first_fit.survival_at(horizon),
    second_survival: second_fit.survival_at(horizon),
    survival_difference: first_fit.survival_at(horizon) -
    second_fit.survival_at(horizon),
    rmst_difference: first_fit.restricted_mean(horizon) -
    second_fit.restricted_mean(horizon),
    logrank_p_value: logrank.p_value,
  }
}

///|
/// Calculate all pairwise contrasts.
pub fn survival_all_group_contrasts(
  data : SurvivalData,
  groups : Array[Int],
  horizon : Double,
) -> Array[SurvivalGroupContrast] {
  let result : Array[SurvivalGroupContrast] = []
  let unique : Array[Int] = []
  for group in groups {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  unique.sort()
  for i in 0..<unique.length() {
    for j in (i + 1)..<unique.length() {
      result.push(
        survival_group_contrast(data, groups, unique[i], unique[j], horizon),
      )
    }
  }
  result
}

///|
/// Compute a standardized group survival curve.
pub fn survival_group_curve(
  data : SurvivalData,
  groups : Array[Int],
  group : Int,
  times : Array[Double],
) -> Array[CurvePoint] {
  let fit = kaplan_meier(survival_group_extract(data, groups, group))
  let result : Array[CurvePoint] = []
  for time in times {
    result.push({ time, value: fit.survival_at(time) })
  }
  result
}

///|
/// Compute a weighted group mean of a numeric value.
pub fn survival_group_weighted_mean(
  values : Array[Double],
  groups : Array[Int],
  weights : Array[Double],
  group : Int,
) -> Double {
  let selected_values : Array[Double] = []
  let selected_weights : Array[Double] = []
  let count = if values.length() < groups.length() {
    values.length()
  } else {
    groups.length()
  }
  for i in 0..<count {
    if groups[i] == group {
      selected_values.push(values[i])
      selected_weights.push(if i < weights.length() { weights[i] } else { 1.0 })
    }
  }
  survival_ext_weighted_mean(selected_values, selected_weights)
}

///|
/// Compute a standardized mean difference between two groups.
pub fn survival_group_standardized_difference(
  values : Array[Double],
  groups : Array[Int],
  first_group : Int,
  second_group : Int,
) -> Double {
  let first : Array[Double] = []
  let second : Array[Double] = []
  for i in 0..<values.length() {
    if i < groups.length() && groups[i] == first_group {
      first.push(values[i])
    } else if i < groups.length() && groups[i] == second_group {
      second.push(values[i])
    }
  }
  let pooled = double_sqrt(
    (survival_ext_variance(first) + survival_ext_variance(second)) / 2.0,
  )
  survival_ext_safe_ratio(
    survival_ext_mean(first) - survival_ext_mean(second),
    pooled,
  )
}

///|
/// Return group indices for a target.
pub fn survival_group_indices(groups : Array[Int], target : Int) -> Array[Int] {
  let result : Array[Int] = []
  for i in 0..<groups.length() {
    if groups[i] == target {
      result.push(i)
    }
  }
  result
}

///|
/// Return a group vector's unique levels.
pub fn survival_group_levels(groups : Array[Int]) -> Array[Int] {
  let result : Array[Int] = []
  for group in groups {
    if !result.contains(group) {
      result.push(group)
    }
  }
  result.sort()
  result
}

///|
/// Compute a group-wise event rate through a horizon.
pub fn survival_group_event_rates(
  data : SurvivalData,
  groups : Array[Int],
  horizon : Double,
) -> Array[Double] {
  let levels = survival_group_levels(groups)
  let result : Array[Double] = []
  for group in levels {
    let subset = survival_group_extract(data, groups, group)
    result.push(1.0 - kaplan_meier(subset).survival_at(horizon))
  }
  result
}

///|
/// Export group counts.
pub fn survival_group_counts_to_csv(
  counts : Array[SurvivalGroupCount],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "group,sample_size,events,censored,event_fraction,mean_time\n",
  )
  for count in counts {
    builder.write_string(count.group.to_string())
    builder.write_char(',')
    builder.write_string(count.sample_size.to_string())
    builder.write_char(',')
    builder.write_string(count.events.to_string())
    builder.write_char(',')
    builder.write_string(count.censored.to_string())
    builder.write_char(',')
    builder.write_string(count.event_fraction.to_string())
    builder.write_char(',')
    builder.write_string(count.mean_time.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export group contrasts.
pub fn survival_group_contrasts_to_csv(
  contrasts : Array[SurvivalGroupContrast],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "first_group,second_group,horizon,first_survival,second_survival,survival_difference,rmst_difference,logrank_p_value\n",
  )
  for contrast in contrasts {
    builder.write_string(contrast.first_group.to_string())
    builder.write_char(',')
    builder.write_string(contrast.second_group.to_string())
    builder.write_char(',')
    builder.write_string(contrast.horizon.to_string())
    builder.write_char(',')
    builder.write_string(contrast.first_survival.to_string())
    builder.write_char(',')
    builder.write_string(contrast.second_survival.to_string())
    builder.write_char(',')
    builder.write_string(contrast.survival_difference.to_string())
    builder.write_char(',')
    builder.write_string(contrast.rmst_difference.to_string())
    builder.write_char(',')
    builder.write_string(contrast.logrank_p_value.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
