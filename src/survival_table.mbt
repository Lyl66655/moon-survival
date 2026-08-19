///|
/// A typed row for analysis-ready survival tables.
pub(all) struct SurvivalTableRow {
  row_id : Int
  time : Double
  event : Int
  group : Int
  weight : Double
  table_row : Bool
} derive(Eq, Debug)

///|
/// A column-oriented table with stable row identifiers.
pub(all) struct SurvivalTable {
  rows : Array[SurvivalTableRow]
  group_count : Int
  table_version : String
} derive(Eq, Debug)

///|
/// Errors returned while constructing an analysis table.
pub(all) enum SurvivalTableError {
  EmptyTable
  LengthMismatch
  InvalidTime
  InvalidEvent
  InvalidGroup
  InvalidWeight
  DuplicateRow
} derive(Eq, Debug)

///|
/// Group-level summary used by descriptive dashboards.
pub(all) struct SurvivalTableGroupSummary {
  group : Int
  sample_size : Int
  events : Int
  censoring_rate : Double
  weighted_events : Double
  total_weight : Double
  mean_time : Double
  median_time : Double
} derive(Eq, Debug)

///|
/// A time-window aggregate for operational monitoring.
pub(all) struct SurvivalWindowSummary {
  start : Double
  end : Double
  at_risk : Int
  entries : Int
  exits : Int
  events : Int
  exposure : Double
  event_rate : Double
  window_summary : Bool
} derive(Eq, Debug)

///|
/// Construct a validated analysis table.
pub fn SurvivalTable::try_new(
  times : Array[Double],
  events : Array[Int],
  groups : Array[Int],
  weights : Array[Double],
) -> Result[SurvivalTable, SurvivalTableError] {
  if times.length() == 0 {
    return Err(EmptyTable)
  }
  if times.length() != events.length() ||
    times.length() != groups.length() ||
    times.length() != weights.length() {
    return Err(LengthMismatch)
  }
  let rows : Array[SurvivalTableRow] = []
  let mut maximum_group = -1
  for i in 0..<times.length() {
    if times[i].is_nan() || times[i].is_inf() || times[i] < 0.0 {
      return Err(InvalidTime)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent)
    }
    if groups[i] < 0 {
      return Err(InvalidGroup)
    }
    if weights[i].is_nan() || weights[i].is_inf() || weights[i] <= 0.0 {
      return Err(InvalidWeight)
    }
    if groups[i] > maximum_group {
      maximum_group = groups[i]
    }
    rows.push({
      row_id: i,
      time: times[i],
      event: events[i],
      group: groups[i],
      weight: weights[i],
      table_row: true,
    })
  }
  rows.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else if a.group < b.group {
      -1
    } else {
      a.row_id.compare(b.row_id)
    }
  })
  Ok({ rows, group_count: maximum_group + 1, table_version: "1" })
}

///|
/// Construct a unit-weight, single-group table.
pub fn survival_table_from_data(data : SurvivalData) -> SurvivalTable {
  if data.length() == 0 {
    return survival_table_empty()
  }
  let times : Array[Double] = []
  let events : Array[Int] = []
  let groups : Array[Int] = []
  let weights : Array[Double] = []
  for observation in data.to_array() {
    times.push(observation.time)
    events.push(observation.event)
    groups.push(0)
    weights.push(1.0)
  }
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Convenience constructor for a table with one group and unit weights.
pub fn SurvivalTable::from_arrays(
  times : Array[Double],
  events : Array[Int],
) -> SurvivalTable {
  let groups = Array::make(times.length(), 0)
  let weights = Array::make(times.length(), 1.0)
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Return a typed empty result for filters that match no rows.
fn survival_table_empty() -> SurvivalTable {
  { rows: [], group_count: 0, table_version: "1" }
}

///|
/// Number of rows.
pub fn SurvivalTable::length(self : SurvivalTable) -> Int {
  self.rows.length()
}

///|
/// Number of represented groups.
pub fn SurvivalTable::group_count(self : SurvivalTable) -> Int {
  self.group_count
}

///|
/// Access a sorted row.
pub fn SurvivalTable::get(
  self : SurvivalTable,
  index : Int,
) -> SurvivalTableRow {
  self.rows[index]
}

///|
/// Copy all rows.
pub fn SurvivalTable::to_array(self : SurvivalTable) -> Array[SurvivalTableRow] {
  let result : Array[SurvivalTableRow] = []
  for row in self.rows {
    result.push(row)
  }
  result
}

///|
/// Return row ids in sorted analysis order.
pub fn SurvivalTable::row_ids(self : SurvivalTable) -> Array[Int] {
  let result : Array[Int] = []
  for row in self.rows {
    result.push(row.row_id)
  }
  result
}

///|
/// Filter rows by group.
pub fn SurvivalTable::filter_group(
  self : SurvivalTable,
  group : Int,
) -> SurvivalTable {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let groups : Array[Int] = []
  let weights : Array[Double] = []
  for row in self.rows {
    if row.group == group {
      times.push(row.time)
      events.push(row.event)
      groups.push(0)
      weights.push(row.weight)
    }
  }
  if times.length() == 0 {
    return survival_table_empty()
  }
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Filter by a closed time interval.
pub fn SurvivalTable::filter_time(
  self : SurvivalTable,
  start : Double,
  end : Double,
) -> SurvivalTable {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let groups : Array[Int] = []
  let weights : Array[Double] = []
  for row in self.rows {
    if row.time >= start && row.time <= end {
      times.push(row.time)
      events.push(row.event)
      groups.push(row.group)
      weights.push(row.weight)
    }
  }
  if times.length() == 0 {
    return survival_table_empty()
  }
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Return rows satisfying a minimum weight threshold.
pub fn SurvivalTable::filter_weight(
  self : SurvivalTable,
  minimum : Double,
) -> SurvivalTable {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let groups : Array[Int] = []
  let weights : Array[Double] = []
  for row in self.rows {
    if row.weight >= minimum {
      times.push(row.time)
      events.push(row.event)
      groups.push(row.group)
      weights.push(row.weight)
    }
  }
  if times.length() == 0 {
    return survival_table_empty()
  }
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Convert a table to the core unweighted survival data.
pub fn SurvivalTable::to_survival_data(self : SurvivalTable) -> SurvivalData {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for row in self.rows {
    times.push(row.time)
    events.push(row.event)
  }
  SurvivalData::new(times, events)
}

///|
/// Convert one group to core survival data.
pub fn SurvivalTable::to_group_data(
  self : SurvivalTable,
  group : Int,
) -> SurvivalData {
  self.filter_group(group).to_survival_data()
}

///|
/// Return weighted group summaries.
pub fn SurvivalTable::group_summaries(
  self : SurvivalTable,
) -> Array[SurvivalTableGroupSummary] {
  let result : Array[SurvivalTableGroupSummary] = []
  for group in 0..<self.group_count {
    let times : Array[Double] = []
    let mut events = 0
    let mut weighted_events = 0.0
    let mut total_weight = 0.0
    let mut sum_time = 0.0
    for row in self.rows {
      if row.group == group {
        times.push(row.time)
        events += row.event
        weighted_events += row.weight * row.event.to_double()
        total_weight += row.weight
        sum_time += row.time
      }
    }
    times.sort()
    let median = if times.length() == 0 {
      0.0
    } else {
      times[times.length() / 2]
    }
    result.push({
      group,
      sample_size: times.length(),
      events,
      censoring_rate: if times.length() == 0 {
        0.0
      } else {
        1.0 - events.to_double() / times.length().to_double()
      },
      weighted_events,
      total_weight,
      mean_time: if times.length() == 0 {
        0.0
      } else {
        sum_time / times.length().to_double()
      },
      median_time: median,
    })
  }
  result
}

///|
/// Count unique rows that are under observation at a query time.
pub fn SurvivalTable::at_risk(self : SurvivalTable, time : Double) -> Int {
  let mut count = 0
  for row in self.rows {
    if row.time >= time {
      count += 1
    }
  }
  count
}

///|
/// Calculate weighted event totals through a horizon.
pub fn SurvivalTable::events_through(
  self : SurvivalTable,
  horizon : Double,
) -> Double {
  let mut total = 0.0
  for row in self.rows {
    if row.event == 1 && row.time <= horizon {
      total += row.weight
    }
  }
  total
}

///|
/// Calculate weighted exposure through a horizon.
pub fn SurvivalTable::exposure_through(
  self : SurvivalTable,
  horizon : Double,
) -> Double {
  let mut total = 0.0
  for row in self.rows {
    total += row.weight * (if row.time < horizon { row.time } else { horizon })
  }
  total
}

///|
/// Build fixed-width time windows for operational reporting.
pub fn SurvivalTable::windows(
  self : SurvivalTable,
  width : Double,
  end : Double,
) -> Array[SurvivalWindowSummary] {
  let result : Array[SurvivalWindowSummary] = []
  if width <= 0.0 || end <= 0.0 {
    return result
  }
  let mut start = 0.0
  while start < end {
    let stop = if start + width > end { end } else { start + width }
    let mut at_risk = 0
    let mut entries = 0
    let mut exits = 0
    let mut events = 0
    let mut exposure = 0.0
    for row in self.rows {
      if row.time >= start {
        at_risk += 1
      }
      if row.time >= start && row.time < stop {
        exits += 1
        events += row.event
      }
      if row.time > start && row.time <= stop {
        entries += 1
      }
      let observed = if row.time < stop {
        row.time - start
      } else {
        stop - start
      }
      if observed > 0.0 {
        exposure += row.weight * observed
      }
    }
    result.push({
      start,
      end: stop,
      at_risk,
      entries,
      exits,
      events,
      exposure,
      event_rate: if exposure <= survival_epsilon {
        0.0
      } else {
        events.to_double() / exposure
      },
      window_summary: true,
    })
    start = stop
  }
  result
}

///|
/// Add a second table, preserving group labels through an offset.
pub fn survival_table_append(
  left : SurvivalTable,
  right : SurvivalTable,
) -> SurvivalTable {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let groups : Array[Int] = []
  let weights : Array[Double] = []
  for row in left.rows {
    times.push(row.time)
    events.push(row.event)
    groups.push(row.group)
    weights.push(row.weight)
  }
  for row in right.rows {
    times.push(row.time)
    events.push(row.event)
    groups.push(row.group + left.group_count)
    weights.push(row.weight)
  }
  match SurvivalTable::try_new(times, events, groups, weights) {
    Ok(table) => table
    Err(_) => panic()
  }
}

///|
/// Compute a weighted mean of a numeric row expression.
pub fn SurvivalTable::weighted_time(self : SurvivalTable) -> Double {
  let mut total = 0.0
  let mut denominator = 0.0
  for row in self.rows {
    total += row.time * row.weight
    denominator += row.weight
  }
  if denominator <= survival_epsilon {
    0.0
  } else {
    total / denominator
  }
}

///|
/// Compute the weighted event fraction.
pub fn SurvivalTable::weighted_event_fraction(self : SurvivalTable) -> Double {
  let mut events = 0.0
  let mut total = 0.0
  for row in self.rows {
    events += row.weight * row.event.to_double()
    total += row.weight
  }
  if total <= survival_epsilon {
    0.0
  } else {
    events / total
  }
}

///|
/// Validate row order and finite values.
pub fn SurvivalTable::is_valid(self : SurvivalTable) -> Bool {
  if self.rows.length() == 0 || self.group_count <= 0 {
    return false
  }
  for row in self.rows {
    if row.row_id < 0 ||
      row.time < 0.0 ||
      row.weight <= 0.0 ||
      (row.event != 0 && row.event != 1) ||
      row.group < 0 {
      return false
    }
  }
  true
}

///|
/// Export the table as CSV.
pub fn SurvivalTable::to_csv(self : SurvivalTable) -> String {
  let builder = StringBuilder::new()
  builder.write_string("row_id,time,event,group,weight\n")
  for row in self.rows {
    builder.write_string(row.row_id.to_string())
    builder.write_char(',')
    builder.write_string(row.time.to_string())
    builder.write_char(',')
    builder.write_string(row.event.to_string())
    builder.write_char(',')
    builder.write_string(row.group.to_string())
    builder.write_char(',')
    builder.write_string(row.weight.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export grouped descriptive statistics.
pub fn survival_table_summaries_to_csv(
  summaries : Array[SurvivalTableGroupSummary],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "group,sample_size,events,censoring_rate,weighted_events,total_weight,mean_time,median_time\n",
  )
  for summary in summaries {
    builder.write_string(summary.group.to_string())
    builder.write_char(',')
    builder.write_string(summary.sample_size.to_string())
    builder.write_char(',')
    builder.write_string(summary.events.to_string())
    builder.write_char(',')
    builder.write_string(summary.censoring_rate.to_string())
    builder.write_char(',')
    builder.write_string(summary.weighted_events.to_string())
    builder.write_char(',')
    builder.write_string(summary.total_weight.to_string())
    builder.write_char(',')
    builder.write_string(summary.mean_time.to_string())
    builder.write_char(',')
    builder.write_string(summary.median_time.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export operational windows.
pub fn survival_windows_to_csv(
  windows : Array[SurvivalWindowSummary],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "start,end,at_risk,entries,exits,events,exposure,event_rate\n",
  )
  for window in windows {
    builder.write_string(window.start.to_string())
    builder.write_char(',')
    builder.write_string(window.end.to_string())
    builder.write_char(',')
    builder.write_string(window.at_risk.to_string())
    builder.write_char(',')
    builder.write_string(window.entries.to_string())
    builder.write_char(',')
    builder.write_string(window.exits.to_string())
    builder.write_char(',')
    builder.write_string(window.events.to_string())
    builder.write_char(',')
    builder.write_string(window.exposure.to_string())
    builder.write_char(',')
    builder.write_string(window.event_rate.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
