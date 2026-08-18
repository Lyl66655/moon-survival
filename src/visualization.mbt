///|
/// Render a numeric curve as a compact Unicode-free sparkline for terminals.
pub fn sparkline(values : Array[Double], width? : Int = 60) -> String {
  if values.length() == 0 || width <= 0 {
    ""
  } else {
    let mut minimum = values[0]
    let mut maximum = values[0]
    for value in values {
      if value < minimum {
        minimum = value
      }
      if value > maximum {
        maximum = value
      }
    }
    let symbols = " .:-=+*#%@"
    let builder = StringBuilder::new()
    let points = if values.length() < width { values.length() } else { width }
    for i in 0..<points {
      let index = if points == 1 {
        0
      } else {
        i * (values.length() - 1) / (points - 1)
      }
      let ratio = if maximum - minimum <= survival_epsilon {
        0.5
      } else {
        (values[index] - minimum) / (maximum - minimum)
      }
      let symbol_index = (ratio * 9.0).to_int()
      builder.write_char(symbols.get_char(symbol_index).unwrap_or(' '))
    }
    builder.to_string()
  }
}

///|
/// Render a Kaplan-Meier curve as an ASCII plot with axes labels.
pub fn KMResult::to_ascii(self : KMResult, width? : Int = 60) -> String {
  let values = []
  for step in self.steps {
    values.push(step.survival)
  }
  "survival " + sparkline(values, width~) + "\n"
}

///|
/// Convert event-table rows to a Markdown table.
pub fn event_table_to_markdown(rows : Array[EventTableRow]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "| Time | At risk | Events | Censored |\n|---:|---:|---:|---:|\n",
  )
  for row in rows {
    builder.write_string(
      "| \{row.time} | \{row.n_at_risk} | \{row.n_events} | \{row.n_censored} |\n",
    )
  }
  builder.to_string()
}

///|
/// Create a Markdown section from a competing-risk estimate.
pub fn CIFResult::to_markdown(self : CIFResult) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "| Time | Survival | Cause | Incidence |\n|---:|---:|---:|---:|\n",
  )
  for step in self.steps {
    for cause in 0..<self.n_causes {
      builder.write_string(
        "| \{step.time} | \{step.survival} | \{cause + 1} | \{step.cumulative_incidence[cause]} |\n",
      )
    }
  }
  builder.to_string()
}

///|
/// Format a vector of doubles as a single-line CSV field.
pub fn doubles_to_csv(values : Array[Double]) -> String {
  let builder = StringBuilder::new()
  for i in 0..<values.length() {
    if i > 0 {
      builder.write_string(",")
    }
    builder.write_string("\{values[i]}")
  }
  builder.to_string()
}
