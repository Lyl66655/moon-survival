///|
/// A subject-by-endpoint outcome matrix.
pub(all) struct SurvivalEndpointRow {
  subject : Int
  outcomes : Array[Int]
  earliest_event : Double
  endpoint_row : Bool
} derive(Eq, Debug)

///|
/// An endpoint matrix with named endpoints.
pub(all) struct SurvivalEndpointMatrix {
  endpoints : Array[String]
  rows : Array[SurvivalEndpointRow]
  endpoint_matrix : Bool
} derive(Eq, Debug)

///|
/// Build a binary endpoint matrix from event-code columns.
pub fn survival_endpoint_matrix(
  names : Array[String],
  observations : Array[Array[Int]],
  times : Array[Double],
) -> SurvivalEndpointMatrix {
  let rows : Array[SurvivalEndpointRow] = []
  for i in 0..<observations.length() {
    let time = if i < times.length() { times[i] } else { 0.0 }
    rows.push({
      subject: i,
      outcomes: observations[i],
      earliest_event: time,
      endpoint_row: true,
    })
  }
  { endpoints: names, rows, endpoint_matrix: true }
}

///|
/// Return the number of subjects.
pub fn SurvivalEndpointMatrix::subject_count(
  self : SurvivalEndpointMatrix,
) -> Int {
  self.rows.length()
}

///|
/// Return the number of endpoints.
pub fn SurvivalEndpointMatrix::endpoint_count(
  self : SurvivalEndpointMatrix,
) -> Int {
  self.endpoints.length()
}

///|
/// Count positive outcomes for one endpoint.
pub fn SurvivalEndpointMatrix::event_count(
  self : SurvivalEndpointMatrix,
  endpoint : Int,
) -> Int {
  let mut count = 0
  for row in self.rows {
    if endpoint >= 0 &&
      endpoint < row.outcomes.length() &&
      row.outcomes[endpoint] == 1 {
      count += 1
    }
  }
  count
}

///|
/// Return incidence for one endpoint.
pub fn SurvivalEndpointMatrix::incidence(
  self : SurvivalEndpointMatrix,
  endpoint : Int,
) -> Double {
  if self.subject_count() == 0 {
    0.0
  } else {
    self.event_count(endpoint).to_double() / self.subject_count().to_double()
  }
}

///|
/// Return a composite endpoint flag.
pub fn survival_endpoint_composite(
  matrix : SurvivalEndpointMatrix,
  endpoints : Array[Int],
) -> Array[Int] {
  let result : Array[Int] = []
  for row in matrix.rows {
    let mut positive = 0
    for endpoint in endpoints {
      if endpoint >= 0 &&
        endpoint < row.outcomes.length() &&
        row.outcomes[endpoint] == 1 {
        positive = 1
      }
    }
    result.push(positive)
  }
  result
}

///|
/// Return all endpoint flags for one subject.
pub fn SurvivalEndpointMatrix::subject_outcomes(
  self : SurvivalEndpointMatrix,
  subject : Int,
) -> Array[Int] {
  for row in self.rows {
    if row.subject == subject {
      return row.outcomes
    }
  }
  []
}

///|
/// Build core survival data for a selected endpoint.
pub fn survival_endpoint_to_data(
  matrix : SurvivalEndpointMatrix,
  endpoint : Int,
  times : Array[Double],
) -> SurvivalData {
  let selected_times : Array[Double] = []
  let events : Array[Int] = []
  for i in 0..<matrix.rows.length() {
    selected_times.push(
      if i < times.length() {
        times[i]
      } else {
        matrix.rows[i].earliest_event
      },
    )
    events.push(
      if endpoint >= 0 && endpoint < matrix.rows[i].outcomes.length() {
        matrix.rows[i].outcomes[endpoint]
      } else {
        0
      },
    )
  }
  SurvivalData::new(selected_times, events)
}

///|
/// Compute endpoint co-occurrence counts.
pub fn survival_endpoint_cooccurrence(
  matrix : SurvivalEndpointMatrix,
) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = []
  for first in 0..<matrix.endpoint_count() {
    let row = Array::make(matrix.endpoint_count(), 0)
    for second in 0..<matrix.endpoint_count() {
      for subject in matrix.rows {
        if first < subject.outcomes.length() &&
          second < subject.outcomes.length() &&
          subject.outcomes[first] == 1 &&
          subject.outcomes[second] == 1 {
          row[second] += 1
        }
      }
    }
    result.push(row)
  }
  result
}

///|
/// Compute mutually exclusive endpoint categories.
pub fn survival_endpoint_category(
  matrix : SurvivalEndpointMatrix,
) -> Array[Int] {
  let result : Array[Int] = []
  for row in matrix.rows {
    let mut category = -1
    for endpoint in 0..<row.outcomes.length() {
      if row.outcomes[endpoint] == 1 && category == -1 {
        category = endpoint
      }
    }
    result.push(category)
  }
  result
}

///|
/// Return rows complete across all endpoints.
pub fn survival_endpoint_complete_rows(
  matrix : SurvivalEndpointMatrix,
) -> Array[Int] {
  let result : Array[Int] = []
  for row in matrix.rows {
    let mut complete = true
    for outcome in row.outcomes {
      if outcome != 0 && outcome != 1 {
        complete = false
      }
    }
    if complete {
      result.push(row.subject)
    }
  }
  result
}

///|
/// Calculate the endpoint correlation matrix using binary correlations.
pub fn survival_endpoint_correlation(
  matrix : SurvivalEndpointMatrix,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for first in 0..<matrix.endpoint_count() {
    let row : Array[Double] = []
    for second in 0..<matrix.endpoint_count() {
      let mut same = 0.0
      for subject in matrix.rows {
        if first < subject.outcomes.length() &&
          second < subject.outcomes.length() &&
          subject.outcomes[first] == subject.outcomes[second] {
          same += 1.0
        }
      }
      row.push(
        if matrix.subject_count() == 0 {
          0.0
        } else {
          same / matrix.subject_count().to_double()
        },
      )
    }
    result.push(row)
  }
  result
}

///|
/// Export endpoint matrix rows.
pub fn survival_endpoint_matrix_to_csv(
  matrix : SurvivalEndpointMatrix,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,")
  for i in 0..<matrix.endpoints.length() {
    if i > 0 {
      builder.write_char(',')
    }
    builder.write_string(matrix.endpoints[i])
  }
  builder.write_char('\n')
  for row in matrix.rows {
    builder.write_string(row.subject.to_string())
    for outcome in row.outcomes {
      builder.write_char(',')
      builder.write_string(outcome.to_string())
    }
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate endpoint dimensions and coding.
pub fn survival_endpoint_matrix_is_valid(
  matrix : SurvivalEndpointMatrix,
) -> Bool {
  for row in matrix.rows {
    if row.outcomes.length() != matrix.endpoints.length() {
      return false
    }
    for outcome in row.outcomes {
      if outcome != 0 && outcome != 1 {
        return false
      }
    }
  }
  true
}
