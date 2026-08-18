///|
/// Numerically stable online moments using Welford's recurrence.
pub(all) struct OnlineMoments {
  mut count : Int
  mut mean : Double
  mut m2 : Double
  mut minimum : Double
  mut maximum : Double
} derive(Eq, Debug)

///|
/// Create an empty online summary.
pub fn OnlineMoments::new() -> OnlineMoments {
  { count: 0, mean: 0.0, m2: 0.0, minimum: 0.0, maximum: 0.0 }
}

///|
/// Add one value in constant time.
pub fn OnlineMoments::add(self : OnlineMoments, value : Double) -> Unit {
  if value.is_nan() || value.is_inf() {
    panic()
  }
  if self.count == 0 {
    self.count = 1
    self.mean = value
    self.minimum = value
    self.maximum = value
  } else {
    self.count += 1
    let difference = value - self.mean
    self.mean += difference / self.count.to_double()
    let second_difference = value - self.mean
    self.m2 += difference * second_difference
    if value < self.minimum {
      self.minimum = value
    }
    if value > self.maximum {
      self.maximum = value
    }
  }
}

///|
/// Merge two independent online summaries.
pub fn OnlineMoments::merge(
  self : OnlineMoments,
  other : OnlineMoments,
) -> Unit {
  if other.count == 0 {
    return
  }
  if self.count == 0 {
    self.count = other.count
    self.mean = other.mean
    self.m2 = other.m2
    self.minimum = other.minimum
    self.maximum = other.maximum
    return
  }
  let total = self.count + other.count
  let difference = other.mean - self.mean
  self.m2 += other.m2 +
    difference *
    difference *
    self.count.to_double() *
    other.count.to_double() /
    total.to_double()
  self.mean += difference * other.count.to_double() / total.to_double()
  self.count = total
  if other.minimum < self.minimum {
    self.minimum = other.minimum
  }
  if other.maximum > self.maximum {
    self.maximum = other.maximum
  }
}

///|
/// Sample variance of an online summary.
pub fn OnlineMoments::variance(self : OnlineMoments) -> Double {
  if self.count < 2 {
    0.0
  } else {
    self.m2 / (self.count - 1).to_double()
  }
}

///|
/// Convert the online summary to a stable numeric summary.
pub fn OnlineMoments::summary(self : OnlineMoments) -> NumericSummary {
  {
    count: self.count,
    mean: self.mean,
    variance: self.variance(),
    standard_deviation: double_sqrt(self.variance()),
    minimum: self.minimum,
    first_quartile: self.mean,
    median: self.mean,
    third_quartile: self.mean,
    maximum: self.maximum,
  }
}

///|
/// Online exposure and event-rate tracker for reliability services.
pub(all) struct OnlineEventRate {
  mut exposure : Double
  mut events : Int
  mut censored : Int
} derive(Eq, Debug)

///|
/// Create an empty event-rate tracker.
pub fn OnlineEventRate::new() -> OnlineEventRate {
  { exposure: 0.0, events: 0, censored: 0 }
}

///|
/// Add one observed follow-up time.
pub fn OnlineEventRate::add(
  self : OnlineEventRate,
  time : Double,
  event : Int,
) -> Unit {
  if time < 0.0 || time.is_nan() || time.is_inf() || (event != 0 && event != 1) {
    panic()
  }
  self.exposure += time
  if event == 1 {
    self.events += 1
  } else {
    self.censored += 1
  }
}

///|
/// Current constant hazard estimate.
pub fn OnlineEventRate::rate(self : OnlineEventRate) -> Double {
  if self.exposure <= survival_epsilon {
    0.0
  } else {
    self.events.to_double() / self.exposure
  }
}

///|
/// Current censoring fraction.
pub fn OnlineEventRate::censoring_fraction(self : OnlineEventRate) -> Double {
  let count = self.events + self.censored
  if count == 0 {
    0.0
  } else {
    self.censored.to_double() / count.to_double()
  }
}
