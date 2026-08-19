///|
/// Online cumulative monitoring state.
pub(all) struct SurvivalOnlineMonitor {
  observations : Int
  events : Int
  exposure : Double
  cumulative_hazard : Double
  latest_time : Double
  alert : Bool
  monitor_row : Bool
} derive(Eq, Debug)

///|
/// A monitor update result.
pub(all) struct SurvivalMonitorUpdate {
  time : Double
  event : Int
  hazard_increment : Double
  cumulative_hazard : Double
  event_rate : Double
  alert : Bool
} derive(Eq, Debug)

///|
/// Create an empty online monitor.
pub fn survival_online_monitor() -> SurvivalOnlineMonitor {
  {
    observations: 0,
    events: 0,
    exposure: 0.0,
    cumulative_hazard: 0.0,
    latest_time: 0.0,
    alert: false,
    monitor_row: true,
  }
}

///|
/// Update a monitor with one record.
pub fn SurvivalOnlineMonitor::update(
  self : SurvivalOnlineMonitor,
  time : Double,
  event : Int,
  alert_threshold : Double,
) -> (SurvivalOnlineMonitor, SurvivalMonitorUpdate) {
  let elapsed = if time > self.latest_time {
    time - self.latest_time
  } else {
    0.0
  }
  let exposure = self.exposure + elapsed
  let observations = self.observations + 1
  let events = self.events + event
  let increment = if observations == 0 {
    0.0
  } else {
    event.to_double() / observations.to_double()
  }
  let cumulative = self.cumulative_hazard + increment
  let rate = survival_ext_safe_ratio(events.to_double(), exposure)
  let alert = cumulative >= alert_threshold
  let next_monitor : SurvivalOnlineMonitor = {
    observations,
    events,
    exposure,
    cumulative_hazard: cumulative,
    latest_time: time,
    alert,
    monitor_row: true,
  }
  let update : SurvivalMonitorUpdate = {
    time,
    event,
    hazard_increment: increment,
    cumulative_hazard: cumulative,
    event_rate: rate,
    alert,
  }
  (next_monitor, update)
}

///|
/// Update a monitor from a batch.
pub fn survival_online_monitor_batch(
  monitor : SurvivalOnlineMonitor,
  data : SurvivalData,
  alert_threshold : Double,
) -> SurvivalOnlineMonitor {
  let mut current = monitor
  for observation in data.to_array() {
    let (next, _) = current.update(
      observation.time,
      observation.event,
      alert_threshold,
    )
    current = next
  }
  current
}

///|
/// Estimate current online survival.
pub fn SurvivalOnlineMonitor::survival(self : SurvivalOnlineMonitor) -> Double {
  double_exp(-self.cumulative_hazard)
}

///|
/// Estimate current event rate.
pub fn SurvivalOnlineMonitor::event_rate(
  self : SurvivalOnlineMonitor,
) -> Double {
  survival_ext_safe_ratio(self.events.to_double(), self.exposure)
}

///|
/// Return whether an alert is active.
pub fn SurvivalOnlineMonitor::has_alert(self : SurvivalOnlineMonitor) -> Bool {
  self.alert
}

///|
/// Reset the alert without resetting data.
pub fn SurvivalOnlineMonitor::clear_alert(
  self : SurvivalOnlineMonitor,
) -> SurvivalOnlineMonitor {
  {
    observations: self.observations,
    events: self.events,
    exposure: self.exposure,
    cumulative_hazard: self.cumulative_hazard,
    latest_time: self.latest_time,
    alert: false,
    monitor_row: true,
  }
}

///|
/// Export monitor state.
pub fn survival_online_monitor_to_csv(
  monitor : SurvivalOnlineMonitor,
) -> String {
  "metric,value\nobservations," +
  monitor.observations.to_string() +
  "\nevents," +
  monitor.events.to_string() +
  "\nexposure," +
  monitor.exposure.to_string() +
  "\ncumulative_hazard," +
  monitor.cumulative_hazard.to_string() +
  "\nlatest_time," +
  monitor.latest_time.to_string() +
  "\nalert," +
  monitor.alert.to_string() +
  "\n"
}

///|
/// Validate online state.
pub fn survival_online_monitor_is_valid(
  monitor : SurvivalOnlineMonitor,
) -> Bool {
  monitor.observations >= monitor.events &&
  monitor.events >= 0 &&
  monitor.exposure >= 0.0 &&
  monitor.cumulative_hazard >= 0.0 &&
  monitor.latest_time >= 0.0
}
