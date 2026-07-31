///|
/// A step in the Kaplan-Meier survival curve.
pub(all) struct KMStep {
  time : Double
  n_at_risk : Int
  n_events : Int
  n_censored : Int
  survival : Double
  std_err : Double 
} derive(Eq, Debug)

///|
/// Result of a Kaplan-Meier estimation.
pub(all) struct KMResult {
  steps : Array[KMStep]
} derive(Eq, Debug)

///|
/// Estimates the Kaplan-Meier survival curve from the given data.
pub fn kaplan_meier(data : SurvivalData) -> KMResult {
  let steps = []
  let mut current_time = -1.0
  let mut n_at_risk = data.length()
  let mut current_events = 0
  let mut current_censored = 0
  let mut current_survival = 1.0
  let mut sum_var = 0.0
  
  if n_at_risk == 0 {
    return { steps: [] }
  }

  for i in 0..<data.length() {
    let obs = data.get(i)
    
    if current_time < 0.0 {
      current_time = obs.time
    }
    
    if obs.time > current_time {
      // Finalize previous time step
      if current_events > 0 || current_censored > 0 {
        let prob_survival = 1.0 - (current_events.to_double() / n_at_risk.to_double())
        current_survival = current_survival * prob_survival
        
        // Greenwood's formula for variance
        if current_events > 0 && n_at_risk > current_events {
            sum_var += current_events.to_double() / (n_at_risk.to_double() * (n_at_risk - current_events).to_double())
        }
        steps.push({
          time: current_time,
          n_at_risk,
          n_events: current_events,
          n_censored: current_censored,
          survival: current_survival,
          std_err: current_survival * double_sqrt(sum_var)
        })
      }
      
      n_at_risk = n_at_risk - current_events - current_censored
      current_time = obs.time
      current_events = 0
      current_censored = 0
    }
    
    if obs.event == 1 {
      current_events += 1
    } else {
      current_censored += 1
    }
  }
  
  // Finalize the last time step
  if current_events > 0 || current_censored > 0 {
    let prob_survival = 1.0 - (current_events.to_double() / n_at_risk.to_double())
    current_survival = current_survival * prob_survival
    
    if current_events > 0 && n_at_risk > current_events {
        sum_var += current_events.to_double() / (n_at_risk.to_double() * (n_at_risk - current_events).to_double())
    }

    steps.push({
      time: current_time,
      n_at_risk,
      n_events: current_events,
      n_censored: current_censored,
      survival: current_survival,
      std_err: current_survival * double_sqrt(sum_var)
    })
  }
  
  { steps, }
}

///|
/// Predict the survival probability at a given time `t`.
pub fn KMResult::predict(self : KMResult, t : Double) -> Double {
  let mut s = 1.0
  for i in 0..<self.steps.length() {
    if self.steps[i].time <= t {
      s = self.steps[i].survival
    } else {
      break
    }
  }
  s
}
