///|
/// A single observation for Cox Proportional Hazards model with covariates.
pub(all) struct CoxObservation {
  time : Double
  event : Int
  covariates : Array[Double]
} derive(Eq, Debug)

///|
/// CoxData represents a collection of observations with covariates.
pub(all) struct CoxData {
  observations : Array[CoxObservation]
  n_covariates : Int
} derive(Eq, Debug)

///|
pub fn CoxData::new(times : Array[Double], events : Array[Int], covs : Array[Array[Double]]) -> CoxData {
  if times.length() != events.length() || times.length() != covs.length() {
    panic()
  }
  let n_cov = if covs.length() > 0 { covs[0].length() } else { 0 }
  let observations = []
  for i in 0..<times.length() {
    if covs[i].length() != n_cov {
      panic()
    }
    observations.push({ time: times[i], event: events[i], covariates: covs[i] })
  }
  
  observations.sort_by(fn (a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      b.event.compare(a.event)
    }
  })
  { observations, n_covariates: n_cov }
}

///|
/// Calculates the partial likelihood of the Cox model given coefficients `beta`.
pub fn cox_partial_likelihood(data : CoxData, beta : Array[Double]) -> Double {
  if beta.length() != data.n_covariates {
    panic()
  }
  
  let mut lik = 1.0
  
  for i in 0..<data.observations.length() {
    let obs = data.observations[i]
    if obs.event == 1 {
      let mut bx = 0.0
      for j in 0..<beta.length() {
        bx += beta[j] * obs.covariates[j]
      }
      
      let mut risk_sum = 0.0
      for k in i..<data.observations.length() {
        let risk_obs = data.observations[k]
        let mut bk = 0.0
        for j in 0..<beta.length() {
          bk += beta[j] * risk_obs.covariates[j]
        }
        risk_sum += double_exp(bk)
      }
      
      lik = lik * (double_exp(bx) / risk_sum)
    }
  }
  
  lik
}
