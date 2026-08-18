///|
test "online moments match batch moments" {
  let values = [1.0, 2.0, 3.0, 4.0, 5.0]
  let online = OnlineMoments::new()
  for value in values {
    online.add(value)
  }
  let batch = numeric_summary(values)
  inspect(online.count, content="5")
  inspect(double_is_close(online.mean, batch.mean), content="true")
  inspect(double_is_close(online.variance(), batch.variance), content="true")
  inspect(online.minimum, content="1")
  inspect(online.maximum, content="5")
  inspect(online.summary().standard_deviation > 0.0, content="true")
}

///|
test "online moments merge is order independent" {
  let left = OnlineMoments::new()
  let right = OnlineMoments::new()
  left.add(1.0)
  left.add(2.0)
  right.add(3.0)
  right.add(4.0)
  left.merge(right)
  inspect(left.count, content="4")
  inspect(double_is_close(left.mean, 2.5), content="true")
  inspect(double_is_close(left.variance(), 1.6666666666666667), content="true")
}

///|
test "online event rate tracks exposure" {
  let rate = OnlineEventRate::new()
  rate.add(1.0, 1)
  rate.add(2.0, 0)
  rate.add(3.0, 1)
  inspect(rate.events, content="2")
  inspect(rate.censored, content="1")
  inspect(rate.exposure, content="6")
  inspect(rate.rate(), content="0.3333333333333333")
  inspect(rate.censoring_fraction(), content="0.3333333333333333")
}

///|
test "landmark estimates rebase follow-up" {
  let data = SurvivalData::new([1.0, 3.0, 5.0, 7.0], [1, 1, 0, 1])
  let landmark = landmark_estimate(data, 2.0, 3.0)
  inspect(landmark.eligible, content="3")
  inspect(landmark.subsequent_events, content="2")
  inspect(landmark.survival <= 1.0, content="true")
  inspect(landmark.to_csv().contains("landmark,eligible"), content="true")
  let table = discrete_hazard_table(data)
  inspect(table.length(), content="4")
  inspect(table[0].hazard > 0.0, content="true")
  inspect(table[3].survival_after >= 0.0, content="true")
}

///|
test "landmark empty cohort remains safe" {
  let estimate = landmark_estimate(SurvivalData::new([1.0], [1]), 2.0, 1.0)
  inspect(estimate.eligible, content="0")
  inspect(estimate.survival, content="1")
}
