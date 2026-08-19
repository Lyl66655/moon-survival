///|
/// A prioritized follow-up queue item.
pub(all) struct SurvivalQueueItem {
  subject : Int
  risk : Double
  priority : Int
  due_time : Double
  reason : String
  queue_item : Bool
} derive(Eq, Debug)

///|
/// A capacity-constrained queue allocation.
pub(all) struct SurvivalQueueAllocation {
  items : Array[SurvivalQueueItem]
  capacity : Int
  allocated : Int
  unmet : Int
  allocation_row : Bool
} derive(Eq, Debug)

///|
/// Rank subjects by predicted risk and deadline.
pub fn survival_queue_build(
  subjects : Array[Int],
  risks : Array[Double],
  due_times : Array[Double],
  threshold : Double,
) -> Array[SurvivalQueueItem] {
  let result : Array[SurvivalQueueItem] = []
  let count = if subjects.length() < risks.length() {
    subjects.length()
  } else {
    risks.length()
  }
  for i in 0..<count {
    let due = if i < due_times.length() { due_times[i] } else { 0.0 }
    if risks[i] >= threshold {
      result.push({
        subject: subjects[i],
        risk: clamp_probability(risks[i]),
        priority: 0,
        due_time: due,
        reason: "risk_threshold",
        queue_item: true,
      })
    }
  }
  result.sort_by(fn(a, b) {
    if a.risk > b.risk {
      -1
    } else if a.risk < b.risk {
      1
    } else if a.due_time < b.due_time {
      -1
    } else {
      a.subject.compare(b.subject)
    }
  })
  for i in 0..<result.length() {
    result[i] = {
      subject: result[i].subject,
      risk: result[i].risk,
      priority: i + 1,
      due_time: result[i].due_time,
      reason: result[i].reason,
      queue_item: true,
    }
  }
  result
}

///|
/// Allocate a queue to available capacity.
pub fn survival_queue_allocate(
  items : Array[SurvivalQueueItem],
  capacity : Int,
) -> SurvivalQueueAllocation {
  let bounded = if capacity < 0 { 0 } else { capacity }
  let allocated = if items.length() < bounded {
    items.length()
  } else {
    bounded
  }
  let selected : Array[SurvivalQueueItem] = []
  for i in 0..<allocated {
    selected.push(items[i])
  }
  {
    items: selected,
    capacity: bounded,
    allocated,
    unmet: items.length() - allocated,
    allocation_row: true,
  }
}

///|
/// Recompute priorities after a queue update.
pub fn survival_queue_reprioritize(
  items : Array[SurvivalQueueItem],
  urgency_weight : Double,
) -> Array[SurvivalQueueItem] {
  let result : Array[SurvivalQueueItem] = []
  for item in items {
    let score = item.risk +
      clamp_positive(urgency_weight) / clamp_positive(1.0 + item.due_time)
    result.push({
      subject: item.subject,
      risk: score,
      priority: item.priority,
      due_time: item.due_time,
      reason: "reprioritized",
      queue_item: true,
    })
  }
  result.sort_by(fn(a, b) {
    if a.risk > b.risk {
      -1
    } else if a.risk < b.risk {
      1
    } else {
      a.subject.compare(b.subject)
    }
  })
  result
}

///|
/// Return queue load by due-time windows.
pub fn survival_queue_window_counts(
  items : Array[SurvivalQueueItem],
  width : Double,
  windows : Int,
) -> Array[Int] {
  let result = Array::make(if windows < 0 { 0 } else { windows }, 0)
  if width <= 0.0 {
    return result
  }
  for item in items {
    let index = (item.due_time / width).to_int()
    if index >= 0 && index < result.length() {
      result[index] += 1
    }
  }
  result
}

///|
/// Compute the mean risk of a queue.
pub fn survival_queue_mean_risk(items : Array[SurvivalQueueItem]) -> Double {
  let values : Array[Double] = []
  for item in items {
    values.push(item.risk)
  }
  survival_ext_mean(values)
}

///|
/// Export queue items.
pub fn survival_queue_to_csv(items : Array[SurvivalQueueItem]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,risk,priority,due_time,reason\n")
  for item in items {
    builder.write_string(item.subject.to_string())
    builder.write_char(',')
    builder.write_string(item.risk.to_string())
    builder.write_char(',')
    builder.write_string(item.priority.to_string())
    builder.write_char(',')
    builder.write_string(item.due_time.to_string())
    builder.write_char(',')
    builder.write_string(item.reason)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate queue allocation accounting.
pub fn survival_queue_is_valid(allocation : SurvivalQueueAllocation) -> Bool {
  allocation.capacity >= 0 &&
  allocation.allocated == allocation.items.length() &&
  allocation.unmet >= 0 &&
  allocation.allocated <= allocation.capacity
}
