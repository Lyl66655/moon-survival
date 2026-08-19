///|
/// One immutable analysis-audit event.
pub(all) struct SurvivalAuditEvent {
  sequence : Int
  action : String
  subject : String
  detail : String
  checksum : String
  audit_event : Bool
} derive(Eq, Debug)

///|
/// An append-only audit trail.
pub(all) struct SurvivalAuditTrail {
  events : Array[SurvivalAuditEvent]
  trail_version : String
  audit_trail : Bool
} derive(Eq, Debug)

///|
/// Start an audit trail.
pub fn survival_audit_trail() -> SurvivalAuditTrail {
  { events: [], trail_version: "1", audit_trail: true }
}

///|
/// Append an audit event.
pub fn SurvivalAuditTrail::append(
  self : SurvivalAuditTrail,
  action : String,
  subject : String,
  detail : String,
) -> SurvivalAuditTrail {
  let sequence = self.events.length() + 1
  let content = action +
    "|" +
    subject +
    "|" +
    detail +
    "|" +
    sequence.to_string()
  let events : Array[SurvivalAuditEvent] = []
  for event in self.events {
    events.push(event)
  }
  events.push({
    sequence,
    action,
    subject,
    detail,
    checksum: survival_export_checksum(content),
    audit_event: true,
  })
  { events, trail_version: self.trail_version, audit_trail: true }
}

///|
/// Append a data-ingestion event.
pub fn SurvivalAuditTrail::record_ingest(
  self : SurvivalAuditTrail,
  data : SurvivalData,
) -> SurvivalAuditTrail {
  self.append(
    "ingest",
    "SurvivalData",
    "rows=" +
    data.length().to_string() +
    ",events=" +
    data.event_count().to_string(),
  )
}

///|
/// Append a model-fit event.
pub fn SurvivalAuditTrail::record_fit(
  self : SurvivalAuditTrail,
  model_name : String,
  score : Double,
) -> SurvivalAuditTrail {
  self.append("fit", model_name, "score=" + score.to_string())
}

///|
/// Append an export event.
pub fn SurvivalAuditTrail::record_export(
  self : SurvivalAuditTrail,
  artifact : String,
  checksum : String,
) -> SurvivalAuditTrail {
  self.append("export", artifact, "checksum=" + checksum)
}

///|
/// Return the last audit event.
pub fn SurvivalAuditTrail::last(
  self : SurvivalAuditTrail,
) -> SurvivalAuditEvent? {
  if self.events.length() == 0 {
    None
  } else {
    Some(self.events[self.events.length() - 1])
  }
}

///|
/// Verify sequence and checksum fields.
pub fn SurvivalAuditTrail::is_valid(self : SurvivalAuditTrail) -> Bool {
  if self.trail_version == "" {
    return false
  }
  for i in 0..<self.events.length() {
    if self.events[i].sequence != i + 1 ||
      self.events[i].action == "" ||
      self.events[i].checksum == "" {
      return false
    }
  }
  true
}

///|
/// Export an audit trail.
pub fn SurvivalAuditTrail::to_csv(self : SurvivalAuditTrail) -> String {
  let builder = StringBuilder::new()
  builder.write_string("sequence,action,subject,detail,checksum\n")
  for event in self.events {
    builder.write_string(event.sequence.to_string())
    builder.write_char(',')
    builder.write_string(event.action)
    builder.write_char(',')
    builder.write_string(event.subject)
    builder.write_char(',')
    builder.write_string(event.detail)
    builder.write_char(',')
    builder.write_string(event.checksum)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Compute the trail checksum.
pub fn survival_audit_checksum(trail : SurvivalAuditTrail) -> String {
  survival_export_checksum(trail.to_csv())
}
