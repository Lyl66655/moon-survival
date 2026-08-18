///|
test "csv parser handles quoted fields and empty cells" {
  let table = csv_parse(
    "name,value,note\nalpha,1.5,\"hello, world\"\nbeta,,ok\n",
  ).unwrap()
  inspect(table.headers.length(), content="3")
  inspect(table.rows.length(), content="2")
  inspect(table.rows[0][2], content="hello, world")
  inspect(table.rows[1][1], content="")
  match table.column("value") {
    Some(value) => inspect(value, content="1")
    None => fail("value column was not found")
  }
  inspect(csv_escape("hello,world"), content="\"hello,world\"")
}

///|
test "csv survival and cox loaders validate numeric columns" {
  let text = "time,event,age,score\n1,1,45,0.2\n2,0,51,0.4\n3,1,60,0.8\n"
  let survival = survival_data_from_csv(text, "time", "event").unwrap()
  inspect(survival.length(), content="3")
  let cox = cox_data_from_csv(text, "time", "event", ["age", "score"]).unwrap()
  inspect(cox.n_covariates, content="2")
  inspect(cox.length(), content="3")
  match
    survival_data_from_csv("time,event\nnot-a-number,1\n", "time", "event") {
    Err(InvalidNumber(_)) => ()
    _ => fail("invalid numeric value was accepted")
  }
}

///|
test "csv round trip and malformed input boundaries" {
  let original = SurvivalData::new([1.0, 2.0], [1, 0])
  let encoded = survival_data_to_csv(original)
  let decoded = survival_data_from_csv(encoded, "time", "event").unwrap()
  inspect(decoded.length(), content="2")
  inspect(decoded.get(1).event, content="0")
  match csv_parse("a,b\n1\n") {
    Err(RowWidthMismatch) => ()
    _ => fail("short row was accepted")
  }
  match csv_parse("a,b\n\"unterminated,1\n") {
    Err(UnterminatedQuote) => ()
    _ => fail("unterminated quote was accepted")
  }
}

///|
test "time varying data checks intervals and snapshots" {
  let data = TimeVaryingData::new(
    [0.0, 2.0, 0.0, 1.0],
    [2.0, 5.0, 1.0, 4.0],
    [0, 1, 0, 1],
    [[0.0], [1.0], [2.0], [3.0]],
    [0, 0, 1, 1],
  )
  inspect(data.length(), content="4")
  inspect(data.n_subjects, content="2")
  inspect(data.is_valid_subject(0), content="true")
  inspect(data.is_valid_subject(1), content="true")
  inspect(time_varying_risk_table(data).length(), content="2")
  inspect(data.subject_intervals(0).length(), content="2")
  inspect(landmark_snapshot(data, 0.5).length(), content="2")
  inspect(data.to_survival_data().length(), content="2")
}
