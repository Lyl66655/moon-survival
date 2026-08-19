///|
/// A typed column schema for survival datasets.
pub(all) enum SurvivalColumnType {
  TimeColumn
  EventColumn
  GroupColumn
  WeightColumn
  CovariateColumn
} derive(Eq, Debug)

///|
/// A column declaration.
pub(all) struct SurvivalColumnSchema {
  name : String
  column_type : SurvivalColumnType
  required : Bool
  minimum : Double
  maximum : Double
  schema_column : Bool
} derive(Eq, Debug)

///|
/// Dataset schema with a version and provenance tag.
pub(all) struct SurvivalSchema {
  name : String
  version : String
  columns : Array[SurvivalColumnSchema]
  allow_missing : Bool
  schema_row : Bool
} derive(Eq, Debug)

///|
/// A schema-validation result.
pub(all) struct SurvivalSchemaViolation {
  column : String
  row : Int
  code : String
  message : String
  schema_violation : Bool
} derive(Eq, Debug)

///|
/// Create the core schema.
pub fn survival_core_schema() -> SurvivalSchema {
  {
    name: "survival-data",
    version: "1",
    columns: [
      {
        name: "time",
        column_type: TimeColumn,
        required: true,
        minimum: 0.0,
        maximum: 1.0e300,
        schema_column: true,
      },
      {
        name: "event",
        column_type: EventColumn,
        required: true,
        minimum: 0.0,
        maximum: 1.0,
        schema_column: true,
      },
    ],
    allow_missing: false,
    schema_row: true,
  }
}

///|
/// Add a covariate declaration.
pub fn survival_schema_add_covariate(
  schema : SurvivalSchema,
  name : String,
  minimum : Double,
  maximum : Double,
) -> SurvivalSchema {
  let columns : Array[SurvivalColumnSchema] = []
  for column in schema.columns {
    columns.push(column)
  }
  columns.push({
    name,
    column_type: CovariateColumn,
    required: false,
    minimum,
    maximum,
    schema_column: true,
  })
  {
    name: schema.name,
    version: schema.version,
    columns,
    allow_missing: schema.allow_missing,
    schema_row: true,
  }
}

///|
/// Return a column-type label.
pub fn survival_column_type_name(column_type : SurvivalColumnType) -> String {
  match column_type {
    TimeColumn => "time"
    EventColumn => "event"
    GroupColumn => "group"
    WeightColumn => "weight"
    CovariateColumn => "covariate"
  }
}

///|
/// Validate core SurvivalData against a schema.
pub fn survival_schema_validate(
  schema : SurvivalSchema,
  data : SurvivalData,
) -> Array[SurvivalSchemaViolation] {
  let result : Array[SurvivalSchemaViolation] = []
  if schema.name == "" || schema.columns.length() == 0 {
    result.push({
      column: "schema",
      row: -1,
      code: "INVALID_SCHEMA",
      message: "schema has no name or columns",
      schema_violation: true,
    })
  }
  for i in 0..<data.length() {
    let observation = data.get(i)
    if observation.time < 0.0 ||
      observation.time.is_nan() ||
      observation.time.is_inf() {
      result.push({
        column: "time",
        row: i,
        code: "INVALID_TIME",
        message: "time is invalid",
        schema_violation: true,
      })
    }
    if observation.event != 0 && observation.event != 1 {
      result.push({
        column: "event",
        row: i,
        code: "INVALID_EVENT",
        message: "event is not binary",
        schema_violation: true,
      })
    }
  }
  result
}

///|
/// Validate a numeric covariate matrix against appended covariate columns.
pub fn survival_schema_validate_covariates(
  schema : SurvivalSchema,
  values : Array[Array[Double]],
) -> Array[SurvivalSchemaViolation] {
  let result : Array[SurvivalSchemaViolation] = []
  let covariates : Array[SurvivalColumnSchema] = []
  for column in schema.columns {
    if column.column_type == CovariateColumn {
      covariates.push(column)
    }
  }
  for row_index in 0..<values.length() {
    if values[row_index].length() != covariates.length() {
      result.push({
        column: "covariates",
        row: row_index,
        code: "WIDTH_MISMATCH",
        message: "covariate width does not match schema",
        schema_violation: true,
      })
      continue
    }
    for column in 0..<covariates.length() {
      let value = values[row_index][column]
      if value.is_nan() || value.is_inf() {
        if !schema.allow_missing {
          result.push({
            column: covariates[column].name,
            row: row_index,
            code: "MISSING_VALUE",
            message: "missing covariate",
            schema_violation: true,
          })
        }
      } else if value < covariates[column].minimum ||
        value > covariates[column].maximum {
        result.push({
          column: covariates[column].name,
          row: row_index,
          code: "OUT_OF_RANGE",
          message: "covariate outside declared range",
          schema_violation: true,
        })
      }
    }
  }
  result
}

///|
/// Validate all schema constraints and return a boolean.
pub fn survival_schema_is_valid(
  schema : SurvivalSchema,
  data : SurvivalData,
) -> Bool {
  survival_schema_validate(schema, data).length() == 0
}

///|
/// Count violations by code.
pub fn survival_schema_violation_counts(
  violations : Array[SurvivalSchemaViolation],
) -> Array[String] {
  let result : Array[String] = []
  for violation in violations {
    if !result.contains(violation.code) {
      result.push(violation.code)
    }
  }
  result
}

///|
/// Export schema metadata.
pub fn SurvivalSchema::to_csv(self : SurvivalSchema) -> String {
  let builder = StringBuilder::new()
  builder.write_string("name,version,column,type,required,minimum,maximum\n")
  for column in self.columns {
    builder.write_string(self.name)
    builder.write_char(',')
    builder.write_string(self.version)
    builder.write_char(',')
    builder.write_string(column.name)
    builder.write_char(',')
    builder.write_string(survival_column_type_name(column.column_type))
    builder.write_char(',')
    builder.write_string(column.required.to_string())
    builder.write_char(',')
    builder.write_string(column.minimum.to_string())
    builder.write_char(',')
    builder.write_string(column.maximum.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export schema violations.
pub fn survival_schema_violations_to_csv(
  violations : Array[SurvivalSchemaViolation],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("column,row,code,message\n")
  for violation in violations {
    builder.write_string(violation.column)
    builder.write_char(',')
    builder.write_string(violation.row.to_string())
    builder.write_char(',')
    builder.write_string(violation.code)
    builder.write_char(',')
    builder.write_string(violation.message)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check that schema names and ranges are meaningful.
pub fn survival_schema_definition_is_valid(schema : SurvivalSchema) -> Bool {
  if schema.name == "" || schema.version == "" || schema.columns.length() == 0 {
    return false
  }
  for column in schema.columns {
    if column.name == "" || column.minimum > column.maximum {
      return false
    }
  }
  true
}
