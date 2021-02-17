#' Convert BMI (numeric) to category
#'
#' @param df A dataframe
#' @param bmivar a character string of the BMI variable in the dataframe
#' @return A variable of BMI classified by BMI categories based on MOH Malaysia
#' @export
jd_bmi_cat = function (df, bmivar) {
  df[[bmivar]] = as.numeric(df[[bmivar]])
  df$bmi_cat <- ifelse(df[[bmivar]] < 18.5, "underweight", ifelse(
                between(df[[bmivar]], 18.5, 24.99), "normal", ifelse(
                between(df[[bmivar]], 25, 29.99), "overweight", "obese")))
  df$bmi_cat = as.factor(df$bmi_cat)
  df$bmi_cat = factor(df$bmi_cat, levels = c("underweight", "normal", "overweight", "obese"))
  return(df$bmi_cat)
}
