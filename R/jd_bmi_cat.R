#' Convert BMI (numeric) to category
#'
#' @param df A dataframe, with factor containing BMI variable named as 'bmi'
#' @return A variable of BMI classified by BMI categories based on MOH Malaysia
#' @export
jd_bmi_cat = function(df) {
  df$bmi = as.numeric(df$bmi)
  df$bmi_cat <- ifelse(df$bmi < 18.5, 'underweight', ifelse(between(df$bmi, 18.5, 24.99), 'normal', ifelse(between(df$bmi, 25, 29.99), 'overweight', 'obese')))
  df$bmi_cat = as.factor(df$bmi_cat)
  df$bmi_cat = factor(df$bmi_cat, levels = c('underweight', 'normal', 'overweight', 'obese'))
  return(df$bmi_cat)
}
