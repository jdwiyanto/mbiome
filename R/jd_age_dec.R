#' Convert age (numeric) to categorical variable (by decades)
#'
#' @param df A dataframe
#' @param agevar a character string of the age variable in the dataframe
#' @return A variable of age classified by decades
#' @export
jd_age_dec = function (df, agevar) {
  df[[agevar]] <- as.numeric(df[[agevar]])
  df$age_cat <- ifelse(df[[agevar]] < 11, "<11", ifelse(
                between(df[[agevar]], 11, 20), "11-20", ifelse(
                between(df[[agevar]], 21, 30), "21-30", ifelse(
                between(df[[agevar]], 31, 40), "31-40", ifelse(
                between(df[[agevar]], 41, 50), "41-50", ifelse(
                between(df[[agevar]], 51, 60), "51-60", ifelse(
                between(df[[agevar]], 61, 70), "61-70", ifelse(
                between(df[[agevar]], 71, 80), "71-80", ifelse(
                between(df[[agevar]], 81, 90), "81-90", ">90")))))))))

  df$age_cat = as.factor(df$age_cat)
  df$age_cat <- factor(df$age_cat, levels = c("<11", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", ">90"))
  return(df$age_cat)
}
