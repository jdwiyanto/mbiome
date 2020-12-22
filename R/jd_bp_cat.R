#' Convert BP systolic and diastolic (numeric) to category
#'
#' @param df A dataframe, with factor containing BMI variable named as 'bp_dia' and 'bp_sys'
#' @return A variable of BP classified based on the American Heart Association (https://www.heart.org/en/health-topics/high-blood-pressure/understanding-blood-pressure-readings)
#' @export
jd_bp_cat = function(df) {
  df$bp_dia = df$bp_dia %>% as.numeric
  df$bp_sys = df$bp_sys %>% as.numeric
  bp = ifelse((df$bp_sys < 120 & df$bp_dia < 80), 'normal', ifelse(
    (between(df$bp_sys, 120, 129) & df$bp_dia < 80), 'elevated', ifelse(
      (between(df$bp_sys, 130, 139) | between(df$bp_dia, 80, 89)), 'stage1' ,ifelse(
        (between(df$bp_sys, 140, 180) | between(df$bp_dia, 90, 120)), 'stage2', 'crisis'))))
  bp = as.factor(bp)
  bp = factor(bp, levels = c('normal', 'elevated', 'stage1', 'stage2', 'crisis'))
  return(bp)
}
