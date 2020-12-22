#' Convert age (numeric) to categorical variable (by decades)
#'
#' @param df A dataframe, with factor containing age variable named as 'age'
#' @return A variable of age classified by decades
#' @export
jd_age_dec = function(df) {
df$age <- as.numeric(df$age)
df$age_cat <- ifelse(df$age < 11, '<11',
                     ifelse(between(df$age, 11, 20), '11-20',
                            ifelse(between(df$age, 21, 30), '21-30',
                                   ifelse(between(df$age, 31, 40), '31-40',
                                          ifelse(between(df$age, 41, 50), '41-50',
                                                 ifelse(between(df$age, 51, 60), '51-60',
                                                        ifelse(between(df$age, 61, 70), '61-70',
                                                               ifelse(between(df$age, 71, 80), '71-80',
                                                                      ifelse(between(df$age, 81, 90), '81-90', '>90')))))))))
df$age_cat = as.factor(df$age_cat)
df$age_cat <- factor(df$age_cat, levels = c('<11',
                                            '11-20',
                                            '21-30',
                                            '31-40',
                                            '41-50',
                                            '51-60',
                                            '61-70',
                                            '71-80',
                                            '81-90',
                                            '>90'))
return(df$age_cat)
}
