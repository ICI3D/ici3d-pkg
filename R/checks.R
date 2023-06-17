
# Checkers for use in various ICI3D functions

check_character <- checker(is.character(x), "`%s` must be a 'character'.")
check_scalar <- checker(length(x) == 1, "`%s` must be a scalar (length == 1).")
check_nonemptychar <- checker(all(nchar(x) > 0), "`%s` must be non-empty.")
check_numeric <- checker(is.numeric(x), "`%s` must be numeric.")
check_integer <- checker(is.integer(x) | all(trunc(x) == x), "`%s` must be an integer.")
check_positive <- checker(all(x > 0), "`%s` must be positive.")
check_nonnegative <- checker(all(x >= 0), "`%s` must be nonnegative.")
check_probability <- checker(data.table::between(x, 0, 1), "`%s` must be a probability.")
check_lt <- checker_against(x < ref, "`%s` must be < %s.")
check_lte <- checker_against(x <= ref, "`%s` must be <= %s.")
check_series <- checker(all(diff(x) == 1), "`%s` must be a seq(..., by=1).")
check_among <- checker_against(all(x %in% ref), "all `%s` must be among %s.")
check_contains <- checker_against(all(ref %in% x), "`%s` must contain all %s.")
check_equal <- checker_against(x == ref, "`%s` must be equal to %s.")
