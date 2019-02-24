# CC-BY-3.0 Brian Diggs
# https://stackoverflow.com/a/11054781/6646912
# https://creativecommons.org/licenses/by-sa/3.0/
library("scales")

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}
