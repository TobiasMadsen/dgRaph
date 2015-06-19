# Helper function: returns midpoints of even intervals between from and to
# E.g.
# midpoint(2,10,4)
# 3,5,7,9
.midpoints <- function(from = 1, to = 1, length.out = 1){
  x <- seq(from, to, length.out = length.out + 1)
  (tail(x, -1)+head(x, -1))/2
}