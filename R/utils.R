#################################################
# Helper functions
#################################################

# Returns midpoints of even intervals between from and to
# E.g.
# midpoint(2,10,4)
# 3,5,7,9
.midpoints <- function(from = 1, to = 1, length.out = 1){
  x <- seq(from, to, length.out = length.out + 1)
  (tail(x, -1)+head(x, -1))/2
}

# Common map
# Returns a mapping as small as possible that both dfg's can use
# f and g are the two mappings
# Find a mapping, h, to the smallest set such that if f(i) != f(j)
# then h(i) != h(j) and similarly for g.
.commonMap <- function(potMap1, potMap2){
  enc <- potMap1+length(potMap1)*potMap2 # encode pairs
  match(match(enc,enc), unique(match(enc, enc))) # unique pairs map to same element
}