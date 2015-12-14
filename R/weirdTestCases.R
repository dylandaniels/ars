# Weird test cases

arcWithFlatPart <- Vectorize(function (x) {
  if (x >= 0 && x <= 1) {
    return(sqrt(1 - (x - 1)^2))
  } else {
    return(1)
  }
})

# should fail if unbounded
quadratic <- function (x) {x^2}
