# Weird test cases

arcWithFlatPart <- Vectorize(function (x) {
  if (x >= 0 && x <= 2) {
    return(sqrt(4 - (x - 2)^2))
  } else {
    return(2)
  }
})

# should fail if unbounded
quadratic <- function (x) {x^2}
