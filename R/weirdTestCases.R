# Weird test cases

arcWithFlatPart <- Vectorize(function (x) {
  if (x >= 0 && x <= 1) {
    return(sqrt(1 - (x - 1)^2))
  } else {
    return(1)
  }
})

twoArcsWithFlatPart <- Vectorize(function (x) {
  if (x >= 0 && x <= 1) {
    return(sqrt(1 - (x - 1)^2))
  } else if (x >= 4 && x <= 5) {
    return(sqrt(1 - (x - 4)^2))
  } else {
    return(1)
  }
})
>>>>>>> e128de73db98b273d85bc7b26f12f676955ce2c5

# should fail if unbounded
quadratic <- function (x) {x^2}
