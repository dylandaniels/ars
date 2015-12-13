# Weird test cases

somePartFlat <- function (x) {
  if (x >=0 && x <= 10) {
    return(x^2)
  } else {
    return(100)
  }
}

# should fail if unbounded
quadratic <- function (x) {x^2}
