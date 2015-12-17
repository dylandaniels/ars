# This 'envelope' function works to output the value of u_k(x), which is the upper hull function when exponentiated and
# normalized. This function takes inputs of z(the intersection points), abscissae, the x value to be evaluated, hx(the
# log function value of the abscissae) and dhx(derivative values of log function).
envelope <- function(z, abscissae, x, hx, dhx) {
  #Here, 'index' determines which piecewise function will be used to evaluate u_k(x)
  #if the 'x' equals to last element of 'z', then the index will be length of 'z' minus
  #1,otherwise it will be the the sum of how many items in z are less than or equal to
  #'x'.
  if (x == z[length(z)]) {
    index = length(z) - 1
  } else {
    index <- sum(z <= x)
  }

  #plug in the value, compute u(x) and return the value.
  ux <- hx[index] + (x - abscissae[index]) * dhx[index]
  return(ux)
}

