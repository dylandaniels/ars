#***********document*********************************************************
#This function is used to test whether the abscissae satisfies preliminary requirements.
#we import dhx,abscissae leftbound and rightbound and returns different imformation if
#the abscissae doesn't satisfies different requirements.And it will return no message if
#the abscissae is OK.
#*****************************************************************************
precheck <- function(abscissae, dhx, leftbound, rightbound){
  if (length(unique(abscissae)) != length(abscissae)) {
    stop('Elements of abscissae should be unique')
  }

  if (length(abscissae) < 2) {
    stop('You must provide 2 or more abscissae.')
  }

  if (leftbound > abscissae[1] || rightbound < abscissae[length(abscissae)]) {
    stop('Abscissae should be within boundaries.')
  }

  if ((is.infinite(leftbound) && dhx[1] <= 0) ||  (is.infinite(rightbound) && dhx[length(dhx)] >= 0)) {
    # make this more descript later.
    stop('Invalid abscissae or integral of function is divergent (cannot be normalized to a valid probability distribution)')
  }
}
