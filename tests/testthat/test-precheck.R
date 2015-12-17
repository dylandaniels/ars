context('abscissae pre-check tests')

test_that('precheck() returns the correct output when the elements of the abscessae are not unique', {
  abs <- c(0,10,10) # mock abscessae points
  dhx <- c(-2,5,5)#mock dhx values
  leftbound <- -2#mock left bound
  rightbound <- 12#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'Elements of abscissae should be unique')
})

test_that('precheck() returns the correct output when the length of the abscessae is less than 2', {
  abs <- c(0)# mock abscessae points
  dhx <- c(-2,5)#mock dhx values
  leftbound <- -2#mock left bound
  rightbound <- 12#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'You must provide 2 or more abscissae.')
})

test_that('precheck() returns the correct output when the abscessaes are not within the boundaries', {
  abs <- c(0,10,20) # mock abscessae points
  dhx <- c(2,5,-5)#mock dhx values
  leftbound <- 2#mock left bound
  rightbound <- 30#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'Abscissae should be within boundaries.')
  abs <- c(0,10,20) # mock abscessae points
  dhx <- c(-2,5,5)#mock dhx values
  leftbound <- -2#mock left bound
  rightbound <- 15#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'Abscissae should be within boundaries.')
})

test_that('precheck() returns the correct output when the dhx values do not satisfy the condition in the paper',{
  abs <- c(0,10,20) # mock abscessae points
  dhx <- c(-2,5,5)#mock dhx values
  leftbound <- -Inf#mock left bound
  rightbound <- 30#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'Invalid abscissae or integral of function is divergent (cannot be normalized to a valid probability distribution)',fixed=TRUE)
  abs <- c(0,10,20) # mock abscessae points
  dhx <- c(2,5,5)#mock dhx values
  leftbound <- -2#mock left bound
  rightbound <- Inf#mock right bound
  expect_error(precheck(abs,dhx,leftbound,rightbound), 'Invalid abscissae or integral of function is divergent (cannot be normalized to a valid probability distribution)',fixed=TRUE)
})



