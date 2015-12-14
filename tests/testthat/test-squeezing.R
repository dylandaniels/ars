

test_that('squeezing() returns the correct output when given value is not in the scope of hx', {
  abs=c(0,1,2) # mock abscessae points
  x=-1 # mock x value
  expect_equal(squeezing(hx,abs, x), -Inf)
})

test_that("squeeing() returns the correct output when given values is in the scope of hx", {
  abs=c(1,2,4)
  hx=c(1,4,8)
  x1=1.3
  x2=2.7
  x3=2
  expect_equal(squeezing(hx,abs,x1),1.9)
  expect_equal(squeezing(hx,abs,x2),5.4)
  expect_equal(squeezing(hx,abs,x3),4)
  
})
