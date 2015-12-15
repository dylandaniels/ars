context('test the Mainline function for abnormal input functions')

test_that('test the Mainline function with piecewise input functions but continuous and differentiable everywhere', {
  set.seed(123)
  myfun <- Vectorize(function (x) {
    if (x >= 0 && x <= 1) {
      return(sqrt(1 - (x - 1)^2))
    } else {
      return(1)
    }
  })
  expect_equal(Mainline(10,myfun,initialPoints =c(1,2),leftbound=0,rightbound=5)[1:3],
               c(1.4378876006230712, 2.0448846090584993, 4.7023364214692265))
})

test_that('test the Mainline function with piecewise input functions but some points not differentiable', {
  set.seed(123)
  myfun<- Vectorize(function (x) {
    if (x >= 0 && x <= 1) {
      return(x^2)
    } else {
      return(1)
    }
  })
  expect_equal(Mainline(10,myfun,initialPoints = c(1,2),leftbound=0,rightbound=5)[1:3],
               c(1.4378876006230712, 2.0448846090584993, 4.7023364214692265))
})

test_that('test the Mainline functions when input functions are not continuous somewhere', {
  set.seed(123)
  myfun <- Vectorize(function(x){
    if (x>=0 && x<=3){
      return(x)
    }else if (x>3 && x<=5){
      return(5)
    }else if(x>5){
      return(-x+7)
    }
  })
  expect_equal(Mainline(10,myfun,initialPoints = c(0.5,1),leftbound = 0,rightbound = 6)[1:3],
            c(3.45943279345439825, 5.61826353815976187, 2.70759249039256655)
            )
})

test_that('test the Mainline functions when input functions are not log-concave', {

  myfun <- Vectorize(function(x){
   return((x-5)^2)
  })
  expect_error(Mainline(100,myfun,initialPoints = c(4,6),leftbound = 0,rightbound = 10),
               "The sampling function is not log-concave."
  )
})



