
test_that("constrained Spikes works",{
  expect_error(Spikes(20, rep(1,2001), -0.5, steps =1000), "hyper must be >0. \n")
  expect_error(Spikes(20, rep(1,1501), 1.5, steps =1000), "The number of steps doesn't divide into the discretization of the intensity funtion.\n")
})

end.time = 20
x = rep(1,5001)
X <- cumsum(x)*(end.time/(length(x)-1))

test_that("constrained PDF works",{
  expect_error(PDF(1, last.spike=0.5, hyper=1.5, end.time=20, x, X, ISI.type='ChiSquared'), "Invalid ISI.type. See help documentation for allowable ISI.type.")
  expect_error(PDF(1, last.spike=0.5, hyper=1.5, end.time=20, x, X, ISI.type='Gamma2'), "Your choice of ISI.type requires two inputs for hyper.")
  expect_error(PDF(1, last.spike=0.5, hyper=c(1.5,4), end.time=20, x, X, ISI.type='Gamma'), "Your choice of ISI.type requires one inputs for hyper. For example hyper=1.")
})
