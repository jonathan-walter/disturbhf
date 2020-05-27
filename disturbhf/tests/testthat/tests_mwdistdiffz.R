context("mwdistdiffz")

test_that("test error checking",{
  d1<-data.frame(tt=1:100,yy=rnorm(100))
  d2<-data.frame(tt=1:100,yy=rnorm(100))
  expect_error(mwdistdiffz(testy=as.matrix(d1), refy=d2, wwidth=10),
               "testy must be a data frame containing the columns 'tt' and 'yy'")
  expect_error(mwdistdiffz(testy=d1, refy=as.matrix(d2), wwidth=10),
               "refy must be a data frame containing the columns 'tt' and 'yy'")
  expect_error(mwdistdiffz(testy=d1, refy=d2, wwidth=10, stat="median"),
               "available values of 'stat' are mean and sum")
})

test_that("test warnings",{
  d1<-data.frame(tt=1:30,yy=rnorm(30))
  d2<-data.frame(tt=1:50,yy=rnorm(50))
  expect_warning(mwdistdiffz(testy=d1, refy=d2, wwidth=10),
                 "the test period is very short, results may not be robust")
})