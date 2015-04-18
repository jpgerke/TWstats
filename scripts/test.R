rm(list=ls())
require(testthat)
source("tw.r")
load("../test/TWvals.RData")
load("../test//Pvals.table.Rdata")

eigs = scan("../test/eigs.txt", what=numeric())


###### testing ######
mytest  = TWcalc(eigs, 60)

expect_equal(mytest[[1]], TWvals, tolerance=0.00001)
expect_equal(mytest[[2]], Pvals.table, tolerance=0.00001)
expect_equal(mytest[[3]], Pvals.dist, tolerance=0.00001)
