#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
univariate_product_inverse_monomial := MultivariatePowerSeries:-PowerSeriesObject:-univariate_product_inverse_monomial;
kernelopts(opaquemodules=true):


Try[testnoerror]("test uni inv 1", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='ps0');
Try[testnoerror]("test uni inv 2", PSO:-FromPolynomial(1 + x^4), 'assign'='ps1');
Try[testnoerror]("test uni inv 3", PSO:-Inverse(ps1)-1, 'assign'='ps_inv');

Try[testerror]("test uni inv 4", univariate_product_inverse_monomial(ps0, x));
Try[testerror]("test uni inv 5", univariate_product_inverse_monomial(ps1, x));
Try[testerror]("test uni inv 6", univariate_product_inverse_monomial(ps_inv, x^5));
Try[testnoerror]("test uni inv 7", univariate_product_inverse_monomial(ps_inv, 3*x), 'assign'='ps2');
Try("test uni inv 8", Truncate(ps2,0), 0);
Try("test uni inv 9", Truncate(ps2,8), x^7/3-x^3/3);
Try("test uni inv 10", HomogeneousPart(ps2,11), -1/3*x^11);

Try[testnoerror]("test uni inv 11", univariate_product_inverse_monomial(ps_inv, 2*x^4), 'assign'='ps3');
Try("test uni inv 12", Truncate(ps3,0), -1/2);
Try("test uni inv 13", Truncate(ps3,5),x^4/2-1/2);
Try("test uni inv 14", HomogeneousPart(ps3,8), -1/2*x^8);
Try("test uni inv 15", HomogeneousPart(ps3,12), 1/2*x^12);

Try[testnoerror]("test uni inv 16", univariate_product_inverse_monomial(ps1, 3));
Try[testerror]("test uni inv 17", univariate_product_inverse_monomial(ps1, 0), "invalid input: %1 must be different than zero");

#end test
