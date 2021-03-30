#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));
# kernelopts(opaquemodules=false):
# UPoPSObject := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
# kernelopts(opaquemodules=true):

Try[testnoerror]("test 1", HenselFactorize(UnivariatePolynomialOverPowerSeries(1)), 'assign'='up1');
Try("test 2", numelems(up1), 0);

Try[testerror]("test 4", HenselFactorize(UnivariatePolynomialOverPowerSeries(0)), "invalid input: "
               "expected univariate polynomial over power series with invertible leading "
               "coefficient, but the leading coefficient has constant part 0");

Try[testnoerror]("test 7", HenselFactorize(UnivariatePolynomialOverPowerSeries(RootOf(Z^5-1))), 'assign'='up1');
Try("test 8", numelems(up1), 1);
Try("test 9", Truncate(up1[1], 10), RootOf(Z^5-1));

p := z^5 - 3*z^4 - 23*z^3 + 51*z^2 + 94*z - 120;
Try[testnoerror]("test 10", HenselFactorize(UnivariatePolynomialOverPowerSeries(p, z)), 'assign'='fact');
Try("test 11", numelems(fact), 5);
Try("test 12", expand(mul(Truncate(fact[i]), i = 1 .. 5)), p);

Try[testnoerror]("test 13", HenselFactorize(UnivariatePolynomialOverPowerSeries((z + sqrt(3))*(z^2 + 1), z)), 'assign'='fact');
Try("test 14", numelems(fact), 3);
Try("test 15", expand(mul(Truncate(fact[i]), i = 1 .. 3)), RootOf(_Z^2 - 3, index = 1)*z^2 + z^3 + RootOf(_Z^2 - 3, index = 1) + z);

Try[testnoerror]("test 16", HenselFactorize(UnivariatePolynomialOverPowerSeries((z - 1)*(z - 2)*(z - 3) + x*(z^2 + z), z)), 'assign'='fact');
Try("test 17", numelems(fact), 3);
Try[testnoerror]("test 18", UpdatePrecision(fact[1], 10));
Try[testnoerror]("test 19", UpdatePrecision(fact[2], 10));
Try[testnoerror]("test 20", UpdatePrecision(fact[3], 10));

Try[testnoerror]("test 21", HenselFactorize(UnivariatePolynomialOverPowerSeries((z - 1)*(z - 2)*(z - 3)*(z - 4) + x*(z^3 + z), z)), 'assign'='fact');
Try("test 22", numelems(fact), 4);
Try[testnoerror]("test 23", UpdatePrecision(fact[1], 10));
Try[testnoerror]("test 24", UpdatePrecision(fact[2], 10));
Try[testnoerror]("test 25", UpdatePrecision(fact[3], 10));
Try[testnoerror]("test 26", UpdatePrecision(fact[3], 10));

test_hf := proc(f, d) local up, res; up := HenselFactorize(f); res := Multiply(seq(up)); res := Truncate(Subtract(res, f), d); if res = 0 then return true; else return res; end if; end proc;

Try[testnoerror]("test 27", UnivariatePolynomialOverPowerSeries(x^2 + 2*x + 1, x), 'assign'='f');
Try("test 28", test_hf(f, 20), true);

Try[testnoerror]("test 29", UnivariatePolynomialOverPowerSeries(-x*z + z^2 - 1, z), 'assign'='f');
Try("test 30", test_hf(f, 20), true);

Try[testnoerror]("test 31", UnivariatePolynomialOverPowerSeries(z^2 + (2 - x)*z - x, z), 'assign'='f');
Try("test 32", test_hf(f, 20), true);

Try[testnoerror]("test 33", UnivariatePolynomialOverPowerSeries(2 + (-1 - x)*z + (-2 + x)*z^2 + z^3, z), 'assign'='f');
Try("test 34", test_hf(f, 20), true);

Try[testnoerror]("test 35", UnivariatePolynomialOverPowerSeries(3*x^2*z^2 + z^3 + x*z - 2*z^2 + x - z + 2 , z), 'assign'='f');
Try("test 36", test_hf(f, 20), true);

Try[testnoerror]("test 37", UnivariatePolynomialOverPowerSeries(z^3 + z^2 + y, z), 'assign'='f');
Try("test 38", test_hf(f, 20), true);

Try[testnoerror]("test 39", UnivariatePolynomialOverPowerSeries(8 + x + (y^2 + 12)*z + (y + 6)*z^2 + z^3, z), 'assign'='f');
Try("test 40", test_hf(f, 20), true);

Try[testnoerror]("test 41", UnivariatePolynomialOverPowerSeries(y^2 + x^2 + (y + 1)*z^2 + z^3, z), 'assign'='f');
Try("test 42", test_hf(f, 20), true);

Try[testnoerror]("test 43", UnivariatePolynomialOverPowerSeries(y + 1/2*z^2 + z^3, z), 'assign'='f');
Try("test 44", test_hf(f, 20), true);

Try[testnoerror]("test 45", UnivariatePolynomialOverPowerSeries(z^3 + (3 + x)*z^2 - (y^2 + 6)*z - 8*(x + y + 1), z), 'assign'='f');
Try("test 46", test_hf(f, 20), true);

Try[testnoerror]("test 47", UnivariatePolynomialOverPowerSeries(z^4 - (y^2 + 26)*z^2 + 25*(x + 1), z), 'assign'='f');
Try("test 48", test_hf(f, 20), true);

Try[testnoerror]("test 49", UnivariatePolynomialOverPowerSeries(z^4 + (8 - x)*z^3 + 18*z^2 - 27, z), 'assign'='f');
Try("test 50", test_hf(f, 20), true);

Try[testnoerror]("test 51", UnivariatePolynomialOverPowerSeries(z^4 + (4 + x)*z^3 - 7*z^2 - 10*z, z), 'assign'='f');
Try("test 52", test_hf(f, 20), true);

Try[testnoerror]("test 53", UnivariatePolynomialOverPowerSeries(z^2 + x - 9*z, z), 'assign'='f');
Try("test 54", test_hf(f, 20), true);

Try[testnoerror]("test 55", UnivariatePolynomialOverPowerSeries(z^2 + x - 9*z, z), 'assign'='f');
Try("test 56", test_hf(f, 5), true);
Try("test 57", test_hf(f, 10), true);
Try("test 58", test_hf(f, 15), true);


#end test

