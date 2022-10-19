#test

with(TestTools):
with(MultivariatePowerSeries):

## Case 1
pso1 := x / PowerSeries(1 + x + y);
monomial1 := x;
expected1 := 1 / PowerSeries(1 + x + y);

# First with Precision(pso1) = 0
Try[type]("1.0", FactorOutMonomial(pso1, monomial1), 'object', 'assign' = "res1");
Try("1.1", ApproximatelyEqual(res1, expected1, 10, 'force'), true);

# then with Precision(pso1) ~ 10
Try[verify,'interval'('closed')]("1.2", Precision(pso1), 10 .. 15);
Try[type]("1.3", FactorOutMonomial(pso1, monomial1), 'object', 'assign' = "res1b");
Try("1.4", ApproximatelyEqual(res1b, expected1, 20, 'force'), true);

## Case 2
pso2 := PowerSeries(d -> x^d/d!, analytic=exp(x)) - 1 - x;
monomial2 := x^2;
expected2 := PowerSeries(d -> x^d/(d+2)!, 'variables' = {x});

# First with Precision(pso2) = 0
Try[type]("2.0", FactorOutMonomial(pso2, monomial2), 'object', 'assign' = "res2");
Try("2.1", ApproximatelyEqual(res2, expected2, 10, 'force'), true);

# then with Precision(pso2) ~ 10
Try[verify,'interval'('closed')]("2.2", Precision(pso2), 10 .. 15);
Try[type]("2.3", FactorOutMonomial(pso2, monomial2), 'object', 'assign' = "res2b");
Try("2.4", ApproximatelyEqual(res2b, expected2, 20, 'force'), true);

## Case 3
pso3 := x^3*y / PowerSeries(1 + x + y + z);
monomial3 := x^2*y;
expected3 := x / PowerSeries(1 + x + y + z);

# First with Precision(pso3) = 0
Try[type]("3.0", FactorOutMonomial(pso3, monomial3), 'object', 'assign' = "res3");
Try("3.1", ApproximatelyEqual(res3, expected3, 10, 'force'), true);

# then with Precision(pso3) ~ 10
Try[verify,'interval'('closed')]("3.2", Precision(pso3), 10 .. 15);
Try[type]("3.3", FactorOutMonomial(pso3, monomial3), 'object', 'assign' = "res3b");
Try("3.4", ApproximatelyEqual(res3b, expected3, 20, 'force'), true);

#end test
