#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
LSO := MultivariatePowerSeries:-LaurentSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", LSO(pso, ord, ordCV, rays, e), 'assign'='lso1');

pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 2", LSO(pso, mp, e), 'assign'='lso2');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1)*y^2, v=y]; 
e := [x=3, y=2];

Try[testnoerror]("test 3", LSO(pso, mp, e), 'assign'='lso3');

ord := [x, y];
pso := PowerSeries(1);
ordCV := [];
rays := [];

Try[testnoerror]("test 4", LSO(pso, ord, ordCV, rays), 'assign'='lsoOne');

ord := [x, y];
pso := PowerSeries(0);

Try[testnoerror]("test 5", LSO(pso, ord, ordCV, rays), 'assign'='lsoZero');

ord := [x, y];
pso := PowerSeries(2);

Try[testnoerror]("test 6", LSO(pso, ord, ordCV, rays), 'assign'='lsoCons');

# Test Multiply 1
Try[testnoerror]("test 7", LSO:-BinaryMultiply(lsoZero, lsoZero), 'assign'='lsMult');

Try("test 8", lsMult:-ApproximatelyEqual(lsMult, lsoZero), true);

Try[testnoerror]("test 9", LSO:-BinaryMultiply(lsoOne, lsoZero), 'assign'='lsMult');

Try("test 10", lsMult:-ApproximatelyEqual(lsMult, lsoZero), true);

Try[testnoerror]("test 11", LSO:-BinaryMultiply(lsoOne, lsoCons), 'assign'='lsMult');

Try("test 12", lsMult:-ApproximatelyEqual(lsMult, lsoCons), true);

Try[testnoerror]("test 13", LSO:-BinaryMultiply(lso1, lsoOne), 'assign'='lsMult');

Try("test 14", lsMult:-ApproximatelyEqual(lsMult, lso1, 5), true);

# Test Multiply 2
Try[testnoerror]("test 15", LSO:-BinaryMultiply(lso1, lso2), 'assign'='lsMult');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)*(-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5);
mp := [x=x, y=y];
e := [x=-7, y=-6];
lsocheck := LSO(poly, mp, e);
Try("test 16", lsMult:-ApproximatelyEqual(lsMult, lsocheck, 4, mode=absolute), true);

Try[testnoerror]("test 17", LSO:-BinaryMultiply(lso1, lso3), 'assign'='lsMult');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-5, y=-2];
lsocheck := LSO(poly, mp, e);
Try("test 18", lsMult:-ApproximatelyEqual(lsMult, lsocheck, 4, mode=absolute), true);

Try[testnoerror]("test 19", LSO:-BinaryMultiply(lso2, lso3), 'assign'='lsMult');

poly := (-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5)*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-2, y=-4];
lsocheck := LSO(poly, mp, e);
Try("test 20", lsMult:-ApproximatelyEqual(lsMult, lsocheck, 7, mode=absolute), true);

#end test