#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
LSO := MultivariatePowerSeries:-LaurentSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
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

# Test Addition 1
Try[testnoerror]("test 7", LSO:-BinaryAdd(lsoZero, lsoZero), 'assign'='lsAdd');

Try("test 8", lsAdd:-ApproximatelyEqual(lsAdd, lsoZero), true);

Try[testnoerror]("test 9", LSO:-BinaryAdd(lsoOne, lsoZero), 'assign'='lsAdd');

Try("test 10", lsAdd:-ApproximatelyEqual(lsAdd, lsoOne), true);

Try[testnoerror]("test 11", LSO:-BinaryAdd(lsoOne, lsoCons), 'assign'='lsAdd');

pso := PowerSeries(3);
lsoCheck :=LSO(pso, ord, ordCV, rays);
Try("test 12", lsAdd:-ApproximatelyEqual(lsAdd, lsoCheck), true);

Try[testnoerror]("test 13", LSO:-BinaryAdd(lso1, lsoOne), 'assign'='lsAdd');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^5*y^2;
mp := [x=x, y=y];
e := [x=-5, y=-2];
lsocheck := LSO(poly, mp, e);
Try("test 14", lsAdd:-ApproximatelyEqual(lsAdd, lsocheck, 5, mode=absolute), true);

# Test Addition 2
Try[testnoerror]("test 15", LSO:-BinaryAdd(lso1, lso2), 'assign'='lsAdd');

poly := y^2*(-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^3*(-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5);
mp := [x=x, y=y];
e := [x=-5, y=-4];
lsocheck := LSO(poly, mp, e);
Try("test 16", lsAdd:-ApproximatelyEqual(lsAdd, lsocheck, 7, mode=absolute), true);

Try[testnoerror]("test 17", LSO:-BinaryAdd(lso1, lso3), 'assign'='lsAdd');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^5*y^2*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-5, y=-2];
lsocheck := LSO(poly, mp, e);
Try("test 18", lsAdd:-ApproximatelyEqual(lsAdd, lsocheck, 7, mode=absolute), true);

Try[testnoerror]("test 19", LSO:-BinaryAdd(lso2, lso3), 'assign'='lsAdd');

poly := (-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5)+ x^2*y^4*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-2, y=-4];
lsocheck := LSO(poly, mp, e);
Try("test 20", lsAdd:-ApproximatelyEqual(lsAdd, lsocheck, 7, mode=absolute), true);

#end test