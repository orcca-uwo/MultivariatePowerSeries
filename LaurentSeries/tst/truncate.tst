#test 500
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

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 7", LSO(pso, ord, ordCV, rays), 'assign'='lso4');

# Test truncate
poly := expand(x^(-5)*y^3*(-x^(10)*y^(-5) + x^(8)*y^(-4) -x^(6)*y^(-3) +x^(4)*y^(-2) -x^2*y^(-1)+1 ));
Try("test 8", LSO:-Truncate(lso1, 9, mode=absolute), poly);

poly := expand(x^3*y^(-4)*(-x^(-5)*y^(10) + x^(-4)*y^(8) -x^(-3)*y^(6) +x^(-2)*y^(4) -x^(-1)*y^2 +1));
Try("test 9", LSO:-Truncate(lso2, 9, mode=absolute), poly);

Try("test 10", LSO:-Truncate(lso3, 40, mode=absolute), expand(x^3*y^2*(2+2*(x^(-1)*y^2+y))));

Try("test 11", LSO:-Truncate(lsoOne, 40, mode=absolute), 1);

Try("test 12", LSO:-Truncate(lsoZero, 40, mode=absolute), 0);

Try("test 13", LSO:-Truncate(lsoCons, 40, mode=absolute), 2);

Try[testnoerror]("test 14", LSO:-Truncate(lso4, 40, mode=absolute), 'assign'='fun');

#end test