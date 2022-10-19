#test 500
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso1');

pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 2", PuSO(pso, mp, e), 'assign'='puso2');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1)*y^2, v=y]; 
e := [x=3, y=2];

Try[testnoerror]("test 3", PuSO(pso, mp, e), 'assign'='puso3');

ord := [x, y];
pso := PowerSeries(1);
ordCV := [];
rays := [];

Try[testnoerror]("test 4", PuSO(pso, ord, ordCV, rays), 'assign'='pusoOne');

ord := [x, y];
pso := PowerSeries(0);

Try[testnoerror]("test 5", PuSO(pso, ord, ordCV, rays), 'assign'='pusoZero');

ord := [x, y];
pso := PowerSeries(2);

Try[testnoerror]("test 6", PuSO(pso, ord, ordCV, rays), 'assign'='pusoCons');

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 7", PuSO(pso, ord, ordCV, rays), 'assign'='puso4');

# Test truncate
poly := expand(x^(-5)*y^3*(-x^(10)*y^(-5) + x^(8)*y^(-4) -x^(6)*y^(-3) +x^(4)*y^(-2) -x^2*y^(-1)+1 ));
Try("test 8", PuSO:-Truncate(puso1, 9, mode=absolute), poly);

poly := expand(x^3*y^(-4)*(-x^(-5)*y^(10) + x^(-4)*y^(8) -x^(-3)*y^(6) +x^(-2)*y^(4) -x^(-1)*y^2 +1));
Try("test 9", PuSO:-Truncate(puso2, 9, mode=absolute), poly);

Try("test 10", PuSO:-Truncate(puso3, 40, mode=absolute), expand(x^3*y^2*(2+2*(x^(-1)*y^2+y))));

Try("test 11", PuSO:-Truncate(pusoOne, 40, mode=absolute), 1);

Try("test 12", PuSO:-Truncate(pusoZero, 40, mode=absolute), 0);

Try("test 13", PuSO:-Truncate(pusoCons, 40, mode=absolute), 2);

Try[testnoerror]("test 14", PuSO:-Truncate(puso4, 40, mode=absolute), 'assign'='fun');

#########################
## Puiseux Series Test ##
#########################
# Truncate absolute
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1/2]];

Try[testnoerror]("test 15", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso5');

poly := expand(x^(-5)*y^3*(x^(12)*y^(-6/2)-x^(10)*y^(-5/2) + x^(8)*y^(-4/2) -x^(6)*y^(-3/2) +x^(4)*y^(-2/2) -x^2*y^(-1/2)+1 ));

Try("test 16", PuSO:-Truncate(puso5, 9, mode=absolute), poly);

#end test