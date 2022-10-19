#test 200
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

######
# Laurent Series Test
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

# Test Multiply 1
Try[testnoerror]("test 7", PuSO:-BinaryMultiply(pusoZero, pusoZero), 'assign'='pusMult');

Try("test 8", pusMult:-ApproximatelyEqual(pusMult, pusoZero), true);

Try[testnoerror]("test 9", PuSO:-BinaryMultiply(pusoOne, pusoZero), 'assign'='pusMult');

Try("test 10", pusMult:-ApproximatelyEqual(pusMult, pusoZero), true);

Try[testnoerror]("test 11", PuSO:-BinaryMultiply(pusoOne, pusoCons), 'assign'='pusMult');

Try("test 12", pusMult:-ApproximatelyEqual(pusMult, pusoCons), true);

Try[testnoerror]("test 13", PuSO:-BinaryMultiply(puso1, pusoOne), 'assign'='pusMult');

Try("test 14", pusMult:-ApproximatelyEqual(pusMult, puso1, 5), true);

# Test Multiply 2
Try[testnoerror]("test 15", PuSO:-BinaryMultiply(puso1, puso2), 'assign'='pusMult');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)*(-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5);
mp := [x=x, y=y];
e := [x=-7, y=-6];
pusocheck := PuSO(poly, mp, e);
Try("test 16", pusMult:-ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);

Try[testnoerror]("test 17", PuSO:-BinaryMultiply(puso1, puso3), 'assign'='pusMult');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-5, y=-2];
pusocheck := PuSO(poly, mp, e);
Try("test 18", pusMult:-ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);

Try[testnoerror]("test 19", PuSO:-BinaryMultiply(puso2, puso3), 'assign'='pusMult');

poly := (-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5)*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-2, y=-4];
pusocheck := PuSO(poly, mp, e);
Try("test 20", pusMult:-ApproximatelyEqual(pusMult, pusocheck, 7, mode=absolute), true);

#########################
## Puiseux Series Test ##
#########################
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/2,0], [1,-1/3]];

Try[testnoerror]("test 21", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso4');

ord := [x, y];
pso := PowerSeries(2+u*v+u^2*v+u*v^2+v^3);
ordCV := [u, v];
e := [];
rays := [[1,0], [1,-1/2]];

Try[testnoerror]("test 22", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso5');

ord := [x, y];
pso := PowerSeries(1+u*v+v);
ordCV := [u, v];
e := [];
rays := [[1/2,0], [1,-1/3]];

Try[testnoerror]("test 23", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso6');

Try[testnoerror]("test 24", PuSO:-BinaryMultiply(puso5, puso6), 'assign'='pusMult');
Try[testnoerror]("test 24.1", puso5*puso6, 'assign'='pusMultS');

p := normal(eval(2+u*v+u^2*v+u*v^2+v^3, [u=x,v=x*y^(-1/2)]));
q :=normal(eval(1+u*v+v, [u=x^(1/2), v=x*y^(-1/3)]));

Try[verify,normal]("test 25", p*q- PuSO:-Truncate(pusMult, 50), 0);
Try[verify,normal]("test 25.1", p*q- PuSO:-Truncate(pusMultS, 50), 0);

Try[testnoerror]("test 26", PuSO:-BinaryMultiply(puso4, puso6), 'assign'='pusMult');
Try[testnoerror]("test 26.1", puso4*puso6, 'assign'='pusMultS');

poly := (1-(u*v)+ (u*v)^2 -(u*v)^3+ (u*v)^4 -(u*v)^5)*(1+u*v+v);
mp := [u=x^(1/2), v=x*y^(-1/3)];
e := [x=-5, y=3];
pusocheck := PuSO(poly, mp, e);
Try("test 27", pusMult:-ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);
Try("test 27.1", pusMult:-ApproximatelyEqual(pusMultS, pusocheck, 4, mode=absolute), true);

Try[testnoerror]("test 28", PuSO:-BinaryMultiply(puso4, puso4), 'assign'='pusMult');

poly := (1-(u*v)+ (u*v)^2 -(u*v)^3+ (u*v)^4 -(u*v)^5)^2;
mp := [u=x^(1/2), v=x*y^(-1/3)];
e := [x=-10, y=6];
pusocheck := PuSO(poly, mp, e);
Try("test 29", pusMult:-ApproximatelyEqual(pusMult, pusocheck, 4, mode=absolute), true);

# Exp power series.
ord := [x, y];
pso := PowerSeries(d -> (u+v)^d/d!, analytic=exp(u+v));
ordCV := [u, v];
rays := [[1/4,0], [3/5,-2/5]];
e := [];

Try[testnoerror]("test 30", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus_exp1');

ord := [x, y];
pso := PowerSeries(d -> (v)^d/d!, analytic=exp(v));
ordCV := [v];
rays := [[3/5,-2/5]];
e := [];

Try[testnoerror]("test 31", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus_exp2');

ord := [x, y];
pso := PowerSeries(d -> (u)^d/d!, analytic=exp(u));
ordCV := [u];
rays := [[1/4,0]];
e := [];

Try[testnoerror]("test 32", PuSO(pso, ord, ordCV, rays, e), 'assign'='pus_exp3');

Try[testnoerror]("test 33", PuSO:-BinaryMultiply(pus_exp2, pus_exp3), 'assign'='pusMult');

Try("test 34", PuSO:-ApproximatelyEqual(pusMult, pus_exp1, 5, mode=absolute), true);

#end test