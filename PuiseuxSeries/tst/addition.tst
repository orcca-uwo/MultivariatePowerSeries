#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
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

Try[testnoerror]("test 4", PuSO(pso), 'assign'='pusoOne');

pso := PowerSeries(0);

Try[testnoerror]("test 5", PuSO(pso), 'assign'='pusoZero');

pso := PowerSeries(2);

Try[testnoerror]("test 6", PuSO(pso), 'assign'='pusoCons');

# Test Addition 1
Try[testnoerror]("test 7", PuSO:-BinaryAdd(pusoZero, pusoZero), 'assign'='pusAdd');
Try[testnoerror]("test 7.1", pusoZero+pusoZero, 'assign'='pusAddS');

Try("test 8", pusAdd:-ApproximatelyEqual(pusAdd, pusoZero), true);
Try("test 8.1", pusAdd:-ApproximatelyEqual(pusAddS, pusoZero), true);

Try[testnoerror]("test 9", PuSO:-BinaryAdd(pusoOne, pusoZero), 'assign'='pusAdd');
Try[testnoerror]("test 9.1", pusoOne+ pusoZero, 'assign'='pusAddS');

Try("test 10", pusAdd:-ApproximatelyEqual(pusAdd, pusoOne), true);
Try("test 10.1", pusAdd:-ApproximatelyEqual(pusAddS, pusoOne), true);

Try[testnoerror]("test 11", PuSO:-BinaryAdd(pusoOne, pusoCons), 'assign'='pusAdd');

pso := PowerSeries(3);
pusoCheck :=PuSO(pso);
Try("test 12", pusAdd:-ApproximatelyEqual(pusAdd, pusoCheck), true);

Try[testnoerror]("test 13", PuSO:-BinaryAdd(puso1, pusoOne), 'assign'='pusAdd');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^5*y^2;
mp := [x=x, y=y];
e := [x=-5, y=-2];
pusocheck := PuSO(poly, mp, e);
Try("test 14", pusAdd:-ApproximatelyEqual(pusAdd, pusocheck, 5, mode=absolute), true);

# Test Addition 2
Try[testnoerror]("test 15", PuSO:-BinaryAdd(puso1, puso2), 'assign'='pusAdd');

poly := y^2*(-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^3*(-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5);
mp := [x=x, y=y];
e := [x=-5, y=-4];
pusocheck := PuSO(poly, mp, e);
Try("test 16", pusAdd:-ApproximatelyEqual(pusAdd, pusocheck, 7, mode=absolute), true);

Try[testnoerror]("test 17", PuSO:-BinaryAdd(puso1, puso3), 'assign'='pusAdd');

poly := (-x^(10)+x^8*y-x^6*y^2+ x^4*y^3 -x^2*y^4+y^5)+ x^5*y^2*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-5, y=-2];
pusocheck := PuSO(poly, mp, e);
Try("test 18", pusAdd:-ApproximatelyEqual(pusAdd, pusocheck, 7, mode=absolute), true);

Try[testnoerror]("test 19", PuSO:-BinaryAdd(puso2, puso3), 'assign'='pusAdd');

poly := (-y^(10)+ x*y^8-x^2*y^6+ x^3*y^4-x^4*y^2+x^5)+ x^2*y^4*(2*x^3*y^3+2*x^2*y^4+2*x^3*y^2);
mp := [x=x, y=y];
e := [x=-2, y=-4];
pusocheck := PuSO(poly, mp, e);
Try("test 20", pusAdd:-ApproximatelyEqual(pusAdd, pusocheck, 7, mode=absolute), true);

#########################
## Puiseux Series Test ##
#########################
pso := PowerSeries(1/(1+u));
mp := [u=x^(-1/3)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 21", PuSO(pso, mp, e), 'assign'='puso4');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1/2)*y, v=y];
e := [x=3, y=2];

Try[testnoerror]("test 22", PuSO(pso, mp, e), 'assign'='puso5');

Try[testnoerror]("test 23", PuSO:-BinaryAdd(puso4, puso5), 'assign'='pusAdd');
Try[testnoerror]("test 23.1", puso4 + puso5, 'assign'='pusAddS');

pso := 1 - u^2+u^4-u^6 + u^8 +2*v^6 +2*u^3*v^4+2*v^7;
mp := [u=x^(-1/6)*y, v=y];
e := [x=3, y=-4];
pusocheck := PuSO(pso, mp, e);
Try("test 24", PuSO:-ApproximatelyEqual(pusAdd, pusocheck, 7, mode=absolute), true);
Try("test 24.1", PuSO:-ApproximatelyEqual(pusAddS, pusocheck, 7, mode=absolute), true);

#########################
## Extra Tests         ##
#########################
## PuSO plus PSO plus algebraic
pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1/2)*y, v=y];
e := [x=3/2, y=5];

Try[testnoerror]("test extra add 1", PuiseuxSeries(pso, mp, e), 'assign'='puso1');

e := [x=-5, y=3/2];
Try[testnoerror]("test extra add 2", PuiseuxSeries(pso, mp, e), 'assign'='puso2');

e := [x=5, y=-3/2];
Try[testnoerror]("test extra add 3", PuiseuxSeries(pso, mp, e), 'assign'='puso3');

Try[testnoerror]("test extra add 4", PuiseuxSeries(1), 'assign'='pusoONE');
Try[testnoerror]("test extra add 5", PuiseuxSeries(x), 'assign'='pusoX');

Try[testnoerror]("test extra add 6", pusoONE + puso1, 'assign'='pusoAdd');
ana := 1+ (2*y+2/x^(1/2)*y+2)*x^(3/2)*y^5;
Try[verify,normal]("test extra add 7", GetAnalyticExpression(pusoAdd), ana);

Try[testnoerror]("test extra add 8", pusoONE + puso2, 'assign'='pusoAdd');
ana := 1+ (2*y+2/x^(1/2)*y+2)/x^5*y^(3/2);
Try[verify,normal]("test extra add 9", GetAnalyticExpression(pusoAdd), ana);

Try[testnoerror]("test extra add 10", pusoONE + puso3, 'assign'='pusoAdd');
ana := 1+ (2*y+2/x^(1/2)*y+2)*x^5/y^(3/2);
Try[verify,normal]("test extra add 11", GetAnalyticExpression(pusoAdd), ana);

Try[testnoerror]("test extra add 12", pusoX + puso1, 'assign'='pusoAdd');
ana := x+ (2*y+2/x^(1/2)*y+2)*x^(3/2)*y^5;
Try[verify,normal]("test extra add 13", GetAnalyticExpression(pusoAdd), ana);

Try[testnoerror]("test extra add 14", pusoX + puso2, 'assign'='pusoAdd');
ana := x+ (2*y+2/x^(1/2)*y+2)/x^5*y^(3/2);
Try[verify,normal]("test extra add 15", GetAnalyticExpression(pusoAdd), ana);

Try[testnoerror]("test extra add 16", pusoX + puso3, 'assign'='pusoAdd');
ana := x+ (2*y+2/x^(1/2)*y+2)*x^5/y^(3/2);
Try[verify,normal]("test extra add 17", GetAnalyticExpression(pusoAdd), ana);

## subtraction.
Try[testnoerror]("test sub ", PuSO:-BinarySub(puso1, puso1), 'assign'='pusSub');
Try[verify,normal]("test sub ", GetAnalyticExpression(pusSub), 0);

Try[testnoerror]("test sub ", PuSO:-BinarySub(puso2, puso2), 'assign'='pusSub');
Try[verify,normal]("test sub ", GetAnalyticExpression(pusSub), 0);

Try[testnoerror]("test sub ", PuSO:-BinarySub(puso3, puso3), 'assign'='pusSub');
Try[verify,normal]("test sub ", GetAnalyticExpression(pusSub), 0);

Try[testnoerror]("test sub ", PuSO:-BinarySub(puso4, puso4), 'assign'='pusSub');
Try[verify,normal]("test sub ", GetAnalyticExpression(pusSub), 0);

Try[testnoerror]("test sub ", PuSO:-BinarySub(pusoOne, pusoOne), 'assign'='pusSub');
Try[verify,normal]("test sub ", GetAnalyticExpression(pusSub), 0);

#end test