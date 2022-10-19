#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
kernelopts(opaquemodules=true):

ord := [x, y];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso');

#error
rays := [[1,-1]];
Try[testerror]("test 2", PuSO(pso, ord, ordCV, rays, e));

#ok
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
rays := [[1,0], [1,-1]];
Try[testnoerror]("test 3", PuSO(pso, ord, ordCV, rays, e));

#error
ordCV := [x, y];
Try[testerror]("test 4", PuSO(pso, ord, ordCV, rays, e));

#ok
rays := [[1,-1]];
pso := PowerSeries(1/(1+x));
Try[testnoerror]("test 5", PuSO(pso, ord, rays, e), 'assign'='puso');

#error
rays := [[1,0,1]];
pso := PowerSeries(1/(1+x));
Try[testerror]("test 6", PuSO(pso, ord, rays, e), 'assign'='puso');

#ok
e := [x=-5, y=3];
rays := [[0,1]];
pso := PowerSeries(1/(1+x));
Try[testnoerror]("test 7", PuSO(pso, ord, rays, e), 'assign'='puso');

#error
rays := [[0,-1, 3]];
pso := PowerSeries(1/(1+x));
Try[testerror]("test 8", PuSO(pso, ord, rays, e));

#error
ord := [x, y, z];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
e := [x=-5, y=3];
rays := [[1,0,0], [0,-1, 0]];
Try[testerror]("test 9", PuSO(pso, ord, ordCV, rays, e));

#[1,0,0] > [0,0,0] true
Try("test 10", PuSO:-GrevLexGreater(puso, [1,0,0], [0,0,0]), true);

#[1,0,0] < [0,0,0] false
Try("test 11", PuSO:-GrevLexGreater(puso, [0,0,0], [1,0,0]), false);

#[-1,1] > [0,0] false
Try("test 12", PuSO:-GrevLexGreater(puso, [-1,1], [0,0]), false);
Try("test 13", PuSO:-Positive(puso, [-1,1]), false);

#[1,-1] > [0,0] true
Try("test 14", PuSO:-GrevLexGreater(puso, [1,-1], [0,0]), true);
Try("test 15", PuSO:-Positive(puso, [1,-1]), true);

#[1,-1,1] > [0,0,1] true
Try("test 16", PuSO:-GrevLexGreater(puso, [1,-1,1], [0,0,1]), true);

#[1,-1,1] > [0,0,0] true
Try("test 17", PuSO:-Positive(puso, [1,-1,1]), true);

#[-1,-1,-1] > [0,0,0] false
Try("test 18", PuSO:-Positive(puso, [-1,-1,-1]), false);

#####################################
pso := PowerSeries(1/(1+u*v));
mp := [u=x^(-1)*y^2, v=y];
e := [x=3, y=-4];

Try[testnoerror]("test 19", PuSO(pso, mp, e), 'assign'='puso');

mp := [u=x^(-1)*y*z, v=y^3*z^(-2)];
e := [y=-4, x=3];
Try[testnoerror]("test 20", PuSO(pso, mp, e), 'assign'='puso');

mp := [x=x^(-1)*y^2, y=y];
Try[testerror]("test 21", PuSO(pso, mp, e), 'assign'='puso');

e := [x=3, y=-4, w=4];
Try[testerror]("test 21", PuSO(pso, mp, e), 'assign'='puso');

#######################################
q := u*v/(u^2 +u^2*v);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 22", PuSO(q, mp, e), 'assign'='puso');

ana := normal(x^3*y/(1+x*y^(-1)));
Try("test 23", normal(PuSO:-GetAnalyticExpression(puso)), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x*y^(-1)];
e := [x=3,y=1];
Try[testnoerror]("test 24", PuSO(p, mp, e), 'assign'='pusocheck');

Try("test 25", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);
##
q := u*v/(u^2*(1 +v));
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 26", PuSO(q, mp, e), 'assign'='puso');

ana := normal(x^3*y/(1+x*y^(-1)));
Try("test 27", normal(PuSO:-GetAnalyticExpression(puso)), ana);

Try("test 28", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);

##
q := u/(u^2*v +u^2*v^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 29", PuSO(q, mp, e), 'assign'='puso');

ana := normal(x*y^3/(1+x*y^(-1)));
Try("test 30", normal(PuSO:-GetAnalyticExpression(puso)), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x*y^(-1)];
e := [x=1,y=3];
Try[testnoerror]("test 31", PuSO(p, mp, e), 'assign'='pusocheck');

Try("test 32", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);
##
q := u/(u + v +u^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 33", PuSO(q, mp, e), 'assign'='puso');

ana := normal(x^3*y^3/(1+y+x*y));
Try("test 34", normal(PuSO:-GetAnalyticExpression(puso)), ana);

p := 1-y+y^2-x*y;
mp := [x=x, y=y];
e := [x=3,y=3];
Try[testnoerror]("test 31", PuSO(p, mp, e), 'assign'='pusocheck');

Try("test 35", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);

##
q := 1/(u*v + u^2*v +u^3*v*w);
mp := [u=x, v=x*y^(-1), w=x^2*y^(-1)]; 

Try[testnoerror]("test 36", PuSO(q, mp), 'assign'='puso');

ana := normal(y/(x^2+x^3+x^6*y^(-1)));
Try("test 37", normal(PuSO:-GetAnalyticExpression(puso)), ana);

p := 1-u+u^2-u^2*w-u^3 -u^5 -3*u^4*w + u^4 + 2*u^3 *w ;
mp := [u= x, w=x^2*y^(-1)];
e := [x=-2,y=1];
Try[testnoerror]("test 38", PuSO(p, mp, e), 'assign'='pusocheck');

Try("test 39", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);

#########################
## Puiseux Series Test ##
#########################
# Definition
q := u/(u^2*v +u^2*v^2);
mp := [u=x^(1/2), v=x^(1/3)*y^(-1/6)]; 
e := [x=3, y=2];

Try[testnoerror]("test 40", PuSO(q, mp, e), 'assign'='puso');

ana := normal(x^(13/6)*y^(13/6)/(1+x^(1/3)*y^(-1/6)));
Try("test 41", normal(PuSO:-GetAnalyticExpression(puso)), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x^(1/3)*y^(-1/6)];
e := [x=13/6,y=13/6];
Try[testnoerror]("test 42", PuSO(p, mp, e), 'assign'='pusocheck');

Try("test 43", PuSO:-ApproximatelyEqual(puso, pusocheck, 4, mode=absolute), true);


# CreatePuSOFromPSO test.
pso1 := PowerSeries(1+x*y);
pso2 := PowerSeries((x+y)/(1+x*y));
pso3 := PowerSeries((2+x^2)/(1+x*y+ x^2+y^3));

pol1 := 1+x*y;
pol2 := 1+x*y+ x*y^2- x^2+y^3 -y^4;

e := [x=-13/6,y=13/6];

Try[testnoerror]("test 44", PuSO(pso1, e), 'assign'='puso1');
Try[testnoerror]("test 45", PuSO(pso2, e), 'assign'='puso2');
Try[testnoerror]("test 46", PuSO(pso3, e), 'assign'='puso3');

Try[testnoerror]("test 47", PuSO(pol1, e), 'assign'='puso4');
Try[testnoerror]("test 48", PuSO(pol2, e), 'assign'='puso5');

Try("test 49", PuSO:-ApproximatelyEqual(puso1, puso4, 50), true);

Try[testerror]("test 50", PuSO(pol2, [[1,2,3]]));
Try[testerror]("test 51", PuSO(pol2, [[1,2],[1,1/2]]));

e := [u=-13/6,y=13/6];
Try[testerror]("test 50", PuSO(pol2, e));

#end test