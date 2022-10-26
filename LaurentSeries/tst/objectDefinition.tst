#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
LSO := MultivariatePowerSeries:-LaurentSeriesObject:
kernelopts(opaquemodules=true):

ord := [x, y];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", LSO(pso, ord, ordCV, rays, e), 'assign'='lso');

#error
rays := [[1,-1]];
Try[testerror]("test 2", LSO(pso, ord, ordCV, rays, e));

#ok
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
rays := [[1,0], [1,-1]];
Try[testnoerror]("test 3", LSO(pso, ord, ordCV, rays, e));

#error
ordCV := [x, y];
Try[testerror]("test 4", LSO(pso, ord, ordCV, rays, e));

#ok
rays := [[1,-1]];
pso := PowerSeries(1/(1+x));
Try[testnoerror]("test 5", LSO(pso, ord, rays, e), 'assign'='lso');

#error
rays := [[1,0,1]];
pso := PowerSeries(1/(1+x));
Try[testerror]("test 6", LSO(pso, ord, rays, e), 'assign'='lso');

#ok
e := [x=-5, y=3];
rays := [[0,1]];
pso := PowerSeries(1/(1+x));
Try[testnoerror]("test 7", LSO(pso, ord, rays, e), 'assign'='lso');

#error
rays := [[0,-1, 3]];
pso := PowerSeries(1/(1+x));
Try[testerror]("test 8", LSO(pso, ord, rays, e));

#error
ord := [x, y, z];
pso := PowerSeries(1/(1+x*y));
ordCV := [x, y];
e := [x=-5, y=3];
rays := [[1,0,0], [0,-1, 0]];
Try[testerror]("test 9", LSO(pso, ord, ordCV, rays, e));

#[1,0,0] > [0,0,0] true
Try("test 10", LSO:-GrevLexGreater(lso, [1,0,0], [0,0,0]), true);

#[1,0,0] < [0,0,0] false
Try("test 11", LSO:-GrevLexGreater(lso, [0,0,0], [1,0,0]), false);

#[-1,1] > [0,0] false
Try("test 12", LSO:-GrevLexGreater(lso, [-1,1], [0,0]), false);
Try("test 13", LSO:-Positive(lso, [-1,1]), false);

#[1,-1] > [0,0] true
Try("test 14", LSO:-GrevLexGreater(lso, [1,-1], [0,0]), true);
Try("test 15", LSO:-Positive(lso, [1,-1]), true);

#[1,-1,1] > [0,0,1] true
Try("test 16", LSO:-GrevLexGreater(lso, [1,-1,1], [0,0,1]), true);

#[1,-1,1] > [0,0,0] true
Try("test 17", LSO:-Positive(lso, [1,-1,1]), true);

#[-1,-1,-1] > [0,0,0] false
Try("test 18", LSO:-Positive(lso, [-1,-1,-1]), false);

#####################################
pso := PowerSeries(1/(1+u*v));
mp := [u=x^(-1)*y^2, v=y];
e := [x=3, y=-4];

Try[testnoerror]("test 19", LSO(pso, mp, e), 'assign'='lso');

mp := [u=x^(-1)*y*z, v=y^3*z^(-2)];
e := [y=-4, x=3];
Try[testnoerror]("test 20", LSO(pso, mp, e), 'assign'='lso');

mp := [x=x^(-1)*y^2, y=y];
Try[testerror]("test 21", LSO(pso, mp, e), 'assign'='lso');

e := [x=3, y=-4, w=4];
Try[testerror]("test 21", LSO(pso, mp, e), 'assign'='lso');

#######################################
q := u*v/(u^2 +u^2*v);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 22", LSO(q, mp, e), 'assign'='lso');

ana := normal(x^3*y/(1+x*y^(-1)));
Try("test 23", normal(LSO:-GetAnalyticExpression(lso)), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x*y^(-1)];
e := [x=3,y=1];
Try[testnoerror]("test 24", LSO(p, mp, e), 'assign'='lsocheck');

Try("test 25", LSO:-ApproximatelyEqual(lso, lsocheck, 4, mode=absolute), true);
##
q := u*v/(u^2*(1 +v));
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 26", LSO(q, mp, e), 'assign'='lso');

ana := normal(x^3*y/(1+x*y^(-1)));
Try("test 27", normal(LSO:-GetAnalyticExpression(lso)), ana);

Try("test 28", LSO:-ApproximatelyEqual(lso, lsocheck, 4, mode=absolute), true);

##
q := u/(u^2*v +u^2*v^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 29", LSO(q, mp, e), 'assign'='lso');

ana := normal(x*y^3/(1+x*y^(-1)));
Try("test 30", normal(LSO:-GetAnalyticExpression(lso)), ana);

p := 1-u+u^2-u^3+u^4;
mp := [u= x*y^(-1)];
e := [x=1,y=3];
Try[testnoerror]("test 31", LSO(p, mp, e), 'assign'='lsocheck');

Try("test 32", LSO:-ApproximatelyEqual(lso, lsocheck, 4, mode=absolute), true);
##
q := u/(u + v +u^2);
mp := [u=x, v=x*y^(-1)]; 
e := [x=3, y=2];

Try[testnoerror]("test 33", LSO(q, mp, e), 'assign'='lso');

ana := normal(x^3*y^3/(1+y+x*y));
Try("test 34", normal(LSO:-GetAnalyticExpression(lso)), ana);

p := 1-y+y^2-x*y;
mp := [x=x, y=y];
e := [x=3,y=3];
Try[testnoerror]("test 31", LSO(p, mp, e), 'assign'='lsocheck');

Try("test 35", LSO:-ApproximatelyEqual(lso, lsocheck, 4, mode=absolute), true);

##
q := 1/(u*v + u^2*v +u^3*v*w);
mp := [u=x, v=x*y^(-1), w=x^2*y^(-1)]; 

Try[testnoerror]("test 36", LSO(q, mp), 'assign'='lso');

ana := normal(y/(x^2+x^3+x^6*y^(-1)));
Try("test 37", normal(LSO:-GetAnalyticExpression(lso)), ana);

p := 1-u+u^2-u^2*w-u^3 -u^5 -3*u^4*w + u^4 + 2*u^3 *w ;
mp := [u= x, w=x^2*y^(-1)];
e := [x=-2,y=1];
Try[testnoerror]("test 38", LSO(p, mp, e), 'assign'='lsocheck');

Try("test 39", LSO:-ApproximatelyEqual(lso, lsocheck, 4, mode=absolute), true);

###########
# Without constant - rational constructor.
q := u/(u + v + w +u^2*w);
mp := [u=x^(1), v=x^(12)*y^(8), w=x^(2)*y^(-1)]; 
e := [x=3, y=-2];

Try[testnoerror]("test 40", LSO(q, mp, e), 'assign'='lso13');
#end test