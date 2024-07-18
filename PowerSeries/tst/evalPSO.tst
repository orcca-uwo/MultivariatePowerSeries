#test
with(TestTools):

Try[testnoerror]("test eval 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup
Try[testnoerror]("test eval 1", PSO:-Zero(), 'assign'='psZero');
Try[testnoerror]("test eval 2", PSO:-One(), 'assign'='psOne');
Try[testnoerror]("test eval 3", PSO:-Constant(5), 'assign'='psConst');
Try[testnoerror]("test eval 4", PSO:-FromPolynomial(2 + 1/3*(x + y)), 'assign'='psP');
Try[testnoerror]("test eval 5", PSO:-Inverse(psP), 'assign'='ps_inv');
Try("test eval 6", PSO:-Truncate(eval(w+1, w=psZero)), 1);
Try("test eval 7", PSO:-Truncate(eval(w+1, w=psOne)), 2);
Try("test eval 8", PSO:-Truncate(eval(w^2+1, w=psConst)), 26);
Try("test eval 9", PSO:-Truncate(eval(w^2+1, w=psP)), 5);
Try("test eval 10", PSO:-Truncate(eval(w^2+1, w=psP), 20), 5+4/3*x+4/3*y+1/9*x^2+2/9*x*y+1/9*y^2);
Try("test eval 11", PSO:-Truncate(eval(w+1, w=ps_inv)), 3/2);
Try("test eval 12", PSO:-Truncate(eval(w+1, w=ps_inv), 1), 3/2-1/12*x-1/12*y);

Try("test eval 13", type(eval(w^3+1, w=psZero), PSO), true);
Try("test eval 14", type(eval(w^3+1, w=ps_inv), PSO), true);

#end test
