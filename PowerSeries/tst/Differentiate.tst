#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
Differentiate := MultivariatePowerSeries:-PowerSeriesObject:-Differentiate:
kernelopts(opaquemodules=true):

Try[testnoerror]("test 1", PowerSeries(1+x), 'assign'='ps1');
Try[testnoerror]("test 2", Differentiate(ps1,x),'assign'='show');
Try[testnoerror]("test 3", PowerSeries(exp(x)), 'assign'='eps');
Try[testnoerror]("test 4", PowerSeries(sin(x)),'assign'='sinps');
Try("test 5", HomogeneousPart(Differentiate(eps,x),3),1/6*x^3);
Try("test 6", HomogeneousPart(Differentiate(sinps,x),2),-1/2*x^2);



#end test
