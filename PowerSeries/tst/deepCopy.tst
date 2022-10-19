#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

# Setup
Try[testnoerror]("test 1", PSO:-GeometricSeries([x, y]), 'assign'='psGeo');

# Test DeepCopy 
Try[testnoerror]("test 2", PSO:-DeepCopy(psGeo), 'assign'='psGeo2');
Try("test 3", PSO:-ApproximatelyEqual(psGeo, psGeo2, ':-force'), true);
Try[testnoerror]("test 4", PSO:-HomogeneousPart(psGeo, 10));
Try("test 5", PSO:-GetPrecision(psGeo), 10);
Try("test 6", PSO:-GetPrecision(psGeo2), 1);

kernelopts(opaquemodules=false):
Try("test 7", upperbound(psGeo:-hpoly), 10);
kernelopts(opaquemodules=true):

kernelopts(opaquemodules=false):
Try("test 8", upperbound(psGeo2:-hpoly), 1);
kernelopts(opaquemodules=true):

#end test

