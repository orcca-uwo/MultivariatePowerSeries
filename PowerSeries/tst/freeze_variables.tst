#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

a := freeze(y+1);
vars :=[a, y, x]; 

Try[testerror]("test 0", PowerSeries(1/a), 'assign'='ps1');
Try[testnoerror]("test 1", PowerSeries(1/(a+1)), 'assign'='ps1');
Try("test 2", Variables(ps1), {y+1});

Try("test 3", HomogeneousPart(ps1, 1), -y-1);
Try("test 4", HomogeneousPart(ps1, 2), (y+1)^2);
Try("test 5", HomogeneousPart(ps1, 3), -(y+1)^3);

Try("test 6", GetAnalyticExpression(ps1), 1/((y+1)+1));

Try("test 7", Truncate(ps1, 3), -(y+1)^3+(y+1)^2-y);

Try[testnoerror]("test 8", ps1+ PowerSeries(1/(x+1)), 'assign'='ps2');
Try("test 9", Variables(ps2), {y+1, x});

Try("test 10", HomogeneousPart(ps2, 1), -y-1-x);
Try("test 11", HomogeneousPart(ps2, 2), (y+1)^2+x^2);
Try("test 12", HomogeneousPart(ps2, 3), -(y+1)^3-x^3);

Try("test 13", GetAnalyticExpression(ps2), 1/(y+2)+1/(x+1));

Try("test 14", Truncate(ps2, 3), -x^3-(y+1)^3+x^2+(y+1)^2-x-y+1);

Try[testnoerror]("test 15", ps1+ PowerSeries(1/(y+1)), 'assign'='ps3');
Try("test 16", Variables(ps3), {y+1, y});

Try("test 17", HomogeneousPart(ps3, 1), -y-1-y);
Try("test 18", HomogeneousPart(ps3, 2), (y+1)^2+y^2);
Try("test 19", HomogeneousPart(ps3, 3), -(y+1)^3-y^3);

Try("test 20", GetAnalyticExpression(ps3), 1/(y+2)+1/(y+1));

Try("test 21", Truncate(ps3, 3), -y^3-(y+1)^3+y^2+(y+1)^2-y-y+1);

Try[testnoerror]("test 22", UnivariatePolynomialOverPowerSeries([ps1, ps2, ps3], 'z'), 'assign'='up1');
Try("test 23", Variables(up1), {y+1, y, x, z});

Try("test 24", GetAnalyticExpression(up1), 1/(y+2)+(1/(y+2)+1/(x+1))*z+(1/(y+2)+1/(y+1))*z^2);

Try("test 25", Truncate(up1, 1), (1-2*y)*z^2+(1-y-x)*z-y);

#end test


