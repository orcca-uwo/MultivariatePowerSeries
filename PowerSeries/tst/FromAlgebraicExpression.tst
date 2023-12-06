#test

with(TestTools):

# MPL-30044
Try[type]("1.0", MultivariatePowerSeries:-PowerSeries(sin(x * y)), 'object', 'assign' = 'ps0');
Try[type]("1.1", ps0:-FromAlgebraicExpression(sin(x * y), ':-method' = ':-series'), 'object',
          'assign' = 'ps_s');
Try[type]("1.2", ps0:-FromAlgebraicExpression(sin(x * y), ':-method' = ':-diff'), 'object',
          'assign' = 'ps_d');
Try("1.3", Truncate(ps_s, 10), x*y-1/6*x^3*y^3+1/120*x^5*y^5);
Try("1.4", Truncate(ps_d, 10), x*y-1/6*x^3*y^3+1/120*x^5*y^5);


#end test
