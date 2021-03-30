#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
kernelopts(opaquemodules=true):

Try[testnoerror]("test 1", PSO:-FromProcedure((d -> 0), 0, check = true, expand = true), 'assign'='ps');
Try("test 2", PSO:-ApproximatelyZero(ps, 20), true);
Try[testnoerror]("test 3", PSO:-FromProcedure((d -> 0), 0, check = false, expand = false), 'assign'='ps');
Try("test 4", PSO:-ApproximatelyZero(ps, 20), true);

foo_one := proc(d) return ifelse(d = 0, 1, 0); end proc;
Try[testnoerror]("test 5", PSO:-FromProcedure(foo_one, 1, check = true, expand = true), 'assign'='ps');
Try("test 6", PSO:-HomogeneousPart(ps, 0), 1);
Try("test 7", PSO:-HomogeneousPart(ps, 1), 0);
Try("test 8", PSO:-HomogeneousPart(ps, 5), 0);
Try("test 9", PSO:-HomogeneousPart(PSO:-Negate(ps), 20), 0);
Try[testnoerror]("test 10", PSO:-FromProcedure(foo_one, 1, check = false, expand = false), 'assign'='ps');
Try("test 11", PSO:-HomogeneousPart(ps, 0), 1);
Try("test 12", PSO:-HomogeneousPart(ps, 1), 0);
Try("test 13", PSO:-HomogeneousPart(ps, 5), 0);
Try("test 14", PSO:-HomogeneousPart(PSO:-Negate(ps), 0), -1);

foo_x := proc(d) return x^d; end proc;
Try[testnoerror]("test 15", PSO:-FromProcedure(foo_x, undefined, {x}, check = true, expand = true), 'assign'='ps');
Try("test 16", PSO:-HomogeneousPart(ps, 0), 1);
Try("test 17", PSO:-HomogeneousPart(ps, 1), x);
Try("test 18", PSO:-HomogeneousPart(ps, 5), x^5);
Try("test 19", PSO:-HomogeneousPart(PSO:-Negate(ps), 20), -x^20);
Try("test 20", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 30), 0);
Try("test 21", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 40), 2*x^40);
Try[testnoerror]("test 22", PSO:-FromProcedure(foo_x, undefined, {x}, check = false, expand = false), 'assign'='ps');
Try("test 23", PSO:-HomogeneousPart(ps, 0), 1);
Try("test 24", PSO:-HomogeneousPart(ps, 1), x);
Try("test 25", PSO:-HomogeneousPart(ps, 5), x^5);
Try("test 26", PSO:-HomogeneousPart(PSO:-Negate(ps), 20), -x^20);
Try("test 27", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 30), 0);
Try("test 28", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 40), 2*x^40);

foo_x := proc(d) return (-1)^d * (1/2)^d * y^d; end proc;
Try[testnoerror]("test 29", PSO:-FromProcedure(foo_x, undefined, {x, y}, check = true, expand = true), 'assign'='ps');
Try("test 30", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 10), 0);
Try("test 31", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 20), (1/2)^19 * y^20);
Try[testnoerror]("test 32", PSO:-FromProcedure(foo_x, undefined, {x, y}, check = false, expand = false), 'assign'='ps');
Try("test 33", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 10), 0);
Try("test 34", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 20), (1/2)^19 * y^20);

foo_sqrt_x := proc(d) return (-1)^d * (sqrt(2))^d * y^d; end proc;
Try[testnoerror]("test 35", PSO:-FromProcedure(foo_sqrt_x, undefined, {x, y}, check = true, expand = true), 'assign'='ps');
Try("test 36", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 9), 0);
Try("test 37", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 10), 64 * y^10);

foo_root_y := proc(d) return (-1)^d * (RootOf(x^2+1))^d * y^d; end proc;
Try[testnoerror]("test 35", PSO:-FromProcedure(foo_root_y, undefined, {y}, check = true, expand = true), 'assign'='ps');
Try("test 36", PSO:-HomogeneousPart(PSO:-BinaryAdd(PSO:-Negate(ps), ps), 9), 0);
Try("test 37", PSO:-HomogeneousPart(PSO:-BinaryAdd(ps, ps), 10), Algebraic:-Expand(2 * RootOf(x^2+1)^10 * y^10));

foo_x_err := proc(d) return x; end proc;
Try[testerror]("test 38", PSO:-FromProcedure(foo_x_err, undefined, {x}, check = true, expand = true),
               "the generator procedure was expected to return a polynomial in {x} of total degree "
               "0, but returned this term: x");
foo_xy_err := proc(d) return x+y; end proc;
Try[testerror]("test 39", PSO:-FromProcedure(foo_xy_err, undefined, {x, y}, check = true, expand = true),
               "the generator procedure was expected to return a polynomial in {x, y} of total degree "
               "0, but returned these terms: x+y");

#end test
