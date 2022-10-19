#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
LSO := MultivariatePowerSeries:-LaurentSeriesObject:
MakeRaysCompatible := MultivariatePowerSeries:-LaurentSeriesObject:-MakeRaysCompatible:
BinaryMakeCompatible := MultivariatePowerSeries:-LaurentSeriesObject:-BinaryMakeCompatible;
MakeOrdCompatible := MultivariatePowerSeries:-LaurentSeriesObject:-MakeOrdCompatible:
ExtendLaurentSeriesObject := MultivariatePowerSeries:-LaurentSeriesObject:-ExtendLaurentSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1,0], [1,-1]];

Try[testnoerror]("test 1", LSO(pso, ord, ordCV, rays, e), 'assign'='lso1');

pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2];
e := [x=3, y=-4];

Try[testnoerror]("test 2", LSO(pso, mp, e), 'assign'='lso2');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1)*y^2, v=y]; 
e := [x=3, y=2];

Try[testnoerror]("test 3", LSO(pso, mp, e), 'assign'='lso3');

ord := [x, y];
ordCV := [];
rays := [];
pso := PowerSeries(2);

Try[testnoerror]("test 4", LSO(pso, ord, ordCV, rays), 'assign'='lsoCons');

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 5", LSO(pso, ord, ordCV, rays, e), 'assign'='lso4');

pso := PowerSeries(1/(1+v1*v2));
e := [x3=-5, x1=3];
rays := [[2,-2,0,0,0], [1,0,-1,0,0]];
ordCV := [v1, v2];

Try[testnoerror]("test 6", LSO(pso, ord, ordCV, rays, e), 'assign'='lso5');

#################################################################################
#Test BinaryMakeCompatible
Try[testnoerror]("test 7", BinaryMakeCompatible(lso1, lso2), 'assign'='P');

Try[testnoerror]("test 8", BinaryMakeCompatible(lso4, lso5), 'assign'='P');

Try[testnoerror]("test 9", BinaryMakeCompatible(lso3, lsoCons), 'assign'='P');

Try[testerror]("test 10", BinaryMakeCompatible(lso4, lso1), 'assign'='P');


################################################################################
#Test MakeRaysCompatible
# We make the rays compatible and check if the "conversion" or the rays are 
# grevlex positive
newOrd := lso1:-GetLaurentSeriesOrder(lso1);
rays := [op(lso1:-GetRays(lso1)), op(lso2:-GetRays(lso2))];
Try[testnoerror]("test 11", MakeRaysCompatible(lso1, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 12", andseq(LSO:-Positive(lso1, r), r in result), true);

newOrd := lso4:-GetLaurentSeriesOrder(lso4);
rays := [op(lso1:-GetRays(lso4)), op(lso2:-GetRays(lso5))];
Try[testnoerror]("test 13", MakeRaysCompatible(lso4, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 14", andseq(LSO:-Positive(lso4, r), r in result), true);

newOrd := LSO:-GetLaurentSeriesOrder(lso1);
rays := [op(lso1:-GetRays(lso1)), op(LSO:-GetRays(lsoCons))];
Try[testnoerror]("test 15", MakeRaysCompatible(lso1, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 16", andseq(LSO:-Positive(lso1, r), r in result), true);

#########################################################################################
#Test MakeOrdCompatible
# Setup
ord1 := [x, y, z];
ord2 := [x, w, z];
ord3 := [y, z, x];
ord4 := [x, z];
ord5 := [a, x, c, z, b];

# test
Try("test 17", MakeOrdCompatible(LSO, ord1, ord1), [x,y,z]);

Try[testerror]("test 18", MakeOrdCompatible(LSO, ord1, ord3));

Try("test 19", MakeOrdCompatible(LSO, ord1, ord2), [x, w, y, z]);

Try("test 20", MakeOrdCompatible(LSO, ord1, ord5), [a, x, c, y, z, b]);

Try("test 21", MakeOrdCompatible(LSO, [a, b], [c, d]), [a, b, c, d]);

Try("test 22", MakeOrdCompatible(LSO, [a, c, b, z], [a, d, z]), [a, c, b, d, z]);

Try("test 23", MakeOrdCompatible(LSO, [c, d], [a, b]), [a, b, c, d]);

Try[testerror]("test 24", MakeOrdCompatible(LSO, ord3, ord5));

#############################################################################################
#Test ExtendLaurentSeriesObject
# Setup
pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2*z];
e := [x=3, y=-4];

Try[testnoerror]("test 25", LSO(pso, mp, e), 'assign'='lso6');

# Test
ord := MakeOrdCompatible(LSO, [x,y], ord1);
Try[testnoerror]("test 26", ExtendLaurentSeriesObject(lso1, ord), 'assign'='lsoExt');

Try("test 27", LSO:-ApproximatelyEqual(lso1, lsoExt, 7), true);

ord := MakeOrdCompatible(LSO, [x,y], ord4);
Try[testnoerror]("test 28", ExtendLaurentSeriesObject(lso1, ord), 'assign'='lsoExt');

Try("test 29", LSO:-ApproximatelyEqual(lso1, lsoExt, 7), true);

ord := MakeOrdCompatible(LSO, [x,y,z], ord5);
Try[testnoerror]("test 30", ExtendLaurentSeriesObject(lso6, ord), 'assign'='lsoExt');

Try("test 31", LSO:-ApproximatelyEqual(lso6, lsoExt, 7), true);
#end test