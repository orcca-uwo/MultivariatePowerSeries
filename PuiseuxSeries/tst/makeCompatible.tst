#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PuSO := MultivariatePowerSeries:-PuiseuxSeriesObject:
MakeRaysCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeRaysCompatible:
BinaryMakeCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-BinaryMakeCompatible;
MakeOrdCompatible := MultivariatePowerSeries:-PuiseuxSeriesObject:-MakeOrdCompatible:
ExtendPuiseuxSeriesObject := MultivariatePowerSeries:-PuiseuxSeriesObject:-ExtendPuiseuxSeriesObject:
kernelopts(opaquemodules=true):

# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
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
ordCV := [];
rays := [];
pso := PowerSeries(2);

Try[testnoerror]("test 4", PuSO(pso, ord, ordCV, rays), 'assign'='pusoCons');

ord := [x1, x2, x3, x4, x5];
pso := PowerSeries(1/(1+u1*u2*u3*u4));
e := [x1=-5, x2=2];
rays := [[-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,5,-1], [-3,4,0,-1,0]];
ordCV := [u1, u2, u3, u4];

Try[testnoerror]("test 5", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso4');

pso := PowerSeries(1/(1+v1*v2));
e := [x3=-5, x1=3];
rays := [[2,-2,0,0,0], [1,0,-1,0,0]];
ordCV := [v1, v2];

Try[testnoerror]("test 6", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso5');

#################################################################################
#Test BinaryMakeCompatible
Try[testnoerror]("test 7", BinaryMakeCompatible(puso1, puso2), 'assign'='P');

Try[testnoerror]("test 8", BinaryMakeCompatible(puso4, puso5), 'assign'='P');

Try[testnoerror]("test 9", BinaryMakeCompatible(puso3, pusoCons), 'assign'='P');

Try[testerror]("test 10", BinaryMakeCompatible(puso4, puso1), 'assign'='P');


################################################################################
#Test MakeRaysCompatible
# We make the rays compatible and check if the "conversion" or the rays are 
# grevlex positive
newOrd := puso1:-GetPuiseuxSeriesOrder(puso1);
rays := [op(puso1:-GetRays(puso1)), op(puso2:-GetRays(puso2))];
Try[testnoerror]("test 11", MakeRaysCompatible(puso1, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 12", andseq(PuSO:-Positive(puso1, r), r in result), true);

newOrd := puso4:-GetPuiseuxSeriesOrder(puso4);
rays := [op(puso1:-GetRays(puso4)), op(puso2:-GetRays(puso5))];
Try[testnoerror]("test 13", MakeRaysCompatible(puso4, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 14", andseq(PuSO:-Positive(puso4, r), r in result), true);

newOrd := PuSO:-GetPuiseuxSeriesOrder(puso1);
rays := [op(puso1:-GetRays(puso1)), op(PuSO:-GetRays(pusoCons))];
Try[testnoerror]("test 15", MakeRaysCompatible(puso1, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 16", andseq(PuSO:-Positive(puso1, r), r in result), true);

#########################################################################################
#Test MakeOrdCompatible
# Setup
ord1 := [x, y, z];
ord2 := [x, w, z];
ord3 := [y, z, x];
ord4 := [x, z];
ord5 := [a, x, c, z, b];

# test
Try("test 17", MakeOrdCompatible(PuSO, ord1, ord1), [x,y,z]);

Try[testerror]("test 18", MakeOrdCompatible(PuSO, ord1, ord3));

Try("test 19", MakeOrdCompatible(PuSO, ord1, ord2), [x, w, y, z]);

Try("test 20", MakeOrdCompatible(PuSO, ord1, ord5), [a, x, c, y, z, b]);

Try("test 21", MakeOrdCompatible(PuSO, [a, b], [c, d]), [a, b, c, d]);

Try("test 22", MakeOrdCompatible(PuSO, [a, c, b, z], [a, d, z]), [a, c, b, d, z]);

Try("test 23", MakeOrdCompatible(PuSO, [c, d], [a, b]), [a, b, c, d]);

Try[testerror]("test 24", MakeOrdCompatible(PuSO, ord3, ord5));

#############################################################################################
#Test ExtendPuiseuxSeriesObject
# Setup
pso := PowerSeries(1/(1+u));
mp := [u=x^(-1)*y^2*z];
e := [x=3, y=-4];

Try[testnoerror]("test 25", PuSO(pso, mp, e), 'assign'='puso6');

# Test
ord := MakeOrdCompatible(PuSO, [x,y], ord1);
Try[testnoerror]("test 26", ExtendPuiseuxSeriesObject(puso1, ord), 'assign'='pusoExt');

Try("test 27", PuSO:-ApproximatelyEqual(puso1, pusoExt, 7), true);

ord := MakeOrdCompatible(PuSO, [x,y], ord4);
Try[testnoerror]("test 28", ExtendPuiseuxSeriesObject(puso1, ord), 'assign'='pusoExt');

Try("test 29", PuSO:-ApproximatelyEqual(puso1, pusoExt, 7), true);

ord := MakeOrdCompatible(PuSO, [x,y,z], ord5);
Try[testnoerror]("test 30", ExtendPuiseuxSeriesObject(puso6, ord), 'assign'='pusoExt');

Try("test 31", PuSO:-ApproximatelyEqual(puso6, pusoExt, 7), true);

########################
(*
## Puiseux Test.
# Setup
ord := [x, y];
pso := PowerSeries(1/(1+u*v));
ordCV := [u, v];
e := [x=-5, y=3];
rays := [[1/3,0], [1/2,-1/2]];

Try[testnoerror]("test 32", PuSO(pso, ord, ordCV, rays, e), 'assign'='puso7');

pso := PowerSeries(2 + 2*(u + v));
mp := [u=x^(-1/2)*y^2, v=y^(1/4)]; 
e := [x=3, y=2];

Try[testnoerror]("test 33", PuSO(pso, mp, e), 'assign'='puso8');

newOrd := PuSO:-GetPuiseuxSeriesOrder(puso7);
rays := [op(PuSO:-GetRays(puso7)), op(PuSO:-GetRays(puso8))];
Try[testnoerror]("test 34", MakeRaysCompatible(puso7, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test 35", andseq(PuSO:-Positive(puso7, r), r in result), true);

Try[testnoerror]("test 36", BinaryMakeCompatible(puso7, puso8), 'assign'='P');
*)

#########################
# More test for MakeRaysCompatible
#Test MakeRaysCompatible
# We make the rays compatible and check if the "conversion" or the rays are 
# grevlex positive. We also check that the vectors in result are integer.
newOrd := puso4:-GetPuiseuxSeriesOrder(puso4);
rays := [[-2,1,8,-1,0], [3,2,1,2,3], [5,-2,-4,5,-1], [-3,5,0,-1,0],
		 [-2,1,2,-1,0], [3,2,1,2,3], [2,-2,-4,7,-1], [-3,9,0,-1,0]];
Try[testnoerror]("test extra 1", MakeRaysCompatible(puso4, rays, newOrd), 'assign'='newRays');

M := Matrix(newRays);
raysAsVectors := map(convert, rays, ':-Vector');
result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in raysAsVectors)];
Try("test extra 2", andseq(PuSO:-Positive(puso4, r), r in result), true);

Try("test extra 3", andseq(type(convert(r,list),list(nonnegint)), r in result), true);

#end test


