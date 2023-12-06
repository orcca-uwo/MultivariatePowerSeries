#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

# Test Weierstrass Preparation
Try[testnoerror]("test 01", UnivariatePolynomialOverPowerSeries(1), 'assign'='upOne');
Try[testnoerror]("test 02", WeierstrassPreparation(UnivariatePolynomialOverPowerSeries(1)), 'assign'='wp1');
Try[testnoerror]("test 03", Truncate(wp1[1], 10), 'assign'='p1');
Try("test 04", evalb(p1 = 1) , true);
Try[testnoerror]("test 05", Truncate(wp1[2], 10), 'assign'='alpha1');
Try("test 06", evalb(alpha1 = 1) , true);

Try[testnoerror]("test 07", UnivariatePolynomialOverPowerSeries(RootOf(x^2+1)), 'assign'='up');
Try[testnoerror]("test 08", WeierstrassPreparation(up), 'assign'='wp1');
Try[testnoerror]("test 09", Truncate(wp1[1], 10), 'assign'='p');
Try("test 10", evalb(expand(p) = 1) , true);
Try[testnoerror]("test 11", Truncate(wp1[2], 10), 'assign'='alpha');
Try("test 12", evalb(Algebraic:-Expand(alpha) = RootOf(x^2+1)) , true);

Try[testnoerror]("test 13", UnivariatePolynomialOverPowerSeries(x^2*y^2 + sqrt(2)*x^2, x), 'assign'='up');
Try[testnoerror]("test 14", WeierstrassPreparation(up), 'assign'='wp1');
Try[testnoerror]("test 15", Truncate(wp1[1], 10), 'assign'='p');
Try("test 16", evalb(expand(p) = x^2) , true);
Try[testnoerror]("test 17", Truncate(wp1[2], 10), 'assign'='alpha');
Try("test 18", evalb(evala(Expand(alpha)) = y^2+2^(1/2)) , true);

createWPInput1 := proc(d::nonnegint)
local upoly, i;
    upoly := Array(0 .. d);
    upoly[0] := PowerSeries(x);
    upoly[1] := PowerSeries(y);
    for i from 2 to d - 1 do
        upoly[i] := PowerSeries(1);
    end do;
    upoly[d] := Inverse(PowerSeries(1 + x + y));
    return UnivariatePolynomialOverPowerSeries(upoly, 'z');
end proc;

createWPInput2 := proc(d::nonnegint)
local upoly, i;
    upoly := Array(0 .. d);
    upoly[0] := PowerSeries(x);
    upoly[1] := PowerSeries(y);
    for i from 2 to d - 1 do
        upoly[i] := PowerSeries(1);
    end do;
    upoly[d] := (GeometricSeries([x, y]));
    return UnivariatePolynomialOverPowerSeries(upoly, 'z');
end proc;

createWPInput4 := proc(d::nonnegint)
local upoly, i;
    upoly := Array(0 .. d);
    upoly[0] := PowerSeries(x);
    for i to d - 2 do
        upoly[i] := PowerSeries(1);
    end do;
    upoly[d - 1] := PowerSeries(y);
    upoly[d] := PowerSeries(1 + x + y);
    return UnivariatePolynomialOverPowerSeries(upoly, 'z');
end proc;

VerifyTools:-AddVerification(
    weierstrass = proc(pa, f, prec, ver := expand, $)
    local deg_p, deg_alpha, g, h, p, alpha;
        p, alpha := op(pa);
        h := Truncate(BinaryMultiply(p, alpha), prec);
        return verify(Truncate(BinaryMultiply(p, alpha), prec),
                      Truncate(f, prec),
                      ver);
    end proc);

Try[testnoerror]("test 19", createWPInput1(3), 'assign'='up3');
Try[testnoerror]("test 20", WeierstrassPreparation(up3), 'assign'='wp3');
Try[verify,weierstrass(5)]("test 21", [wp3], up3);
Try[verify,weierstrass(10)]("test 22", [wp3], up3);

Try[testnoerror]("test 23", createWPInput1(5), 'assign'='up5');
Try[testnoerror]("test 24", WeierstrassPreparation(up5), 'assign'='wp5');
Try[verify,weierstrass(3)]("test 25", [wp5], up5);
Try[verify,weierstrass(7)]("test 26", [wp5], up5);

Try[testnoerror]("test 27", createWPInput4(3), 'assign'='up3');
Try[testnoerror]("test 28", WeierstrassPreparation(up3), 'assign'='wp3');
Try[verify,weierstrass(5)]("test 29", [wp3], up3);
Try[verify,weierstrass(10)]("test 30", [wp3], up3);

Try[testnoerror]("test 31", createWPInput4(5), 'assign'='up5');
Try[testnoerror]("test 32", WeierstrassPreparation(up5), 'assign'='wp5');
Try[verify,weierstrass(4)]("test 33", [wp5], up5);
Try[verify,weierstrass(8)]("test 34", [wp5], up5);

Try[testnoerror]("test 35", createWPInput2(5), 'assign'='up6');
Try[testnoerror]("test 36", WeierstrassPreparation(up6), 'assign'='wp6');
Try[verify,weierstrass(4)]("test 37", [wp6], up6);
Try[verify,weierstrass(8)]("test 38", [wp6], up6);

#end test
