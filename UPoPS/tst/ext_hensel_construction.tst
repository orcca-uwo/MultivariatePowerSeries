#test 150
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
UPOP := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:
firstNonzeroTerm := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-firstNonzeroTerm;
exponentList := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-exponentList;
NewtonLine := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-NewtonLine;
NewtonPolynomial := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-NewtonPolynomial;
initialFactors := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-initialFactors;
LagrangeInterpolationPolynomials := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-LagrangeInterpolationPolynomials;
MainIteration := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-MainIteration;
HenselGenerator := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-HenselGenerator;
ExtendedHenselConstructionUnivariate := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-ExtendedHenselConstructionUnivariate;
ConvertToHenselUPoP := MultivariatePowerSeries:-UnivariatePolynomialOverPowerSeriesObject:-ConvertToHenselUPoP;
force_precision := MultivariatePowerSeries:-PowerSeriesObject:-force_precision:
real_precision := MultivariatePowerSeries:-PowerSeriesObject:-real_precision:
kernelopts(opaquemodules=true):

ModuloSk := proc(p, k::nonnegint, 
                     d::nonnegint, 
                     delta_hat::integer, 
                     d_hat::nonnegint,
                     v::name, u::name, 
                     multivariate::truefalse := false, $)
    local l;
    local G := [seq(v^l*u^(k+(d-l)*delta_hat), l=0..d)];
    local T :=  plex(v,u);
    local p_eval := simplify(p, power, symbolic);

    if multivariate then 
        local r := p_eval;
        for local q in G do 
            r := evala(Prem(r, q, v));
        end do;

        return r;
    else 
        Groebner:-TrailingTerm(coeff(p_eval, v,0),plex(v,u));
        local p_mod := Groebner:-NormalForm(p_eval, G, T);

        return p_mod;
    end if;
end proc;

# To check that up - mul(up_ehc)=0 mod S_currentJ. Note that 
# currentJ is the next iteration to be computed. So, currentJ-1
# is what have been computed.
# NOTE: this test is meant to test the EHC before converting the
# output to UPoP with Puiseux series coefficient. 
testEHC:= proc(up, 
                up_ehc, 
                deg::nonnegint, $)
    # TODO: Make sure that each lifted factor has degree one in v.
    local size := numelems(up_ehc);
    for local i from 1 to size do
        up:-UpdatePrecision(up_ehc[i], deg);
    end do;
    up:-UpdatePrecision(up, deg);
    # Required information.
    local c := GetCoefficient(up_ehc[1],Degree(up_ehc[1])-1);
    local state := GetAncestor(c, ':-state');
    local k := state:-currentJ;
    local d_hat := state:-d_hat;
    local delta_hat := state:-delta_hat;
    local v := state:-v;
    local d := Degree(up);
    local u := state:-u;
    # We computing the approximation of the factors up
    # to the iteration k-1.
    local fac;
    local approx_poly := mul(seq(up:-Truncate(fac), fac in up_ehc));
    #approx_poly :=  expand(eval(approx_poly, u = u^(1/d_hat)));
    local up_as_poly := up:-Truncate(up);
    up_as_poly := eval(up_as_poly, u = u^(d_hat));
    
    # The error polynomial.
    local err_poly := simplify(up_as_poly - approx_poly);
    
    # We check that the remaining monomials are zero mod
    # the ideal Sk.
    if err_poly <> 0 then
        local err_poly_eval := eval(err_poly, u=u^d_hat);
 
        local err_poly_reduced := ModuloSk(err_poly_eval, k, d, delta_hat, d_hat, v, u);
        if err_poly_reduced <> 0 then 
         return false;
        end if;
    end if;
    return(true);
end proc:

# To check that up - mul(up_ehc)=0 mod S_currentJ. Note that 
# currentJ is the next iteration to be computed. So, currentJ-1
# is what have been computed.
testEHC1:= proc(up, 
				up_ehc, 
				deg::nonnegint, 
                returnleadingcoefficient::truefalse:= false, 
                multivariate::truefalse := false,
                $)
    # TODO: Make sure that each lifted factor has degree one in v.

    local size := numelems(up_ehc);
    for local i from 1 to size do
    	up:-UpdatePrecision(up_ehc[i], deg);
    end do;
    up:-UpdatePrecision(up, deg);

    # Required information.
    local c;
    if returnleadingcoefficient then
        c := GetCoefficient(up_ehc[2],Degree(up_ehc[2]));
    else
        c := GetCoefficient(up_ehc[1],Degree(up_ehc[1]));
    end if;
    local c_pso := c:-GetPowerSeries(c);
    local state := c_pso:-GetAncestor(c_pso, ':-state');
    local k := state:-currentJ;
    local d_hat := state:-d_hat;
    local delta_hat := state:-delta_hat;
    local v := state:-v;
    local d := Degree(up);
    local u := state:-u;

    # We computing the approximation of the factors up
    # to the iteration k-1.
    local fac;
    local approx_poly := mul(seq(up:-Truncate(fac), fac in up_ehc));

    local up_as_poly := up:-Truncate(up);
 	
 	# The error polynomial.
    local err_poly := simplify(up_as_poly - approx_poly);

    # We check that the remaining monomials are zero mod
    # the ideal Sk.
    if err_poly <> 0 then
        local err_poly_eval := eval(err_poly, u=u^d_hat);
 
        local err_poly_reduced := ModuloSk(err_poly_eval, k, d, 
                                            delta_hat, d_hat, v, u, 
                                            multivariate);
        if err_poly_reduced <> 0 then 
         return false;
        end if;
    end if;

    return(true);
end proc:

testEHC2:= proc(up, 
                up_ehc, 
                deg::nonnegint, 
                returnleadingcoefficient::truefalse:= false, $)

    for local fac in up_ehc do
        if Degree(fac)>1 then
            return false;
        end if;
    end do;

    local poly := up:-GetAnalyticExpression(up);

    if not (type(poly, ':-polynom') and poly <> undefined) then 
        error "case not cover.";
    end;

    local size := numelems(up_ehc);
    for local i from 1 to size do
        up:-UpdatePrecision(up_ehc[i], deg);
    end do;
    up:-UpdatePrecision(up, deg);

    # Required information.
    local c;
    if returnleadingcoefficient then
        c := GetCoefficient(up_ehc[2],Degree(up_ehc[2]));
    else
        c := GetCoefficient(up_ehc[1],Degree(up_ehc[1]));
    end if;
    local c_pso := c:-GetPowerSeries(c);
    local state := GetAncestor(c_pso, ':-state');
    local k := state:-currentJ;

    local d_hat := state:-d_hat;
    local delta_hat := state:-delta_hat;
    local v := state:-v;
    local d := Degree(up);
    local vars := up:-Variables(up) minus {v};
    local u := state:-u;

    # Old algorithm.
    local U := RegularChains:-UnivariatePolynomialOverPowerSeries([op(vars)], v):
    local OutputFlag := 'series';
    local point := 0;
    if state:-multivariate then 
        point := [0, 0];
    end if;
    local old_result:= U:-ExtendedHenselConstruction(poly, point, k-1, OutputFlag);
    local old_result_hat;
    local result_hat;

    local f;
    old_result_hat := [seq(v-op(f), f in old_result)];
    old_result_hat := [seq(simplify(f, power, symbolic), f in old_result_hat)];
    
    local j;
    if returnleadingcoefficient then
        result_hat := [seq(simplify(up:-Truncate(up_ehc[j])), j=2..numelems(up_ehc))];
        old_result_hat := [seq(expand(f), f in old_result_hat)];
        result_hat := [seq(expand(f), f in result_hat)];
    else
        result_hat := [seq(simplify(up:-Truncate(f)), f in up_ehc)];
    end if;

    return verify(result_hat, old_result_hat, ':-multiset');
end proc:

# Multiplicity case.
# To check that up - mul(up_ehc)=0 mod S_currentJ. Note that 
# currentJ is the next iteration to be computed. So, currentJ-1
# is what have been computed.
# returnleadingcoefficient must be false! Only works for univariate upops.
testEHC3:= proc(fac, 
                deg::nonnegint, $)
    # TODO: Make sure that each lifted factor has degree one in v.
    fac:-UpdatePrecision(fac, deg);

    local c := GetCoefficient(fac,Degree(fac));
    local c_pso := c:-GetPowerSeries(c);
    local state := c_pso:-GetAncestor(c_pso, ':-state');
    local hensel_factors := state:-upolys;
    local init_upop := state:-init_upop;

    return testEHC(init_upop, hensel_factors, 0);
end proc:

# To check that fac^k=fac mod S_1. 
testEHC4:= proc(fac, 
                deg::nonnegint, $)
    # TODO: Make sure that each lifted factor has degree one in v.
    fac:-UpdatePrecision(fac, deg);

    local c := GetCoefficient(fac,Degree(fac));
    local c_pso := c:-GetPowerSeries(c);
    local fac0 := c_pso:-GetAncestor(c_pso, ':-init_fact');
    local state := c_pso:-GetAncestor(c_pso, ':-state');
    local d := Degree(state:-init_upop);
    local d_hat := state:-d_hat;
    local delta_hat := state:-delta_hat;
    local u := state:-u;
    local v := state:-v;
    local hensel_factors := state:-upolys;
    local init_upop := state:-init_upop;

    local fac0_poly := eval(Truncate(fac), u=u^d_hat);
    local err_poly := fac0 - fac0_poly;

    if err_poly <> 0 then
        local err_poly_reduced := ModuloSk(err_poly, 1, d, delta_hat, d_hat, v, u);

        if err_poly_reduced <> 0 then 
         return false;
        end if;
    end if;

    return true;
end proc:

bproc := proc(d) 
	if d  = 10 then
		return y^(10);
	else  
		return 0;
	end if; 
end proc; 
pso := PowerSeries(bproc, variables={y});

bproc1 := proc(d) 
	if d  = 4 then
		return y^(4);
	else  
		return 0;
	end if; 
end proc; 
pso1 := PowerSeries(bproc1, variables={y});

bproc2 := proc(d) 
	if d  = 20 then
		return y^(20);
	else  
		return 0;
	end if; 
end proc; 
pso2 := PowerSeries(bproc2, variables={y});
one := PowerSeries(1);

bpro := proc(d) 
    if d  = 4 then
        return y^(4);
    else  
        return 0;
    end if; 
end proc; 
ps0 := PowerSeries(bpro, variables={y});

bpro1 := proc(d) 
    if d  = 3 then
        return -y^(3);
    else  
        return 0;
    end if; 
end proc; 
ps1 := PowerSeries(bpro1, variables={y});

bpro2 := proc(d) 
    if d  = 2 then
        return -y^(2);
    else  
        return 0;
    end if; 
end proc; 
ps2 := PowerSeries(bpro2, variables={y});

# UPOP definitions.
Try[testnoerror]("test 1", UnivariatePolynomialOverPowerSeries(x^3  + (-2 * y + y^2)*x + y, 'x'), 'assign'='p1');
Try[testnoerror]("test 2", UnivariatePolynomialOverPowerSeries(x^3  + (-2 * y*z + 3*z + 1)*x + z^2*y, 'x'), 'assign'='p2');
Try[testnoerror]("test 3", UnivariatePolynomialOverPowerSeries([pso,pso1,pso2,one], 'x'), 'assign'='p3');
poly := x^3 -x^2*y^2 -x*y^3 + y^4;
Try[testnoerror]("test 4", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p4');
poly := x^5 + x^4*y - 2*x^3*y - 2*x^2*y^2 + x*(y^2 - y^3) + y^3;
Try[testnoerror]("test 5", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p5');
Try[testnoerror]("test 6", UnivariatePolynomialOverPowerSeries(x^3  + y^2*x , 'x'), 'assign'='p6');
poly := x^4 -x^2*y^2 -x*y^3 - y^4;
Try[testnoerror]("test 7", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p7');
Try[testnoerror]("test 8", UnivariatePolynomialOverPowerSeries([ps0,ps1,ps2,one], 'x'), 'assign'='p8');
# Note that p9 is equal to p3.
Try[testnoerror]("test 9", UnivariatePolynomialOverPowerSeries(x^3+y^(20)*x^2+y^4*x+y^10, 'x'), 'assign'='p9');
# Non-monics
Try[testnoerror]("test 10", UnivariatePolynomialOverPowerSeries(y*x^2  + (-2 * y + 1)*x + y, 'x'), 'assign'='p10');
poly := y^3*(y+1)*x^2 + RootOf(x^2+1)*x*y + y^2:
Try[testnoerror]("test 11", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p11');
poly := y*x^2 + x + 1;
Try[testnoerror]("test 12", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p12');
poly := x^2 + 2\*exp(1)\*x\*y + exp(2)\*y^2 + x\*y^2;
Try[testnoerror]("test 13", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p13');
# Multivariate case.
poly := x^3 + (-2*y+z+1)*x + y:
Try[testnoerror]("test 14", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p14');
poly := x^3 + (y - z + z^2)*x^2 - (y + z + y^2 - z^2)*x + (y^2 - z^3);
Try[testnoerror]("test 15", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p15');
poly := x^2  - (y + z + y^2 - z^2)*x + (y^2 - z^3);
Try[testnoerror]("test 16", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p16');
poly := x^4 + (y - z + z^2)*x - (y + z + y^2 - z^2)*x^2 + (y^2 - z^3);
Try[testnoerror]("test 17", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p17');
# Multiplicity
poly :=x^3 + (2*y)*x^2 + (y^2+y^4)*x;
Try[testnoerror]("test 18", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p18');
poly :=x^4 + 3*(y+y^2)*x^3 + 3*(y^2+y^3)*x^2+y^3*x;
Try[testnoerror]("test 19", UnivariatePolynomialOverPowerSeries(poly, 'x'), 'assign'='p19');
bound := 10;

# firstNonzeroTerm tests.
Try("test FNZT 1", firstNonzeroTerm(pso, bound), 10);
Try("test FNZT 2", firstNonzeroTerm(pso1, bound), 4);
Try("test FNZT 3", firstNonzeroTerm(pso2, bound), FAIL);

# exponentList test.
poly := x^3  + (-2 * y + y^2)*x + y;
exp_list := [[3,0],[1,1], [1,2], [0, 1]];
Try[verify,multiset]("test EL 1", exponentList(p1, poly), exp_list);

poly := x^3  + (-2 * y*z + 3*z + 1)*x + z^2*y;
exp_list := [[3,0],[1,2], [1,1], [1,0], [0, 3]];
Try[verify,multiset]("test EL 2", exponentList(p2, poly), exp_list);

poly := x^3 + (y - z + z^2)*x^2 - (y + z + y^2 - z^2)*x + (y^2 - z^3);
exp_list := [[3,0],[2,1], [2,1], [2,2], [1,1], [1,1], [1,2],
             [1,2], [0,2], [0, 3]];
Try[verify,multiset]("test EL 3", exponentList(p15, poly), exp_list);

# NewtonLine tests.
Try[verify, normal]("test NL 1", [NewtonLine(p1, bound)], [3, 1]);
Try("test NL 2", [NewtonLine(p2, bound)], [3, 0]);
Try("test NL 3", [NewtonLine(p3, bound)], [3, 6]);
Try("test NL 4", [NewtonLine(p4, bound)], [3, 4]);
Try("test NL 5", [NewtonLine(p5, bound)], [5, 5/2]);
Try("test NL 6", [NewtonLine(p6, bound)], [3, 3]);
Try("test NL 7", [NewtonLine(p15, bound)], [3, 3/2]);

# NewtonPolynomial test.
Try("test NP 1", [NewtonPolynomial(p1, bound, x, y)], [3, 1, x^3  + y]);
Try("test NP 3", [NewtonPolynomial(p3, bound, x, y)], [3, 6, x^3 +y^4*x]);
Try("test NP 4", [NewtonPolynomial(p4, bound, x, y)], [3, 4, x^3 +y^4]);
Try("test NP 5", [NewtonPolynomial(p5, bound, x, y)], [5, 5/2, x^5-2*x^3*y+x*y^2]);
Try("test NP 6", [NewtonPolynomial(p6, bound, x, y)], [3, 3, x^3  + y^2*x]);
# Note: In the main algorithm we substitute z=z*t and y=y*t so the
# output would be x^3  + (-y-z)*t*x.
Try("test NP 7", [NewtonPolynomial(p15, bound, x, t, multivariate=true)], [3, 3/2, x^3  + (-y(0)-z(0))*t*x]);

# initialFactors tests.
Try[testnoerror]("test IF 1", [initialFactors(p4, 10, x, y)], 'assign'='result');
Try[verify, multiset]("test IF 2", result[3], [[x+1, 1], [x-RootOf(_Z^2-_Z+1), 1], [x-1+RootOf(_Z^2-_Z+1), 1]]);
Try[testnoerror]("test IF 3", [initialFactors(p5, 10, x, y)], 'assign'='result');
Try[verify, multiset]("test IF 4", result[3], [[x, 1], [x+1, 2], [x-1, 2]]);
Try[testnoerror]("test IF 5", [initialFactors(p15, 10, x, t, multivariate=true)], 'assign'='result');
Try[verify, multiset]("test IF 6", result[3], [[x, 1], [x-(y(0)+z(0))^(1/2), 1], [x+(y(0)+z(0))^(1/2), 1]]);

# ExtendedHenselConstruction univariate case.
Try[testnoerror]("test EHC 1", ExtendedHenselConstructionUnivariate(p5, 10), 'assign'='result');
result := [seq(ifelse(numelems(f)>1, ConvertToHenselUPoP(op(f)), op(f)), f in result)];

Try("test EHC 2", testEHC1(p5, result, 0), true);

hensel_fac := [x, x^2 - 2*x*y^(1/2) + y, x^2 + 2*x*y^(1/2) + y];
Try[verify, multiset]("test EHC 4", [seq(Truncate(r), r in result)], hensel_fac);

Try("test EHC 3", testEHC1(p5, result, 1), true);

Try("test EHC 5", testEHC1(p5, result, 2), true);
hensel_fac := [x+y, x^2 - 2*x*y^(1/2) + y, x^2 + 2*x*y^(1/2) + y];
Try[verify, multiset]("test EHC 6", [seq(Truncate(r, 2), r in result)], hensel_fac);

Try("test EHC 7", testEHC1(p5, result, 3), true);
#hensel_fac := [x+y, x^2 + 2*x*y^(1/2) + y -((1/4)*x*y^(3/2)+(1/2)*y^2),
#                        x^2 - 2*x*y^(1/2) + y +((1/4)*x*y^(3/2)-(1/2)*y^2)];
#Try[verify, multiset]("test EHC 8", [seq(normal(Truncate(r,4)), r in result)], hensel_fac);

Try("test EHC 9", testEHC1(p5, result, 4), true);
hensel_fac := [x+y+y^2, x^2 + 2*x*y^(1/2) + y -((1/4)*x*y^(3/2)+(1/2)*y^2) - (1/2*x*y^2+3/4*y^(5/2))-(53/64*x*y^(5/2)+9/8*y^3),
                         x^2 - 2*x*y^(1/2) + y +((1/4)*x*y^(3/2)-(1/2)*y^2)- (1/2*x*y^2-3/4*y^(5/2))+(53/64*x*y^(5/2)-9/8*y^3)];
#Try[verify, multiset]("test EHC 10", [seq(normal(Truncate(r)), r in result)], hensel_fac);

Try[testnoerror]("test EHC 11", p4:-ExtendedHenselConstruction(p4, 10), 'assign'='result');
Try("test EHC 12", testEHC1(p4, result, 0), true);
Try("test EHC 13", testEHC1(p4, result, 1), true);
Try("test EHC 14", testEHC1(p4, result, 5), true);
Try("test EHC 15", testEHC2(p4, result, 5), true);
Try("test EHC 16",testEHC1(p4, result, 10), true);
Try("test EHC 17",testEHC2(p4, result, 10), true);

Try[testnoerror]("test EHC 18", p1:-ExtendedHenselConstruction(p1, 10), 'assign'='result');
Try("test EHC 19", testEHC1(p1, result, 0), true);
Try("test EHC 20", testEHC1(p1, result, 1), true);
Try("test EHC 21", testEHC1(p1, result, 5), true);
Try("test EHC 22", testEHC2(p1, result, 5), true);
Try("test EHC 23",testEHC1(p1, result, 10), true);
Try("test EHC 24",testEHC2(p1, result, 10), true);

Try[testnoerror]("test EHC 25", p6:-ExtendedHenselConstruction(p6, 10), 'assign'='result');
Try("test EHC 26", testEHC1(p6, result, 0), true);
Try("test EHC 27", testEHC1(p6, result, 1), true);
Try("test EHC 28", testEHC1(p6, result, 5), true);
Try("test EHC 29", testEHC2(p6, result, 5), true);
Try("test EHC 30",testEHC1(p6, result, 10), true);
Try("test EHC 31",testEHC2(p6, result, 10), true);

Try[testnoerror]("test EHC 32", p7:-ExtendedHenselConstruction(p7, 10), 'assign'='result');
Try("test EHC 33", testEHC1(p7, result, 0), true);
Try("test EHC 34", testEHC1(p7, result, 1), true);
Try("test EHC 35", testEHC1(p7, result, 5), true);
Try("test EHC 36", testEHC2(p7, result, 5), true);
Try("test EHC 37",testEHC1(p7, result, 10), true);
Try("test EHC 38",testEHC2(p7, result, 10), true);

Try[testnoerror]("test EHC 39", p8:-ExtendedHenselConstruction(p8, 10), 'assign'='result');
Try("test EHC 40", testEHC1(p8, result, 0), true);
Try("test EHC 41", testEHC1(p8, result, 1), true);
Try("test EHC 42", testEHC1(p8, result, 5), true);
Try("test EHC 43",testEHC1(p8, result, 10), true);

Try[testnoerror]("test EHC 44", p9:-ExtendedHenselConstruction(p9, 10), 'assign'='result');
Try("test EHC 45", testEHC2(p9, result, 0), true);
Try("test EHC 46", testEHC2(p9, result, 1), true);
Try("test EHC 47", testEHC2(p9, result, 9), true);
Try("test EHC 48",testEHC2(p9, result, 20), true);

Try[testnoerror]("test EHC 49", p3:-ExtendedHenselConstruction(p3, 10), 'assign'='result');
Try("test EHC 50", testEHC1(p3, result, 0), true);
Try("test EHC 51", testEHC1(p3, result, 1), true);
Try("test EHC 52", testEHC1(p3, result, 8), true);
Try("test EHC 53",testEHC1(p3, result, 30), true);

Try[testnoerror]("test EHC 54", p10:-ExtendedHenselConstruction(p10, 10), 'assign'='result');
#Try("test EHC 55", testEHC1(p10, result, 0, true), true);
#Try("test EHC 56", testEHC1(p10, result, 1, true), true);
Try("test EHC 57", testEHC2(p10, result, 8, true), true);
Try("test EHC 58", testEHC2(p10, result, 30, true), true);

Try[testnoerror]("test EHC 59", p12:-ExtendedHenselConstruction(p12, 10), 'assign'='result');
Try("test EHC 60", testEHC2(p12, result, 8, true), true);

# ExtendedHenselConstructionUnivariate case with algebraic variable.
Try[testnoerror]("test EHCA 1", p11:-ExtendedHenselConstruction(p11, 10), 'assign'='result');
Try("test EHCA 2", testEHC1(p11, result, 1, true, true), true);
Try("test EHCA 3", testEHC1(p11, result, 3, true, true), true);
Try("test EHCA 4", testEHC1(p11, result, 5, true, true), true);

# ExtendedHenselConstructionUnivariate Algebraic variable. Note this is also an
# example of multiplicity on one factor!
Try[testnoerror]("test EHCA 5", p13:-ExtendedHenselConstruction(p13, 10), 'assign'='result');
Try("test EHCA 6", testEHC1(p13, result, 1, false, true), true);
Try("test EHCA 7", testEHC1(p13, result, 3, false, true), true);
Try("test EHCA 8", testEHC1(p13, result, 5, false, true), true);

# ExtendedHenselConstruction univariate case with multiplicity.
Try[testnoerror]("test EHCM 1", ExtendedHenselConstruction(p5, 10), 'assign'='result');
Try("test EHCM 2", andseq(testEHC3(fac, 0), fac in result), true);
Try("test EHCM 3", andseq(testEHC3(fac, 5), fac in result), true);
Try("test EHCM 4", andseq(testEHC3(fac, 10), fac in result), true);
Try("test EHCM 5", andseq(testEHC3(fac, 15), fac in result), true);
Try("test EHCM 6", andseq(testEHC3(fac, 20), fac in result), true);

Try[testnoerror]("test EHCM 7", ExtendedHenselConstruction(p18, 10), 'assign'='result');
Try("test EHCM 8", andseq(testEHC3(fac, 5), fac in result), true);
Try("test EHCM 9", andseq(testEHC3(fac, 10), fac in result), true);

Try[testnoerror]("test EHCM 10", ExtendedHenselConstruction(p19, 10, t), 'assign'='result');
Try("test EHCM 11", andseq(testEHC3(fac, 5), fac in result), true);
Try("test EHCM 12", andseq(testEHC3(fac, 10), fac in result), true);

# ExtendedHenselConstruction multivariate.
Try[testnoerror]("test EHCM 7", p14:-ExtendedHenselConstruction(p14, 10, t), 'assign'='result');
Try("test EHCM 8", testEHC1(p14, result, 1, false, true), true);
Try("test EHCM 9", testEHC1(p14, result, 3, false, true), true);
Try("test EHCM 10", testEHC1(p14, result, 5, false, true), true);

Try[testnoerror]("test EHCM 11", p15:-ExtendedHenselConstruction(p15, 10, t), 'assign'='result');
Try("test EHCM 12", testEHC1(p15, result, 1, false, true), true);
Try("test EHCM 13", testEHC1(p15, result, 3, false, true), true);
Try("test EHCM 14", testEHC1(p15, result, 5, false, true), true);

Try[testnoerror]("test EHCM 15", p16:-ExtendedHenselConstruction(p16, 10, t), 'assign'='result');
Try("test EHCM 16", testEHC1(p16, result, 1, false, true), true);
Try("test EHCM 17", testEHC1(p16, result, 3, false, true), true);
Try("test EHCM 18", testEHC1(p16, result, 5, false, true), true);

Try[testnoerror]("test EHCM 19", p17:-ExtendedHenselConstruction(p17, 10, t), 'assign'='result');
Try("test EHCM 20", testEHC1(p17, result, 1, false, true), true);
Try("test EHCM 21", testEHC1(p17, result, 3, false, true), true);
Try("test EHCM 22", testEHC1(p17, result, 5, false, true), true);

#############################################
r1 := x1*x2^2 + x2 + 1;
Try[testnoerror]("test homg comp coeff 1", UnivariatePolynomialOverPowerSeries(r1, 'x2'), 'assign'='p');

Try[testnoerror]("test homg comp coeff 2", ExtendedHenselConstruction(p, 
                                returnleadingcoefficient=false, output=raw), 'assign'='result');

fac := result[1][1];
Try[testnoerror]("test homg comp coeff 3", GetCoefficient(fac, 0), 'assign'='a');
Try[testnoerror]("test homg comp coeff 3.1",Negate(a));
Try[testnoerror]("test homg comp coeff 4", degree(HomogeneousPart(a,0)), 0);
Try[testnoerror]("test homg comp coeff 5", degree(HomogeneousPart(a,1)), 1);
Try[testnoerror]("test homg comp coeff 6", degree(HomogeneousPart(a,4)), 4);
Try[testnoerror]("test homg comp coeff 7", degree(HomogeneousPart(a,5)), 5);
Try[testnoerror]("test homg comp coeff 8", degree(HomogeneousPart(a,7)), 7);

############################################
# force_precision
Try[testnoerror]("test force precision 1", Inverse(PowerSeries(1/(x+1))), 'assign'='ps');
Try[testnoerror]("test force precision 2", HomogeneousPart(ps,5), 'assign'='prec5');
Try[testnoerror]("test force precision 3",Precision(ps), 5);
Try[testnoerror]("test force precision 4",force_precision(ps,4));
Try[testnoerror]("test force precision 5", HomogeneousPart(ps,5), prec5);
Try[testnoerror]("test force precision 6",Precision(ps), 4);
Try[testerror]("test force precision 7",force_precision(ps,6), "forced degree should be no larger than the real "
                                                "precision, %1, but received %2", 5, 6);

#end test
