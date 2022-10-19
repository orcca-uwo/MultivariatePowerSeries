#test 400

with(TestTools):
with(MultivariatePowerSeries):

pso_exp := PowerSeries(d -> (x+y)^d/d!, analytic=exp(x+y));
pso_rat := 1/PowerSeries(1 + x + y + z);
pso_pol := PowerSeries(1 + x^2 + z^3*y);

pso_exp_shift := PowerSeries(d -> exp(1) * (x+y)^d/d!, analytic=exp(1+x+y));
pso_rat_shift := 1/PowerSeries(2 + x + y + z);
pso_pol_shift := PowerSeries(1 + (x+1)^2 + z^3*y);
pso_pol_shift_y := PowerSeries(1 + x^2 + z^3*(y+1));
pso_exp_shift_both := PowerSeries(d -> exp(2) * (x+y)^d/d!, analytic=exp(2+x+y));
pso_rat_shift_both := 1/PowerSeries(3 + x + y + z);
pso_pol_shift_both := PowerSeries(1 + (x+1)^2 + z^3*(y+1));


poly_x := x^2 + x^3;
poly_x_y := x^2 + y^3;

VerifyTools:-AddVerification(
    up_to_homdeg = proc(p1, p2, d := 20, $)
    local _p1 := ifelse(type(p1, ':-object'), Truncate(p1, d),
                        add(select(t -> degree(t) <= d, convert(expand(p1), ':-list', ':-`+`'))));
    local _p2 := ifelse(type(p2, ':-object'), Truncate(p2, d),
                        add(select(t -> degree(t) <= d, convert(expand(p2), ':-list', ':-`+`'))));
        return evalb(_p1 = _p2);
    end proc);

# Special cases and errors.
Try[verify,up_to_homdeg(20)]("0.0", Substitute([], pso_exp), pso_exp);
Try[testerror]("0.1", Substitute({x=poly_x, x=poly_x_y}, pso_exp), "these left hand sides occur more than once: {x}");
Try[testerror]("0.2", Substitute({x=poly_x, y=y^2/z}, pso_exp), "invalid substitution, [y = y^2/z]");

# Substitution of single polynomials into power series.
Try[verify,up_to_homdeg(20)]("1.0.0", Substitute(x = poly_x,   pso_exp), subs(x = poly_x,   Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("1.0.1", Substitute(x = poly_x_y, pso_exp), subs(x = poly_x_y, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("1.1.0", Substitute(x = poly_x,   pso_rat), subs(x = poly_x,   Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("1.1.1", Substitute(x = poly_x_y, pso_rat), subs(x = poly_x_y, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("1.2.0", Substitute(x = poly_x,   pso_pol), subs(x = poly_x,   Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("1.2.1", Substitute(x = poly_x_y, pso_pol), subs(x = poly_x_y, Truncate(pso_pol, 20)));

Try[verify,up_to_homdeg(20)]("1.3.0", Substitute(x = 1 + poly_x,   pso_exp), subs(x = poly_x,   Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("1.3.1", Substitute(x = 1 + poly_x_y, pso_exp), subs(x = poly_x_y, Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("1.4.0", Substitute(x = 1 + poly_x,   pso_rat), subs(x = poly_x,   Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("1.4.1", Substitute(x = 1 + poly_x_y, pso_rat), subs(x = poly_x_y, Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("1.5.0", Substitute(x = 1 + poly_x,   pso_pol), subs(x = poly_x,   Truncate(pso_pol_shift, 20)));
Try[verify,up_to_homdeg(20)]("1.5.1", Substitute(x = 1 + poly_x_y, pso_pol), subs(x = poly_x_y, Truncate(pso_pol_shift, 20)));

# Substitution of multiple polynomials into power series.
Try[verify,up_to_homdeg(20)]("2.0", Substitute({x = poly_x, y = poly_x_y}, pso_exp),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("2.1", Substitute({x = poly_x, y = poly_x_y}, pso_rat),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("2.2", Substitute({x = poly_x, y = poly_x_y}, pso_pol),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_pol, 20)));

Try[verify,up_to_homdeg(20)]("2.3", Substitute({x = 1 + poly_x, y = poly_x_y}, pso_exp),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("2.4", Substitute({x = 1 + poly_x, y = poly_x_y}, pso_rat),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("2.5", Substitute({x = 1 + poly_x, y = poly_x_y}, pso_pol),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_pol_shift, 20)));

Try[verify,up_to_homdeg(20)]("2.6", Substitute({x = poly_x, y = 1 + poly_x_y}, pso_exp),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("2.7", Substitute({x = poly_x, y = 1 + poly_x_y}, pso_rat),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("2.8", Substitute({x = poly_x, y = 1 + poly_x_y}, pso_pol),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_pol_shift_y, 20)));

# Substitute zero for some variables.
Try[verify,up_to_homdeg(20)]("3.0", Substitute(x = 0, pso_exp), subs(x = 0, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("3.1", Substitute(x = 0, pso_rat), subs(x = 0, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("3.2", Substitute(x = 0, pso_pol), subs(x = 0, Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("3.3", Substitute({x = 0, y =     poly_x  }, pso_exp),
                             subs({x = 0, y = poly_x}, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("3.4", Substitute({x = 0, y =     poly_x_y}, pso_rat),
                             subs({x = 0, y = poly_x_y}, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("3.5", Substitute({x = 0, y =     poly_x  }, pso_pol),
                             subs({x = 0, y = poly_x}, Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("3.6", Substitute({x = 0, y = 1 + poly_x  }, pso_exp),
                             subs({x = 0, y = poly_x}, Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("3.7", Substitute({x = 0, y = 1 + poly_x  }, pso_rat),
                             subs({x = 0, y = poly_x}, Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("3.8", Substitute({x = 0, y = 1 + poly_x_y}, pso_pol),
                             subs({x = 0, y = poly_x_y}, Truncate(pso_pol_shift_y, 20)));

# Rename some variables.
Try[verify,up_to_homdeg(20)]("4.0", Substitute(x = y, pso_exp), subs(x = y, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("4.1", Substitute(x = y, pso_rat), subs(x = y, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("4.2", Substitute(x = y, pso_pol), subs(x = y, Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("4.3", Substitute({x = y, y =     poly_x  }, pso_exp),
                             subs({x = y, y = poly_x}, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("4.4", Substitute({x = y, y =     poly_x_y}, pso_rat),
                             subs({x = y, y = poly_x_y}, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("4.5", Substitute({x = y, y =     poly_x  }, pso_pol),
                             subs({x = y, y = poly_x}, Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("4.6", Substitute({x = y, y = 1 + poly_x  }, pso_exp),
                             subs({x = y, y = poly_x}, Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("4.7", Substitute({x = y, y = 1 + poly_x  }, pso_rat),
                             subs({x = y, y = poly_x}, Truncate(pso_rat_shift, 20)));
Try[verify,up_to_homdeg(20)]("4.8", Substitute({x = y, y = 1 + poly_x_y}, pso_pol),
                             subs({x = y, y = poly_x_y}, Truncate(pso_pol_shift_y, 20)));

# Substitute power series that are actually polynomials. (Same as the above tests, prefixed by "5.".)
Try[verify,up_to_homdeg(20)]("5.1.0.0", Substitute(x = PowerSeries(poly_x),   pso_exp), subs(x = poly_x,   Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("5.1.2.1", Substitute(x = PowerSeries(poly_x_y), pso_pol), subs(x = poly_x_y, Truncate(pso_pol, 20)));
Try[verify,up_to_homdeg(20)]("5.1.3.0", Substitute(x = PowerSeries(1 + poly_x),   pso_exp), subs(x = poly_x,   Truncate(pso_exp_shift, 20)));
Try[verify,up_to_homdeg(20)]("5.1.5.1", Substitute(x = PowerSeries(1 + poly_x_y), pso_pol), subs(x = poly_x_y, Truncate(pso_pol_shift, 20)));
Try[verify,up_to_homdeg(20)]("5.2.1", Substitute({x = PowerSeries(poly_x), y = PowerSeries(poly_x_y)}, pso_rat),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_rat, 20)));
Try[verify,up_to_homdeg(20)]("5.2.8", Substitute({x = PowerSeries(poly_x), y = 1 + poly_x_y}, pso_pol),
                             subs({x = poly_x, y = poly_x_y}, Truncate(pso_pol_shift_y, 20)));
Try[verify,up_to_homdeg(20)]("5.3.3", Substitute({x = 0, y = PowerSeries(poly_x)}, pso_exp),
                             subs({x = 0, y = poly_x}, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("5.4.0", Substitute(x = PowerSeries(y), pso_exp), subs(x = y, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("5.4.3", Substitute({x = y, y = PowerSeries(poly_x)}, pso_exp),
                             subs({x = y, y = poly_x}, Truncate(pso_exp, 20)));
Try[verify,up_to_homdeg(20)]("5.4.6", Substitute({x = PowerSeries(y), y = PowerSeries(1 + poly_x)}, pso_exp),
                             subs({x = y, y = poly_x}, Truncate(pso_exp_shift, 20)));

# Substitute non-polynomial, non-unit power series
## - single substitution
Try[verify,up_to_homdeg(8)]("6.0.0", Substitute(x = pso_exp - 1, pso_exp), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_exp, 8)));
Try[verify,up_to_homdeg(8)]("6.0.1", Substitute(x = pso_exp - 1, pso_rat), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_rat, 8)));
Try[verify,up_to_homdeg(8)]("6.0.2", Substitute(x = pso_exp - 1, pso_pol), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_pol, 8)));
Try[verify,up_to_homdeg(8)]("6.0.3", Substitute(x = pso_rat - 1, pso_exp), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_exp, 8)));
Try[verify,up_to_homdeg(8)]("6.0.4", Substitute(x = pso_rat - 1, pso_rat), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_rat, 8)));
Try[verify,up_to_homdeg(8)]("6.0.5", Substitute(x = pso_rat - 1, pso_pol), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_pol, 8)));

## - add a poly substitution
Try[verify,up_to_homdeg(8)]("6.1.0", Substitute({x = pso_exp - 1, y = poly_x_y}, pso_exp),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_exp, 8)));
Try[verify,up_to_homdeg(8)]("6.1.1", Substitute({x = pso_exp - 1, y = poly_x_y}, pso_rat),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_rat, 8)));
Try[verify,up_to_homdeg(8)]("6.1.2", Substitute({x = pso_exp - 1, y = poly_x_y}, pso_pol),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_pol, 8)));
Try[verify,up_to_homdeg(8)]("6.1.3", Substitute({x = pso_rat - 1, y = poly_x_y}, pso_exp),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_exp, 8)));
Try[verify,up_to_homdeg(8)]("6.1.4", Substitute({x = pso_rat - 1, y = poly_x_y}, pso_rat),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_rat, 8)));
Try[verify,up_to_homdeg(8)]("6.1.5", Substitute({x = pso_rat - 1, y = poly_x_y}, pso_pol),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_pol, 8)));

## - two non-poly, non-unit power series substitutions
Try[verify,up_to_homdeg(8)]("6.2.0", Substitute({x = pso_exp - 1, y = pso_rat - 1}, pso_exp),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_exp, 8)));
Try[verify,up_to_homdeg(8)]("6.2.1", Substitute({x = pso_exp - 1, y = pso_rat - 1}, pso_rat),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_rat, 8)));
Try[verify,up_to_homdeg(8)]("6.2.2", Substitute({x = pso_exp - 1, y = pso_rat - 1}, pso_pol),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_pol, 8)));

# Substitute non-polynomial, unit power series
## - single substitution
Try[verify,up_to_homdeg(8)]("7.0.0", Substitute(x = pso_exp, pso_exp), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_exp_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.0.1", Substitute(x = pso_exp, pso_rat), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_rat_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.0.2", Substitute(x = pso_exp, pso_pol), subs(x = Truncate(pso_exp, 8) - 1, Truncate(pso_pol_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.0.3", Substitute(x = pso_rat, pso_exp), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_exp_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.0.4", Substitute(x = pso_rat, pso_rat), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_rat_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.0.5", Substitute(x = pso_rat, pso_pol), subs(x = Truncate(pso_rat, 8) - 1, Truncate(pso_pol_shift, 8)));
 
## - add a poly substitution
Try[verify,up_to_homdeg(8)]("7.1.0", Substitute({x = pso_exp, y = poly_x_y}, pso_exp),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_exp_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.1.1", Substitute({x = pso_exp, y = poly_x_y}, pso_rat),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_rat_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.1.2", Substitute({x = pso_exp, y = poly_x_y}, pso_pol),
                            subs({x = Truncate(pso_exp, 8) - 1, y = poly_x_y}, Truncate(pso_pol_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.1.3", Substitute({x = pso_rat, y = poly_x_y}, pso_exp),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_exp_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.1.4", Substitute({x = pso_rat, y = poly_x_y}, pso_rat),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_rat_shift, 8)));
Try[verify,up_to_homdeg(8)]("7.1.5", Substitute({x = pso_rat, y = poly_x_y}, pso_pol),
                            subs({x = Truncate(pso_rat, 8) - 1, y = poly_x_y}, Truncate(pso_pol_shift, 8)));
 
## - two non-poly, non-unit power series substitutions
Try[verify,up_to_homdeg(8)]("7.2.0", Substitute({x = pso_exp, y = pso_rat}, pso_exp),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_exp_shift_both, 8)));
Try[verify,up_to_homdeg(8)]("7.2.1", Substitute({x = pso_exp, y = pso_rat}, pso_rat),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_rat_shift_both, 8)));
Try[verify,up_to_homdeg(8)]("7.2.2", Substitute({x = pso_exp, y = pso_rat}, pso_pol),
                            subs({x = Truncate(pso_exp, 8) - 1, y = Truncate(pso_rat, 8) - 1}, Truncate(pso_pol_shift_both, 8)));

# The 'product of powers' case

## - non-error examples
pso_in := 1/PowerSeries(1 + x*y + x*z);

eqns := {y=z^2/x, z=y^2};
Try[verify,up_to_homdeg(20)]("10.0.0", Substitute(eqns, pso_in), subs(eqns, Truncate(pso_in, 20)));
eqns := {y=z^3/x^2, z=y, x=x^2};
Try[verify,up_to_homdeg(20)]("10.0.1", Substitute(eqns, pso_in), subs(eqns, Truncate(pso_in, 20)));

## - error examples
eqns := {y=z/x};
Try[testerror]("10.1.0", Substitute(eqns, pso_in), "unexpected substitution with nonpositive total degree: %1", y=z/x);
eqns := {y=x^2/z};
Try[type]("10.1.1", Substitute(eqns, pso_in), 'object', 'assign' = 'pso_out');
Try[testerror]("10.1.2", Truncate(pso_out, 4), "non-polynomial result for power series: %1", -x^3/z - x*z);

#end test
