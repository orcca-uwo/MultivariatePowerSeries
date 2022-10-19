#test

with(TestTools):

VerifyTools:-AddVerification(
    ps_against_series = proc(x :: object,
                             y :: algebraic,
                             o := 10,
                             series_extra := 3,
                             $)
    local x_trunc := Truncate(x, o);

    local yy := y;
        for local v in indets(y, ':-name') minus {:-constants} do
            yy := series(yy, v, o + series_extra);
            yy := convert(yy, ':-polynom');
        end do;
    local terms := convert(expand(yy), ':-list', ':-`+`');
        terms := select(t -> degree(t) <= o, terms);

        return evalb(x_trunc = add(terms));
    end proc);

Try[type]("0.0.0", MultivariatePowerSeries:-PowerSeries(x + y + z + w^6), 'object', 'assign' = 'ps1');
Try[testnoerror]("0.0.1", x + y + z + w^6, 'assign' = 'expr1');
Try[type]("0.1.0", MultivariatePowerSeries:-PowerSeries(d -> ifelse(d :: even, 0, (-1)^((d - 1)/2) * x^d / d!), analytic = sin(x)),
          'object', 'assign' = 'ps2');
Try[testnoerror]("0.1.1", sin(x), 'assign' = 'expr2');
Try[type]("0.2.0", MultivariatePowerSeries:-PowerSeries(d -> ifelse(d :: even, 0, (-1)^((d - 1)/2) * x^d / d!), variables = {x}),
          'object', 'assign' = 'ps3');
Try[type]("0.3.0", MultivariatePowerSeries:-PowerSeries(d -> y^d / d!, analytic = exp(y)),
          'object', 'assign' = 'ps4');
Try[testnoerror]("0.3.1", exp(y), 'assign' = 'expr4');
Try[type]("0.4.0", ps2 * ps4, 'object', 'assign' = 'ps5');
Try[testnoerror]("0.4.1", exp(y) * sin(x), 'assign' = 'expr5');
Try[type]("0.5.0", MultivariatePowerSeries:-PowerSeries(1/(x^2 + y - 4)), 'object', 'assign' = 'ps6');
Try[testnoerror]("0.5.1", 1/(x^2 + y - 4), 'assign' = 'expr6');
Try[type]("0.6.0", MultivariatePowerSeries:-PowerSeries(d -> y^d / d!, analytic = 'int'(exp(z), z=-infinity .. y)), 'object', 'assign' = 'ps7');


# Cases that should work:
Try[verify,ps_against_series]("1.0", TaylorShift(ps1, w = 1), subs(w=w+1, expr1));
Try[verify,ps_against_series]("1.1", TaylorShift(ps1, w = sqrt(2)), subs(w=w+sqrt(2), expr1));
Try[verify,ps_against_series]("1.2", TaylorShift(ps2, x = -1), subs(x=x-1, expr2));
Try[verify,ps_against_series]("1.3", TaylorShift(ps5, y = ln(2)), subs(y = y+ln(2), expr5));

# Cases that should error out:
Try[testerror]("2.0", TaylorShift(ps1, [x = 1, x = 2, y = 3]), "these variables occurred multiple times on the left hand side: {x}");
Try[testerror]("2.1", TaylorShift(ps3, x = -1), "cannot compute the Taylor shift of a power series of which the analytic "
               "expression is not known");
Try[testerror]("2.2", TaylorShift(ps6, x = -2), "invalid Taylor shift: tried to shift to a pole of the analytic expression");
Try[testerror]("2.3", TaylorShift(1/(ps4 - exp(1)), y = 1), "invalid Taylor shift: tried to shift to a pole of the analytic expression");

# Test that some trivial transformations yield the same object.
Try("3.0", TaylorShift(ps1, this_variable_does_not_occur = 10), ps1);
Try("3.1", TaylorShift(ps1, w = 0), ps1);
Try("3.2", TaylorShift(ps7, z = 2), ps7);

# Test that we properly deal with multiple shifts
Try[verify,ps_against_series]("4.0", TaylorShift(ps6, [x = 1, y = 2]), subs(x=x+1, y=y+2, expr6));
Try[verify,ps_against_series]("4.1", TaylorShift(ps5, [x = Pi/2, y = ln(2)]), subs(x=x+Pi/2, y=y+ln(2), expr5));

#end test
