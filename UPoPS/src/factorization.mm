##########################################
## # Factorization via Hensel's' Lemma UPoPS
##########################################

# EvaluateAtOrigin
# given a up :: UPoPS in K[[x_1, ..., x_n]][Y] and a name ('Y' is default)
# return univariate polynomial up(0, ..., 0)[Y] in K[Y].
export 
    EvaluateAtOrigin ::static := proc(up :: UnivariatePolynomialOverPowerSeriesObject, $)
        # _Z shouldn't be documented 
        local x := ifelse(up:-vname <> undefined, up:-vname, '_Z');
        local i;
        return add(seq(PowerSeriesObject:-HomogeneousPart(up:-upoly[i], 0)*x^i, i = 0 .. DEGREE(up)));
    end proc;

# Convert a constant to RootOf form. This is for coefficients such as cos(Pi/12).
local
    convert_algebraic :: static := proc(expr :: constant, $)
    local result := convert(expr, ':-RootOf');
        if not type(result, ':-algnum') then
            error "invalid input: expected a univariate polynomial over power series with "
            "coefficients (after EvaluateAtOrigin) that are algebraic numbers, but found "
            "coefficient involving %1%2", expr,
            ifelse(expr = result, "", sprintf(" (= %a)", result));
        end if;

        return result;
    end proc;

local
    make_monic :: static := proc(up :: UnivariatePolynomialOverPowerSeriesObject,
                                 {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                  := ':-automatic'},
                                 $)
    local lc := up:-upoly[DEGREE(up)];
    local lc_constant_part := lc:-HomogeneousPart(lc, 0);
        if lc_constant_part = 0 then
            error "invalid input: expected univariate polynomial over power series with "
            "invertible leading coefficient, but the leading coefficient has constant part 0";
        end if;

        if lc_constant_part = 1 and lc:-GetAnalyticExpression(lc) = 1 then
            return up, (if returnleadingcoefficient = true then
                            [UnivariatePolynomialOverPowerSeriesObject(1, up:-vname)];
                        else
                            [];
                        end if);
        else
            # Make scaled equal to up / lc.
            local scaled := up:-DeepCopy(up);
            local lc_inverse := lc:-Inverse(lc);
            for local i from 0 to DEGREE(up) - 1 do
                scaled:-upoly[i] := lc_inverse:-BinaryMultiply(lc_inverse, scaled:-upoly[i]);
            end do;
            # We know the leading coefficient is exactly equal to 1, now.
            scaled:-upoly[DEGREE(scaled)] := PowerSeriesObject:-One();
            return scaled, (if returnleadingcoefficient = false then
                                [];
                            else
                                [UnivariatePolynomialOverPowerSeriesObject([lc], up:-vname)];
                            end if);
        end if;
    end proc;

# HenselFactorization
# Factorization up :: UPoPS via Hensel's' lemma w.r.t the main variable (up:vname)
export 
    HenselFactorize :: static := proc(_up :: UnivariatePolynomialOverPowerSeriesObject,
                                      {returnleadingcoefficient :: {truefalse, identical(automatic)}
                                       := ':-automatic'},
                                      $)
        local up, extra_factors;
        up, extra_factors := make_monic(_up, _options['returnleadingcoefficient']);

        if DEGREE(up) = 0 then
            return extra_factors;
        end if;

        if up:-vname = undefined then
            error "invalid input: the univariate polynomial over power series is non-constant, but its "
            "main variable is unspecified";
        end if;
        
        # compute roots of the univariate polynomial up when x_1 -> 0, ..., x_n -> 0.
        local up_at_origin := subsindets[':-flat'](
            EvaluateAtOrigin(up), ':-And'(':-constant', ':-Not'(':-complex'(':-numeric')), ':-Not'(':-algnum')),
            convert_algebraic);
        local c := SolveTools:-Polynomial(up_at_origin, up:-vname, ':-dropmultiplicity');
        local r := numelems(c);
        if r = 1 then
            return [op(extra_factors), DeepCopy(up)];
        end if;
        local F := Array(1 .. r); # factors of up
        local f := up;
        for local i from 1 to r do
            local g := TaylorShift(f, c[i]); # g := f(Y + c_i) 
            local p, alpha;
            p, alpha := WeierstrassPreparation(g);
            F[i] := TaylorShift(p, -c[i]); # p(Y - c_i)
            if i+1 < r then 
                # shifting back to alpha(Y - c_i) to generate F[i+1], ..., F[r] in the next iterations
                f := TaylorShift(alpha, -c[i]);
            end if;
        end do;
        return [op(extra_factors), seq(F)];
    end proc;
