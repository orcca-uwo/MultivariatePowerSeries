# We apply the transformation x-(1/n)*a(y) for a univariate 
# polynomial of degree n in the variable x, with coefficients 
# Puiseux series in the variable y. 
# If the option detransform is specified, then the transformation
# x+(1/n)*a(y) is applied.
export TschirnhausenTransformation::static := proc(up::UnivariatePolynomialOverPowerSeriesObject, 
                                                   a_in::{PuiseuxSeriesObject, PowerSeriesObject},
                                                   m_in::nonnegint := 0,
                                                   {detransform :: truefalse := false}, 
                                                   {powerseriesmode::truefalse:=false}, $)
    local l, k, coeff, uno, m;
    local n := DEGREE(up);

    if m_in = 0 then 
        m := n;
    else 
        m := m_in;
    end if;
    
    if detransform= true then 
        uno := 1;
    else 
        uno := -1;
    end if;

    local new_up := up;

    if n < 1 then 
        return up;
    end if;

    local a := a_in;              

    # testeq to test if the expression is zero.
    if testeq(a:-GetAnalyticExpression(a)) then
        return up;
    end if;

    local A := Array(0..n);

    ## We have a polynomial of the form
    ## p:= x^n+a_(n-1)(y)*x^(n-1)+...+ a_1(y)*x^1.
    ## we want to delete the monomial of degree n-1. So
    ## we apply the change of variables given by:
    ## x-(1/n)*a_(n-1)(y).
    ## Then, we expand (x-(1/n)*a_(n-1)(y))^i for i=1,..,n-1
    ## using the binomial theorem, and reorganize the sums.
    ## Let a:= a_{n-1}, then we get: 
    ## Sum_{l=0}^n Sum_{k=0}^l (-1)^k*a_{n-l+k}*binomial(n-l+k, k)*y^(n-l)((1/n)*a)^k.
    local exp_a := Array(0..n);

    if powerseriesmode then 
        exp_a[0] := PowerSeries(1);
    else 
        exp_a[0] := PuiseuxSeries(1);
    end if;

    for k to n do
      exp_a[k] := exp_a[k-1] * a;
    end do;

    for l from 1 to n do
        for k from 0 to l do
            coeff := up:-upoly[n-l+k];

            A[l] := A[l] + (uno/m)^k*coeff*binomial(n-l+k, k)
                            *exp_a[k];
        end do;
    end do;

    # Case l=0.
    A[0] := up:-upoly[n];

    ArrayTools:-Reverse(A, ':-inplace');

    return Object(UnivariatePolynomialOverPowerSeriesObject, A, up:-MainVariable(up));
end proc; 

# To set a nonzero_pso_bound or a nonzero_pso_bound_static.
export SetPuiseuxBound ::static := proc(_self :: UnivariatePolynomialOverPowerSeriesObject,
                                                 bound::nonnegint,
                                                 {instance :: truefalse := true},
                                                 $)
    local old_value;
    
    if instance then
        old_value := _self:-Puiseux_bound;
        _self:-Puiseux_bound := bound;
    else
        old_value := _self:-Puiseux_bound_static;
        _self:-Puiseux_bound_static := bound;
    end if;

    return old_value;
end proc;

local ConvertBackToUPoPuS::static := proc(_self::UnivariatePolynomialOverPowerSeriesObject, 
                                            r::rational, q::integer, 
                                            var::name, vars::list(name), $)
    
    local upoly := _self:-upoly;
    local pso;
    local vars_as_set := convert(vars, ':-set');
    local A := Array(0..(numelems(upoly)-1));

    for local i from 0 to numelems(upoly)-1 do 
        if upoly[i]:-Variables(upoly[i])={} then 
            local c := upoly[i]:-GetAnalyticExpression(upoly[i]);
            if c <> undefined then 
                pso := PowerSeries(d -> ifelse(d=0, c, 0), variables=vars_as_set, analytic=c); 
                A[i] :=PuiseuxSeries(pso, vars, vars, [[1/q]], [vars[1]=(-r)*i]);
            else 
                A[i] := PuiseuxSeries(upoly[i], vars, [vars[1]=(-r)*i]);
            end if;  
        else
            A[i] :=PuiseuxSeries(upoly[i], vars, vars, [[1/q]], [vars[1]=(-r)*i]);
        end if;       
    end do;

    return UnivariatePolynomialOverPowerSeries(A, _self:-vname);
end proc;

local FactorVariable::static := proc(up::UnivariatePolynomialOverPowerSeriesObject,
                                     bnd::nonnegint:= FAIL, $)
    local one_upops := Object(UnivariatePolynomialOverPowerSeriesObject, 
                             [1]);

    local n := DEGREE(up);
    local a := up:-upoly[0];  

    local i :=1;
    while i<=n do
        a := up:-upoly[i];

        local check := testeq(a:-GetAnalyticExpression(a));
        if check then
            i:=i+1;
        elif check=FAIL or a:-GetAnalyticExpression(a)=undefined then 
            local appro := a:-Truncate(a, bnd);

            if appro=0 then 
                error "invalid input: constant term in the main variable of"
                      " %1 may be equal to 0 which is not allowed", _up;
            else 
                break;
            end if;
        else    
            break;
        end if; 
    end do;

    local upoly := Array(0..n-i);
    upoly[0..n-i]:= up:-upoly[i..n];

    local new_up := Object(UnivariatePolynomialOverPowerSeriesObject,
                           upoly, up:-vname); 

    return [i, new_up];

end proc;                                    

# Nowak's formulation of the Newton-Puiseux theorem.
export PuiseuxTheorem::static := proc(_up::UnivariatePolynomialOverPowerSeriesObject, 
                                      bnd::nonnegint := FAIL, 
                                      {returnleadingcoefficient::{truefalse, identical(automatic)}
                                        := ':-automatic'}, 
                                      {useevala::truefalse := false},
                                      $)
    local k, s, my_Puiseux_bound, i, ld_coeff, Puiseux_fac;
    local facs := [0,1];
    local var := _up:-vname;
    local var_coeff := op(_up:-Variables(_up) minus {var});
    local up := _up:-ConvertToPuSOUPoP(_up);
    local a := up:-upoly[0]; 

    # We check the Puiseux's bound.
    if bnd <> FAIL then
        my_Puiseux_bound := bnd;
    elif _up:-Puiseux_bound <> undefined then
        my_Puiseux_bound := _up:-Puiseux_bound;
    else
        my_Puiseux_bound := _up:-Puiseux_bound_static;
    end if;

    # up must be monic, the monomial of degree n-1 must be zero
    # and the constant term must be different than zero.
    # We check that the constant term is non-zero.
    local check := testeq(a:-GetAnalyticExpression(a));
    if check then
        facs := up:-FactorVariable(up, my_Puiseux_bound);
        up := facs[2];

    elif check=FAIL then 
        local appro := a:-Truncate(a, my_Puiseux_bound);
        if appro=0 then 
            error "invalid input: constant term in the main variable of"
              " %1 may be equal to 0 which is not allowed", _up;
        end if;
    end if;

    local n := DEGREE(up);
    # If up is a linear polynomial, there is nothing to do.
    if n <= 1 then 
        return up;
    end if;

    # We make up monic, and save its leading coefficient if required.
    if returnleadingcoefficient = true then 
        local dummy := up:-MakeMonic(up, ':-returnleadingcoefficient'=true);
        up := dummy[1];
        ld_coeff := dummy[2];
    else 
        up := up:-MakeMonic(up, ':-returnleadingcoefficient'=false);  
    end if;

    local c := up:-upoly[n-1];

    # We apply the Tschirnhausen Transformation to delete the 
    # monomial of degree n-1. 
    up := up:-TschirnhausenTransformation(up, c);

    ASSERT(up:-upoly[n-1]:-ApproximatelyZero(up:-upoly[n-1], my_Puiseux_bound));

    # We compute the order of the coefficients as Puiseux series.
    # We know that up:-upoly[n-1] is 0 by construction. So, we do not need to
    # send it to the GetOrder function.
    # Note: that if analytic expression of up:-upoly[n-1] is undefined,
    # GetOrder would always throw an error.
    # We also know that up is monic, so the order of up:-upoly[n] is 0.
    local ords := map(PuiseuxSeriesObject:-GetOrder, up:-upoly[0..n-2], my_Puiseux_bound);
    # We divide the orders by k. 
    # Note: We do not need to include orders of up:-upoly[n-1] nor up:-upoly[n].
    local ords_over_k := [seq(ords[k]/(n-k), k=0..n-2)];
    # We get the minimum.
    local r := min(ords_over_k);

    local rays := [seq(op(up:-upoly[k]:-GetRays(up:-upoly[k])), k=0..n-1)];
    # TODO: Think how to compute q.
    local q := ilcm(seq(seq(denom(s), s in i), i in rays), denom(r));

    # We compute p such that r:= p/q.
    local p := r*q;
    local ONE := PowerSeriesObject:-One();

    # We apply the change of variables var_coeff=w^q, var=U*w^p.
    # Then up = w^{n*p}*Q(w,U), where Q is equal to 
    # Q :=  U^n +b_2(w)*U^{n−2} +· · ·+b_n(w) with
    # b_k(w) = a_k(w^q)*w^{−kp}. 
    local upoly := [seq(PuiseuxSeriesObject:-ChangeOfVariablesForPuiseuxTheorem(
                            up:-upoly[k], q, p, n-k), k=0..n-1), ONE];

    # The Puiseux theorem (Nowak's version) guarantees that 
    # the coefficients in upoly are all PSO. Thus, we create 
    # a UPoPS with PSO coefficients.
    local upops := Object(UnivariatePolynomialOverPowerSeriesObject, upoly, var);

    # We apply the Hensel lemma to upops.
    local  Hensel_fac := upops:-HenselFactorize(upops, _options['useevala']);
    local l;
    if add(DEGREE(l), l in Hensel_fac)<>DEGREE(upops) then
        error "the PuiseuxTheorem command failed while applying"
                " the HenselFactorize algorithm."
                " Try rerunning the command with the"
                " option useevala=true.";
    end if;

    # We convert back to polynomial over Puiseux series.
    # We apply the change of variables w = var_coeff^(1/q), 
    # U=var_coeff^(-r)*var to our factors.
    Hensel_fac := map(UnivariatePolynomialOverPowerSeriesObject:-ConvertBackToUPoPuS, 
                       Hensel_fac, r, q, var, [var_coeff]);
    local extra := Object(PuiseuxSeriesObject, ONE, [var_coeff], [var_coeff=n*r]);

    local extra_fac := Object(UnivariatePolynomialOverPowerSeriesObject,
                                     [extra], var);                                

    Hensel_fac := [extra_fac, op(Hensel_fac)];

    if returnleadingcoefficient<>false then 
        local extra_ld := Object(PuiseuxSeriesObject, ONE, [var_coeff]);

        if returnleadingcoefficient=true then 
            extra_ld := extra_ld*ld_coeff;
        end if;

        local extra_fac_ld := Object(UnivariatePolynomialOverPowerSeriesObject,
                                     [extra_ld], var);                                

        Puiseux_fac := [extra_fac_ld, op(Hensel_fac)];
    else 
        Puiseux_fac :=Hensel_fac;
    end;

    # We revert the TschirnhausenTransformation.
    Puiseux_fac := map(UnivariatePolynomialOverPowerSeriesObject:-TschirnhausenTransformation,
                      Puiseux_fac, c, n, detransform=true);

    # We append facs[1] times var as a UPOPS.
    if facs[1]>=1 then
        local up_var := Object(UnivariatePolynomialOverPowerSeriesObject, 
                            [PowerSeriesObject:-Zero(),ONE], 'var');
        Puiseux_fac := [up_var$facs[1], op(Puiseux_fac)];
    end if;
    
    return Puiseux_fac;
end proc; 
