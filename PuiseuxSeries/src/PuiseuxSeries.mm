module PuiseuxSeriesObject()
option object;

local 
    pso, # Power Series object
    rays::list(list(rational)), #list containing the information
                      # of the exponent of each variable
    e::table, # Exponent of the monomial multiplying the lso
    ord::list, # list containing the order of the Puiseux series variable
    ordCV::list, # list containing the order of the power series variable
    dstyle::PSDS_TYPE, #a display style for this object
    # Degree bounds for the inverse function
    nonzero_pso_bound_static :: static(nonnegint) := 10,
    nonzero_pso_bound :: {nonnegint, identical(undefined)} := undefined,
    smallest_term_bound_static :: static(nonnegint) := 10,
    smallest_term_bound :: {nonnegint, identical(undefined)} := undefined,
    # Bound for the Puiseux theorem
    Puiseux_bound_static :: static(nonnegint) := 10,
    Puiseux_bound :: {nonnegint, identical(undefined)} := undefined; 


local ModuleApply::static := proc()
    return Object(PuiseuxSeriesObject, _passed);
end proc;   

# Constructor of PuSO using rays.
local CreatePuSOFromRays::static := proc(new :: PuiseuxSeriesObject, 
                                 old :: PuiseuxSeriesObject,
                                 pso::{PowerSeriesObject, polynom},
                                 ord::list(name) :=[],
                                 ordCV::list(name) := [],
                                 rays::list(list(rational)) := [],
                                 e::{table, list(name=rational)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $) option overload(':-callseq_only');
    local l, i, j, eq;
    
    if type(pso, ':-polynom') then
        new:-pso := MultivariatePowerSeries:-PowerSeries(pso);
    else
        new:-pso := pso;
    end if;
    new:-dstyle := dstyle;

    if not ord = [] then
        new:-ord := ord;
    else 
        new:-ord := convert(new:-pso:-Variables(new:-pso), ':-list');
        if numelems(new:-ord) > 1 and rays <> [] then
            error "you cannot specify rays if you don't specify the order of the power series variables";
        end if;
    end if;

    # We check if there is an order for the power series variable.
    # If there is no order, then we use the one in pso
    if ordCV = [] then
        if nops(new:-pso:-Variables(new:-pso))<=1 then
            new:-ordCV := [op(new:-pso:-Variables(new:-pso))]; 
        elif ord = [] then
            new:-ordCV := new:-ord;
        else 
            error "an order for the power series variables must be specified";
        end if;
    elif new:-pso:-Variables(new:-pso) = convert(ordCV, ':-set') then
        new:-ordCV := ordCV;
    else
        error "the variables in %1 must be equal to %2", new:-pso, ordCV;
    end if;

    # We check if there are rays for the change of variables
    if rays = [] then
        if numelems(new:-ordCV)=numelems(new:-ord) then
            # If there is a lorder of the power series 
            # variable and no rays, we assume that this are 
            # the rays (1,0,...,0), (0,1,0,...,0),...,(0,...,0,1).
            new:-rays := [seq( [seq(ifelse(j=i, 1, 0), j=1..nops(new:-ord))] ,i=1..nops(new:-ord))];
        elif numelems(new:-ordCV) = 0 then
            new:-rays := [];
        else
            error "the rays must be specified if the number of variables in the Puiseux series variable and the power series variable is different";
        end if;
    elif numelems(new:-ordCV) <> numelems(rays)  then
        error "number of rays must agree with the number of variables in %1 or be empty", new:-pso;
    else
        # We check that the rays have the correct dimension for the 
        # change of variable
        for i in rays do 
            if numelems(i)<>numelems(new:-ord) then
                error "all the rays must have dimension %1", numelems(new:-ord);
            end if;
        end do;

        # We check that the rays are grevlex positive with the given order
        if new:-CheckRays(new, rays) then
            new:-rays := rays;
        else 
            error "all the rays in %1 must be grevlex(%2) positive", rays, new:-ord;
        end if;
    end if;

    # We define the exponents for the monomial that multiplies lso
    new:-e := new:-MakeMonomialTable(new, e);

end proc;

# Constructor of PuSO based in equations. No orders needed.
local CreatePuSOFromEqu::static :=  proc(new :: PuiseuxSeriesObject, 
                                 old :: PuiseuxSeriesObject,
                                 pso::{PowerSeriesObject, polynom},
                                 mp::list(`=`( name, {:-`*`( { name, name^rational } ), name, name^rational})) := [],
                                 e::{table, list(name=rational)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $) option overload(':-callseq_only');
    
    if type(pso, ':-polynom') then
        new:-pso := MultivariatePowerSeries:-PowerSeries(pso);
    else
        new:-pso := pso;
    end if;

    new:-dstyle := dstyle;
    new:-ordCV := map(lhs, mp);

    # We check that there is an equation for each variable in pso
    if not convert(new:-ordCV, ':-set') = new:-pso:-Variables(new:-pso) then
        error "the variables in %1 must be equal to %2", mp, new:-pso:-Variables(new:-pso);
    end if;

    # We convert the rhs of the equations as a list separated by *
    local as_products := map(convert, map(rhs, mp), ':-list', :-`*`);
    # We separate the previous list as a list separated by ^
    local as_powers := map2(map, convert, as_products, ':-list', `^`);
    # We generate the order for the cone
    local ord := map2(op, 1, map(op, as_powers));
    new:-ord := convert(convert(ord, ':-set'), ':-list');

    # We generate the rays based in the order and as_powers
    new:-rays := new:-GenerateRays(new, as_powers);

    # We check that the rays are grevlex positive with the given order
    if not new:-CheckRays(new, new:-rays) then
        error "all the rays in %1 must be grevlex(%2) positive", new:-rays, new:-ord;
    end if;

    # We define the exponents for the monomial that multiplies lso
    new:-e := new:-MakeMonomialTable(new, e);

end proc;

# Constructor of PuSO based in equations and a rational polynomic function. 
# No orders needed.
local CreatePuSOFromRatFunc::static :=  proc(new :: PuiseuxSeriesObject, 
                                 old :: PuiseuxSeriesObject,
                                 f::ratpoly,
                                 mp::list('`=`'( name, {':-`*`'( { name, name^rational } ), name, name^rational})) := undefined,
                                 e::{table, list(name=rational)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $) option overload(':-callseq_only');
    
    local g, e1, h, g_reduced, q, const, V, v, 
          rays_q_CV, r, s, eq, sm;

    new:-dstyle := dstyle;

    local new_mp := mp;

    # We check if there is a list of equations as input.
    if mp = undefined then
        new_mp := [seq(v=v, v in (indets( f, ':-And'(':-name',':-Not'(':-constant')))) )];

    # We check that there is an equation for each variable in f
    elif not convert(map(lhs, new_mp), ':-set') = (indets( f, ':-And'(':-name',':-Not'(':-constant')))) then
        error "the variables in %1 must be equal to the variables in %2", new_mp, f;
    end if;

    # We convert the rhs of the equations as a list separated by ^
    local as_powers :=  new:-ComputeAsPowers(new_mp);
    # We generate the order for the cone
    local ord := map2(op, 1, map(op, as_powers));
    new:-ord := convert(convert(ord, ':-set'), ':-list');

    # We generate the rays based in the order and as_powers.
    # Note that this rays could change for two possible reasons:
    # We simplify one of the variables in f.
    # We do not have a intervible polynomial in the denominator.
    local rays := new:-GenerateRays(new, as_powers);
    # We make the rays integer and save the lcm of the
    # denominators.
    local denominator := lcm(seq(seq(denom(s), s in r) , r in rays));
    
    #########
    # Puiseux Series adjustment.
    if denominator<>1 then
        rays := denominator*rays;
    end if;
    #########

    # We check that the rays are grevlex positive with the given order
    if not new:-CheckRays(new, rays) then
        error "all the rays in %1 must be grevlex(%2) positive", rays, new:-ord;
    end if;

    # We Power each equation to denominator to get integer exponents.
    local local_mp := [seq(lhs(eq)=rhs(eq)^denominator, eq in new_mp)];
    local gcd_h := 1;
    local sm_ray := [];

    # We write f as g_reduced/(gcd_h*q), i.e.,
    # we simplify the fraction and factorize the denominator.
    # Note that we need to use gcd_h in the computation of e.
    g_reduced, q, gcd_h := PuiseuxSeriesObject:-SimplifyFraction(PuiseuxSeriesObject, f);
                                                       
    # Variables of q. Note that this variables could be different 
    # from the variables of f. 
    #V := convert( indets( q, ':-name' ) minus {':-constants'}, ':-list' );
    V := convert( indets( q, ':-And'(':-name',':-Not'(':-constant')) ), ':-list' );
    # Constant term of q.
    const := eval(q, V =~ 0);

    if not type(const, ':-numeric') then
        const := Algebraic:-Expand(subsindets(d, ':-radical', convert, RootOf));
    end if;

    # We check if some variable was eliminated during the 
    # simplification an recompute the rays and local_mp.
    local_mp, rays := new:-CheckVariableSimplification(new, f, 
                                                       g_reduced/q,
                                                       rays,
                                                       local_mp);

    # We compute a local ordCV
    local local_ordCV := map(lhs, local_mp);                                                   

    if const <> 0 then
        # If q <> 0, then we can invert q as a Power Series.
        # We generate an order for the PSO.
        new:-ordCV := local_ordCV;

        # We generate the rays
        new:-rays := rays;
    else
        # We need to compute an appropiate change of variables.
        # We apply the change of variables given by mp.
        local q_CV := subs(local_mp, q);
        local q_CV_as_list := convert(q_CV, ':-list', ':-`+`');

        # We look for the smallest monomial in q_CV and save
        # this information as a ray. We also get all the rays
        # that represent the monomials in q_CV.
        # Note that sm_ray must be use in the computation of e.
        sm_ray, rays_q_CV := PuiseuxSeriesObject:-SmallestGrevlexMonomial(PuiseuxSeriesObject, 
                                                                    q_CV_as_list,
                                                                     new:-ord);
        # We generate a new set of rays, i.e., a new change of variables.
        local rays_minus_sm := map(':-`-`', rays_q_CV, sm_ray);
        new:-rays := new:-MakeRaysCompatible(new, [op(rays), op(rays_minus_sm), sm_ray], 
                                            new:-ord);

        local M := Matrix(new:-rays);
        # We compute new variables for our pso.
        new:-ordCV := new:-ComputeNewOrdCV(new, numelems(new:-rays), local_ordCV);

        # We compute the old variables in terms of the new ones.
        local change_of_variable := new:-PolynomialChangeOfVariables(new, M, new:-ordCV, 
                                                               map(lhs, local_mp), rays);

        # We compute the smallest monomial of q_CV in terms of the new variables.
        local sm_mon := new:-PolynomialChangeOfVariables(new, M, new:-ordCV, [sm], [sm_ray]);
        
        # We compute g_reduced, q and sm_mon after the new change of variables.
        g_reduced := subs(change_of_variable, g_reduced);
        q := subs(change_of_variable, q);
        sm_mon := rhs(sm_mon[1]);

        # We divide by the smallest term of q.
        if not divide(q, sm_mon, 'q') then
            error "while looking for the smallest grevlex monomial in %1", q;
        end if;
    end if;

    # We define the exponents for the monomial that multiplies lso
    new:-e := new:-MakeMonomialTable(new, e);

    # If gcd_h is different than 0, we need to devide e by it.
    if gcd_h <> 1 then
        # We apply the change of variable.
        local gcd_h_cv := subs(new_mp, gcd_h);
        # We see it as a list of the form [name, integer]
        local (constant_factors, gcd_h_as_nonconst_product) 
                := selectremove(type, convert(gcd_h_cv, ':-list', ':-`*`'), ':-constant');
        q := q * mul(constant_factors);
        
        local gcd_h_as_powers := map(convert, gcd_h_as_nonconst_product, ':-list', ':-`^`');

        # We subtract the previous information from the current e.
        for local p in gcd_h_as_powers do
            new:-e[p[1]] := new:-e[p[1]] - p[2]; 
        end do;
    end if;

    # We create a new PSO.
    new:-pso := MultivariatePowerSeries:-PowerSeries(g_reduced/q);

    # We divide sm_ray by denominator.
    sm_ray := (1/denominator)*sm_ray;

    # we need to simplify sm_rays from e.
    if sm_ray <> [] then
        for local i from 1 to numelems(new:-ord) do
            new:-e[new:-ord[i]] := new:-e[new:-ord[i]] - sm_ray[i]; 
        end do;
    end if;

    # We use our denominator to make the new:-rays rational again.
    #########
    # Puiseux Series adjustment.
    if denominator<>1 then
        new:-rays := (1/denominator)*new:-rays;
    end if;
    #########
end proc;
                            
export DeepCopy::static := proc(_self::PuiseuxSeriesObject, $)
    return Object(PuiseuxSeriesObject, _self:-pso, _self:-ord, 
                                       _self:-ordCV, _self:-rays, 
                                       _self:-e);
end proc;

export Variables::static := proc(_self::PuiseuxSeriesObject, $)
    return convert(_self:-ord, ':-set');
end proc;

local CreatePuSO_NoSuitableMethod::static := proc() option overload(':-callseq_only');
    error "invalid input: arguments to %1, %2, do not match any of the accepted calling sequences", procname, [ _passed ] ;
end proc;

# Constructor of the class.
local ModuleCopy::static := overload([CreatePuSOFromRays, CreatePuSOFromEqu, 
                                      CreatePuSOFromRatFunc,
                                      CreatePuSO_NoSuitableMethod]);


# Maybe we should at least attempt to do this in the case where one ray 
# is exactly in the direction of x (so w.l.o.g. u_1=x), and the others 
# are perpendicular to x (i.e. u_j is independent of x for j > 1). In 
# this case, we should really be able to figure this out, and I think 
# this is common enough. In that case we still need a bound for the degree 
# in x, which could come from the analytic expression or it could be given 
# by the user; let's say that bound is d. We'd create d+1 new PuiseuxSeries 
# objects P[0], ..., P[d], corresponding to the coefficients of x^0, ..., x^d. 
# Each would have the same rays as before except the one corresponding to x, 
# each ray projected to the subspace without x.

# So let's consider a term c * prod(u_i^a_i, i=1..k) with sum(a_i, i=1..k) = A. 
# This would contribute to the homogeneous part of P[a_1] of degree A - a_1. 
# In other words, we would just need to write a generator for PowerSeriesObject 
# that returns the coefficient of x^d for fixed x and d

# To convert a PuSO to a UPoP.
export ConvertToUPoPS ::static := proc(_self :: PuiseuxSeriesObject,
                                x :: name,
                                $)
    if not member(x, _self:-ord) then
        return Object(UnivariatePolynomialOverPowerSeriesObject, [_self], x);
    end if;
    
    local algexpr := _self:-GetAnalyticExpression(_self);

    ## TODO: We should be able to do this more general (a "polynomial" in x) 
    ## with: type(algexpr, ':-polynom'(':-anything', x)) .
    if algexpr = undefined or not type(algexpr, ':-polynom'(':-ratpoly', x)) then
        error "attempted to convert a Puiseux series involving %1 to a univariate polynomial over "
        "Puiseux series in %1, but it is not known to be polynomial in %1", x;
    end if;
    
    local i, d := degree(algexpr, x);   
    local coeffs := PolynomialTools:-CoefficientList(algexpr, x);

    coeffs := map(PuiseuxSeries, coeffs);

    return Object(UnivariatePolynomialOverPowerSeriesObject, coeffs, x);
    
end proc;                                     

# We check if some variable was eliminated during simplification.
# If so, we compute a new mp and rays without the information of
# the deleted variable.
local CheckVariableSimplification::static := proc(_self :: PuiseuxSeriesObject,
                                                  f::ratpoly,
                                                  new_f::ratpoly,
                                                  rays::list(list(rational)),
                                                  mp::list('`=`'( name, {':-`*`'( { name, name^integer } ), name, name^integer})) := [],
                                                  $)

    # We check if after simplifying f, we have less variables.
    local eq, as_powers;

    local A := indets( f, ':-name' ) minus {':-constants'};
    local B := indets( new_f, ':-name' ) minus {':-constants'};
    local C := A minus B;
    local local_mp := mp;
    local local_rays := rays;

    # If we have less variables than before, then we do not
    # take the deleted variables into account.
    if C <> {}  then
        local_mp := [seq(ifelse(lhs(eq) in C, NULL, eq), eq in mp)];
        # We generate the rays.
        as_powers := _self:-ComputeAsPowers(local_mp);
        local_rays := _self:-GenerateRays(_self, as_powers);
    end if;

    return local_mp, local_rays;
end proc; 

local ComputeAsPowers::static := proc(mp::list('`=`'( name, {':-`*`'( { name, name^rational } ), name, name^rational})) := [],
                                  $)
    # We convert the rhs of the equations as a list separated by *
    local as_products := map(convert, map(rhs, mp), ':-list', ':-`*`');
    # We separate the previous list as a list separated by ^
    local as_powers := map2(map, convert, as_products, ':-list', ':-`^`');

    return as_powers;
end proc;

# We compute the change of variables given by M. 
local PolynomialChangeOfVariables::static := proc(_self :: PuiseuxSeriesObject,
                                        M::Matrix,  
                                        new_ordCV::list(name), 
                                        old_ordCV::list(name), 
                                        rays::list(list(rational)),
                                        $) 
    
    local r;
    local rays_as_vectors := map(convert, rays, ':-Vector');

    # We do Linear Algebra to compute the change of variables in terms
    # of the new variables and rays.
    local result := [seq(LinearAlgebra:-LinearSolve(M^+, r), r in rays_as_vectors)];

    local cv := [seq(mul(new_ordCV^~r), r in result)];
    cv := old_ordCV=~cv;

    return cv;
end proc;


# We take a f := g/h, a rational polynomial.
# We simplify f and factor out the gcd of the terms of h.
local SimplifyFraction::static := proc(_self::PuiseuxSeriesObject,
                                       f::ratpoly,
                                       $)

    local h_reduced, g_reduced, q, i, c;
    local g := numer(f);
    local h := denom(f);

    # We simplfy g/h.
    gcd(g, h, 'g_reduced', 'h_reduced');
    
    local h_expanded := expand(h_reduced);

    # Terms of h_expanded as a list.
    local T_h := convert(h_expanded, ':-list', ':-`+`');

    # gcd of the terms in the denominator.
    local gcd_h := foldr( gcd, seq(T_h) );

    # We devide h_expanded over gcd_h
    evala(':-Divide'(op(subs(seq(:-constants[i] = c[i], i = 1 .. numelems([:-constants]))
            , [h_expanded, gcd_h])),'q'));

    return g_reduced, q, gcd_h;
end proc;


# Monomial generator.
local MakeMonomialTable::static := proc(_self::PuiseuxSeriesObject, e)
    local l;

    if type(e, ':-table') then
        l := [entries(e, ':-nolist', ':-pairs')];
        if not type(l, ':-list'(':-name'=':-rational')) then
            error "invalid input: %1 expects its %-2 argument, %3, to consist of entries of type %4, but received %5", 
                                    procname, 4, ':-e', ':-name' =':-rational', eval(e, 1);
        end if;
    else
        l := e;
    end if; 

    if convert(map(lhs,l), ':-set') subset convert(_self:-ord, ':-set') then
        return table(':-sparse' = 0, l);
    else
        error "%1 must be a list of equations in %2", eval(e, 1), _self:-ord;
    end if;
end proc;

# Ray generator.
# powlist represents a monomial a_1^(n_1)x...xa_z^(n_z).
# Using a list of the form powlist:=[[a_1,n_1], ..., [a_z,n_z]],
# we first verify that the variable a_i is in _self:-ord,
# if this happens, then we can save the exponent n_i in 
# our ray "result". 
local GenerateRay::static := proc(_self::PuiseuxSeriesObject, powlist, $)
    local n := numelems(_self:-ord);
    local result := Array(1 .. n);

    for local pair in powlist do
        local i := ListTools:-Search(pair[1], _self:-ord);
        result[i] := pair[2];
    end do;

    return convert(result, ':-list');
end proc;

# Rays generator.
local GenerateRays::static := proc(_self::PuiseuxSeriesObject, powlist, $)
    local ls;
    return [seq(_self:-GenerateRay(_self, ls    ), ls in powlist)];
end proc;

# Funtion that checks if the rays are grevlex positive with the given order or not
local CheckRays::static := proc(_self::PuiseuxSeriesObject, rays::list(list(rational)), $)
    local r;

    for r in rays do 
        if not _self:-Positive(_self, r) then
            return false;
        end if;
    end do;

    return true;
end proc; 

# Display the change of variables of a PuSO.
# We compute the change of variables given by the rays
# using the variables ord and ordCV. This is the change of 
# variable that we apply to our pso to get the PuSO.
export ChangeOfVariables::static := proc(_self :: PuiseuxSeriesObject, $)
    local r, i, j;

    r := [seq(_self:-ordCV[j]
            =mul(_self:-ord[i]^(_self:-rays[j][i]), i = 1 .. nops(_self:-ord)), 
                                                    j=1..nops(_self:-ordCV))];

    return r;
end proc;

# SetDisplayStyle.
# to set the dstyle of the input object.
export SetDisplayStyle :: static:= proc(_self :: PuiseuxSeriesObject, s :: PSDS_TYPE, $)
    _self:-dstyle := s; 
end proc;

export Display :: static := proc(_self :: PuiseuxSeriesObject,
                                 user_dstyle :: list := [],
                                 output :: identical("full", "typeset", "fullstring", "string")
                                 := "full",
                                 {updateterms :: name := NULL},
                                 $)
local mystyle;
    if user_dstyle <> [] then
        mystyle := user_dstyle;
    elif _self:-dstyle <> [] then
        mystyle := _self:-dstyle;
    else
        mystyle := PowerSeriesObject:-GetDefaultDisplayStyle();
    end if;

local x;
    return _self:-pso:-Display(
        _self:-pso, mystyle, output,
        ':-objectname' = "PuiseuxSeries", ':-substitutions' = _self:-ChangeOfVariables(_self),
        ':-cofactor' = mul(x^(_self:-e[x]), x in _self:-ord),
        ifelse(updateterms = NULL, NULL, _options['updateterms']));
end proc;


# ModulePrint: print _self :: PuiseuxSeriesObject.
local ModulePrint :: static := proc(_self :: PuiseuxSeriesObject, $)
    if IsWorksheetInterface() then
        return _self:-Display(_self);
    else
        return convert(_self:-Display(_self, [], "fullstring"), ':-name');
    end if;
end proc;

# Module Deconstruct
local ModuleDeconstruct :: static := proc(_self :: PuiseuxSeriesObject, $)
    return 'Object'('PuiseuxSeriesObject', _self:-pso, _self:-ord, _self:-ordCV, _self:-rays,
                    [indices(_self:-e, ':-pairs')]);
end proc;    

# Basic routines
$include "MultivariatePowerSeries/PuiseuxSeries/src/gets.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/grevlex.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/truncate.mm"
# Arithmetic 
$include "MultivariatePowerSeries/PuiseuxSeries/src/makeCompatible.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/multiplication.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/exponentiation.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/multiplicativeInverse.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/addition.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/subtraction.mm"
$include "MultivariatePowerSeries/PuiseuxSeries/src/equality.mm"

end module; 
