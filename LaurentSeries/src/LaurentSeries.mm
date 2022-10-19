local module LaurentSeriesObject()
option object;

local 
	pso, # Power Series object
	rays::list(list(integer)), #list containing the information
					  # of the exponent of each variable
    e::table, # Exponent of the monomial multiplying the lso
    ord::list, # list containing the order of the Laurent series variable
    ordCV::list, # list containing the order of the power series variable
    dstyle::PSDS_TYPE, #a display style for this object
    # Degree bounds for the inverse function
    degree_bound_pso_static :: static(nonnegint) := 10,
    degree_bound_pso :: {nonnegint, identical(undefined)} := undefined,
    degree_bound_extra_static :: static(nonnegint) := 10,
    degree_bound_extra :: {nonnegint, identical(undefined)} := undefined;


local ModuleApply::static := proc()
	return Object(LaurentSeriesObject, _passed);
end proc;	

# Constructor of LSO using rays
local CreateLSOFromRays::static := proc(new :: LaurentSeriesObject, 
                                 old :: LaurentSeriesObject,
                                 pso::{PowerSeriesObject, polynom},
                                 ord::list(name),
                                 ordCV::list(name) := [],
                                 rays::list(list(integer)) := [],
                                 e::{table, list(name=integer)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $) option overload(':-callseq_only');
    local l, i, j, eq;
    
    if type(pso, ':-polynom') then
        new:-pso := MultivariatePowerSeries:-PowerSeries(pso);
    else
        new:-pso := pso;
    end if;
    new:-dstyle := dstyle;
    new:-ord := ord;

    # We check if there is an order for the power series variable.
    # If there is no order, then we use the one in pso
    if ordCV = [] then
        if nops(new:-pso:-Variables(new:-pso))<=1 then
            new:-ordCV := [op(new:-pso:-Variables(new:-pso))]; 
        else 
            error "an order for the power series variable must be specified";
        end if;
    elif new:-pso:-Variables(new:-pso) = convert(ordCV, set) then
        new:-ordCV := ordCV;
    else
        error "the variables in %1 must be equal to %2", new:-pso, ordCV;
    end if;

    # We check if there are rays for the change of variables
    if rays = [] then
        if numelems(new:-ordCV)=numelems(new:-ord) then
            new:-rays := [seq( [seq(ifelse(j=i, 1, 0), j=1..nops(new:-ord))] ,i=1..nops(new:-ord))];
        elif numelems(new:-ordCV) = 0 then
            new:-rays := [];
        else
            error "the rays must be specified if the number of variables in the Laurent series variable and the power series variable is different";
        end if;
    elif numelems(new:-ordCV) <> numelems(rays)  then
        error "number of rays must agree with the number of variables in %1 or be empty", new:-pso;
    else
        # We check that the rays have the correct dimension for the 
        # change of variable
        for i in rays do 
            if numelems(i)<>numelems(ord) then
                error "all the rays must have dimension %1", numelems(ord);
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

# Constructor of LSO based in equations. No orders needed.
local CreateLSOFromEqu::static :=  proc(new :: LaurentSeriesObject, 
                                 old :: LaurentSeriesObject,
                                 pso::{PowerSeriesObject, polynom},
                                 mp::list('`=`'( name, {'`*`'( { name, name^integer } ), name, name^integer})) := [],
                                 e::{table, list(name=integer)} := [],
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
    local as_products := map(convert, map(rhs, mp), ':-list', `*`);
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

# Constructor of LSO based in equations and a rational polynomic function. 
# No orders needed.
local CreateLSOFromRatFunc::static :=  proc(new :: LaurentSeriesObject, 
                                 old :: LaurentSeriesObject,
                                 f::ratpoly,
                                 mp::list('`=`'( name, {'`*`'( { name, name^integer } ), name, name^integer})) := undefined,
                                 e::{table, list(name=integer)} := [],
                                 dstyle::PSDS_TYPE := [],
                                 $) option overload(':-callseq_only');
    
    local g, e1, h, g_reduced, q, const, V, v, rays_q_CV;

    new:-dstyle := dstyle;

    # We check if there is a list of equations as input.
    if mp = undefined then
        local new_mp := [seq(v=v, v in (indets(f, ':-name') minus {':-constants'}) )];
        return LaurentSeriesObject(f, new_mp, e);

    # We check that there is an equation for each variable in f
    elif not convert(map(lhs, mp), ':-set') = (indets( f, ':-name') minus {':-constants'}) then
        error "the variables in %1 must be equal to the variables in %2", mp, f;
    end if;

    # We convert the rhs of the equations as a list separated by ^
    local as_powers :=  new:-ComputeAsPowers(mp);
    # We generate the order for the cone
    local ord := map2(op, 1, map(op, as_powers));
    new:-ord := convert(convert(ord, ':-set'), ':-list');

    # We generate the rays based in the order and as_powers.
    # Note that this rays could change for two possible reasons:
    # We simplify one of the variables in f.
    # We do not have a intervible polynomial in the denominator.
    local rays := new:-GenerateRays(new, as_powers);

    # We check that the rays are grevlex positive with the given order
    if not new:-CheckRays(new, rays) then
        error "all the rays in %1 must be grevlex(%2) positive", rays, new:-ord;
    end if;

    local local_mp := mp;
    local gcd_h := 1;
    local sm_ray := [];

    # We write f as g_reduced/(gcd_h*q), i.e.,
    # we simplify the fraction and factorize the denominator.
    # Note that we need to use gcd_h in the computation of e.
    g_reduced, q, gcd_h := LaurentSeriesObject:-SimplifyFraction(LaurentSeriesObject, 
                                                                f);

    # Variables of q. Note that this variables could be different 
    # from the variables of f. 
    V := convert( indets( q, ':-name' ) minus {':-constants'}, ':-list' );
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
        sm_ray, rays_q_CV := LaurentSeriesObject:-SmallestGrevlexMonomial(LaurentSeriesObject, 
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

    # We create a new PSO.
    new:-pso := MultivariatePowerSeries:-PowerSeries(g_reduced/q);

    # We define the exponents for the monomial that multiplies lso
    new:-e := new:-MakeMonomialTable(new, e);

    # If gcd_h is different than 0, we need to devide e by it.
    if gcd_h <> 1 then
        # We apply the change of variable.
        local gcd_h_cv := subs(mp, gcd_h);
        # We see it as a list of the for [name, integer].
        local gcd_h_as_product := convert(gcd_h_cv, ':-list', ':-`*`');
        local gcd_h_as_powers := map(convert, gcd_h_as_product, ':-list', ':-`^`');

        # We substract the previous information from the current e.
        for local p in gcd_h_as_powers do
            new:-e[p[1]] := new:-e[p[1]] - p[2]; 
        end do;
    end if;

    # we need to simplify sm_rays from e.
    if sm_ray <> [] then
        for local i from 1 to numelems(new:-ord) do
            new:-e[new:-ord[i]] := new:-e[new:-ord[i]] - sm_ray[i]; 
        end do;
    end if;
end proc;

local CreateLSO_NoSuitableMethod::static := proc() option overload(':-callseq_only');
    error "invalid input: arguments to %1, %2, do not match any of the accepted calling sequences", procname, [ _passed ] ;
end proc;

# Constructor of the class
local ModuleCopy::static := overload([CreateLSOFromRays, CreateLSOFromEqu, 
                                      CreateLSOFromRatFunc, CreateLSO_NoSuitableMethod]);

# We check if some variable was eliminated during simplification.
# If so, we compute a new mp and rays without the information of
# the deleted variable.
local CheckVariableSimplification::static := proc(_self :: LaurentSeriesObject,
                                                  f::ratpoly,
                                                  new_f::ratpoly,
                                                  rays::list(list(integer)),
                                                  mp::list('`=`'( name, {'`*`'( { name, name^integer } ), name, name^integer})) := [],
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

local ComputeAsPowers::static := proc(mp::list('`=`'( name, {'`*`'( { name, name^integer } ), name, name^integer})) := [],
                                  $)
    # We convert the rhs of the equations as a list separated by *
    local as_products := map(convert, map(rhs, mp), ':-list', ':-`*`');
    # We separate the previous list as a list separated by ^
    local as_powers := map2(map, convert, as_products, ':-list', ':-`^`');

    return as_powers;
end proc;

# We compute the change of variables given by M 
local PolynomialChangeOfVariables::static := proc(_self :: LaurentSeriesObject,
                                        M::Matrix,  
                                        new_ordCV::list(name), 
                                        old_ordCV::list(name), 
                                        rays::list(list(integer)),
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
local SimplifyFraction::static := proc(_self::LaurentSeriesObject,
                                       f::ratpoly,
                                       $)

    local h_reduced, g_reduced, q;
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
    divide(h_expanded, gcd_h, 'q');

    return g_reduced, q, gcd_h;
end proc;


# Monomial generator
local MakeMonomialTable::static := proc(_self::LaurentSeriesObject, e)
    local l;

    if type(e, ':-table') then
        l := [entries(e, ':-nolist', ':-pairs')];
        if not type(l, ':-list'(':-name'=':-integer')) then
            error "invalid input: %1 expects its %-2 argument, %3, to consist of entries of type %4, but received %5", 
                                    procname, 4, ':-e', ':-name' =':-integer', eval(e, 1);
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

# Ray generator
local GenerateRay::static := proc(_self::LaurentSeriesObject, powlist, $)
    local n := numelems(_self:-ord);
    local result := Array(1 .. n);

    for local pair in powlist do
        local i := ListTools:-Search(pair[1], _self:-ord);
        result[i] := pair[2];
    end do;

    return convert(result, ':-list');
end proc;

# Rays generator
local GenerateRays::static := proc(_self::LaurentSeriesObject, powlist, $)
    local ls;
    return [seq(_self:-GenerateRay(_self, ls    ), ls in powlist)];
end proc;

# Funtion that checks if the rays are grevlex positive with the given order or not
local CheckRays::static := proc(_self::LaurentSeriesObject, rays::list(list(integer)), $)
    local r;

    for r in rays do 
        if not _self:-Positive(_self, r) then
            return false;
        end if;
    end do;

    return true;
end proc; 

# Display the change of variable
export ChangeOfVariables::static := proc(_self :: LaurentSeriesObject, $)
    local r, i, j;

    r := [seq(_self:-ordCV[j]
            =mul(_self:-ord[i]^_self:-rays[j][i], i = 1 .. nops(_self:-ord)), j=1..nops(_self:-ordCV))];

    return r;
end proc;

# SetDisplayStyle
# to set the dstyle of the input object
export SetDisplayStyle :: static:= proc(_self :: LaurentSeriesObject, s :: PSDS_TYPE, $)
    _self:-dstyle := s; 
end proc;

export Display :: static := proc(_self :: LaurentSeriesObject,
                                 user_dstyle :: list := [],
                                 output :: identical("typeset", "string") := "typeset",
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
        _self:-pso, mystyle, ifelse(output = "typeset", "full", "fullstring"),
        ':-objectname' = "LaurentSeries", ':-substitutions' = _self:-ChangeOfVariables(_self),
        ':-cofactor' = mul(x^_self:-e[x], x in _self:-ord));
end proc;


# ModulePrint: print _self :: LaurentSeriesObject 
local ModulePrint :: static := proc(_self :: LaurentSeriesObject, $)
    if IsWorksheetInterface() then
        return _self:-Display(_self);
    else
        return convert(_self:-Display(_self, [], "string"), ':-name');
    end if;
end proc;

# Module Deconstruct
local ModuleDeconstruct :: static := proc(_self :: LaurentSeriesObject, $)
    return 'Object'('LaurentSeriesObject', _self:-pso, _self:-ord, _self:-ordCV, _self:-rays,
                    [indices(_self:-e, ':-pairs')]);
end proc;	 

# Basic routines
$include "MultivariatePowerSeries/LaurentSeries/src/gets.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/grevlex.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/truncate.mm"
# Arithmetic 
$include "MultivariatePowerSeries/LaurentSeries/src/makeCompatible.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/multiplication.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/multiplicativeInverse.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/addition.mm"
$include "MultivariatePowerSeries/LaurentSeries/src/equality.mm"


end module; 
