#Multiplicative Inverse
export Inverse ::static := proc(_self :: PuiseuxSeriesObject,
                                nonzero_pso_bound :: nonnegint := FAIL,
                                smallest_term_bound :: nonnegint := FAIL,
                                {useanalytic :: boolean := true},
                                $)
    local v, newPso, non_zero_term, bnd_pso1, bnd_pso2, 
          exp_non_zero_term, M, min_monomial_exps, mp,
          new_rays, new_ordCV, change_of_variable,
          sm_mon, i, mon_as_exp, r, s;

    local my_nonzero_pso_bound, my_smallest_term_bound;

    # We check bounds.
    if nonzero_pso_bound <> FAIL then
        my_nonzero_pso_bound := nonzero_pso_bound;
    elif _self:-nonzero_pso_bound <> undefined then
        my_nonzero_pso_bound := _self:-nonzero_pso_bound;
    else
        my_nonzero_pso_bound := _self:-nonzero_pso_bound_static;
    end if;

    if smallest_term_bound <> FAIL then
        my_smallest_term_bound := smallest_term_bound;
    elif _self:-smallest_term_bound <> undefined then
        my_smallest_term_bound := _self:-smallest_term_bound;
    else
        my_smallest_term_bound := _self:-smallest_term_bound_static;
    end if;

    # We invert the monomial multiplying our lso.
    local new_e := map(:-`*`, _self:-e, -1);

    # We get the analytic expression of the internal pso.
    local ana := _self:-pso:-GetAnalyticExpression(_self:-pso);

    # We check if _self:-pso is a unit.
    if _self:-pso:-IsUnit(_self:-pso) then
        newPso := MultivariatePowerSeries:-Inverse(_self:-pso);

        return Object(PuiseuxSeriesObject, newPso, _self:-ord, 
                                         _self:-ordCV, _self:-rays, new_e);
    # We check if ana is a rational function.
    elif ana <> undefined and type(ana, ':-ratpoly') and useanalytic then
        # We must create a new PuSO with ana^(-1) as an input.
        mp := _self:-ChangeOfVariables(_self);

        return Object(PuiseuxSeriesObject, ana^(-1), mp, new_e);
    else  
        #########
        # Puiseux Series adjustment.
        local rays := _self:-rays;
        local denominator := lcm(seq(seq(denom(s), s in r) , r in rays));

        # This could be an improvement of multiplying by denominator:
        #rays := [[1/2, 0], [1, -1/3], [1, 0], [1, -1/2]];
        #m := Matrix(rays);
        #k, n := upperbound(m);
        #lcm_so_far := 1;
        #lcms := Array(1..n);
        #for j from n to 1 by -1 do
        #    lcms[j] := lcm(lcm_so_far, seq(denom(m[i, j]), i=1..k));
        #    lcm_so_far := lcms[j];
        #end do;
        #denom_mat := LinearAlgebra:-DiagonalMatrix(lcms);

        # Copy of _self.
        local C := Object(PuiseuxSeriesObject, _self:-pso, 
                                             _self:-ord, 
                                             _self:-ordCV, 
                                            _self:-rays*denominator,
                                             copy(_self:-e));

        #########
        # We look for the first non-zero hom part in pso and then
        # we choose the smallest grevlex term in the LS order between then.
        non_zero_term, bnd_pso1 := C:-LookForNonZeroTermInPso(C, 
                                                              my_nonzero_pso_bound);
        
        # We compute the bound 2.
        # Since HasZeroSumRays(C) = false, we know that if there
        # is a grevlex smaller term than non_zero_term, it must have 
        # degree less than bnd_pso2. 
        if C:-HasZeroSumRays(C) = false then
            bnd_pso2 := infinity;
        else
            bnd_pso2 := bnd_pso1 +  my_smallest_term_bound;
        end if;

        # We look for the grevlex smallest term in lso and also
        # we save all the monomials that we check.
        min_monomial_exps, mon_as_exp := C:-LookForGrevlexSmallestTerm(C, bnd_pso1, 
                                                           bnd_pso2, non_zero_term,
                                                           C:-HasZeroSumRays(C));


        # We compute a set of rays for a cone in which lso is invertible.
        new_rays := C:-RaysFromTheInverse(C, min_monomial_exps, mon_as_exp);

        M := Matrix(new_rays);
        # We compute new variables for our pso.
        new_ordCV := C:-ComputeNewOrdCV(C, numelems(C:-rays), C:-ordCV);

        # We compute the old variables in terms of the new ones.
        change_of_variable := C:-PolynomialChangeOfVariables(C, M, new_ordCV, 
                                                               C:-ordCV, C:-rays);


        # We compute the smallest monomial of pso in terms of the new variables.
        sm_mon := C:-PolynomialChangeOfVariables(C, M, new_ordCV, 
                                                        [sm], [min_monomial_exps]);

        # We make a change of variable in pso.
        newPso := C:-pso:-Substitute(change_of_variable, C:-pso);
        # We factor out sm_mon from newPso.
        newPso := newPso:-FactorOutMonomial(newPso, rhs(sm_mon[1]));
        # We invert newPso.
        newPso := newPso:-Inverse(newPso);

        #########
        # Puiseux Series adjustment.
        if denominator<>1 then
            new_rays := (1/denominator)*new_rays;
            min_monomial_exps := (1/denominator)*min_monomial_exps;
        end if;
        #########

        # We 'divide' new_e by min_monomial_exp
        for i from 1 to numelems(C:-ord) do
            new_e[C:-ord[i]] := new_e[C:-ord[i]]- min_monomial_exps[i] ;
        end do;

        return Object(PuiseuxSeriesObject, newPso, C:-ord, new_ordCV, new_rays, new_e);
    end if;
end proc;

# We compute the rays of the cone from inverse of lso
# after having the minimum element of its support.
local RaysFromTheInverse::static := proc(_self::PuiseuxSeriesObject,
                                         sm_exp_list::list(integer),
                                         mon_as_exp::list(list(integer)),
                                         $)
    local r, c_prima, c_as_vector;
    local mons_mins_min := [seq(r-sm_exp_list, r in mon_as_exp)];
    
    local new_rays_1 := [op(_self:-rays), op(mons_mins_min)];
    local new_rays_2 := Array(1..0);

    local k := numelems(_self:-ordCV);
    local d := add(sm_exp_list) +1;

    local M := Matrix(_self:-rays)^+;
    
    for local c in Iterator:-BoundedComposition([d $ k], d) do
        c_as_vector := ArrayTools:-Alias(c, [k], ':-column');
        c_prima := M.c_as_vector;
        c_prima := convert(c_prima, ':-list');

        if _self:-GrevLexGreater(_self, c_prima, sm_exp_list) then
            new_rays_2 ,= c_prima - sm_exp_list;
        end if;
    end do;

    local new_rays := [seq(new_rays_1), seq(new_rays_2)];
    new_rays := _self:-MakeRaysCompatible(_self, new_rays,
                                          _self:-ord);

    return new_rays;
end proc;


# Get the exponent of a pso monomial as a multi-index
local GetDegreeAsMultiindex::static := proc(_self::PuiseuxSeriesObject, 
                                            mon::polynom,
                                            ordCV::list(name),
                                             $)
    # We get the exponent of mon as a vector.
    local exp_as_vector := map2(map2, degree, mon, ordCV);

    return exp_as_vector;
end proc;

# We get the total degree of a monomial in pso after applying the change of variable 
# given by rays.
local GetDegreeLSMonomial::static := proc(_self::PuiseuxSeriesObject, 
                                          mon::polynom, 
                                          ordCV::list(name)
                                          ,$)
    local exp_mon_as_vector;

    exp_mon_as_vector := Vector(_self:-GetDegreeAsMultiindex(_self, mon, ordCV));

    local M := Matrix(_self:-rays)^+;

    return add(M.exp_mon_as_vector);
end proc;

# To set a nonzero_pso_bound or a nonzero_pso_bound_static.
export SetNonzeroPowerSeriesDegreeBound ::static := proc(_self :: PuiseuxSeriesObject,
                                                  bound::nonnegint,
                                                  {instance :: truefalse := true},
                                                  $)
    local old_value;
    
    if instance then
        old_value := _self:-nonzero_pso_bound;
        _self:-nonzero_pso_bound := bound;
    else
        old_value := _self:-nonzero_pso_bound_static;
        _self:-nonzero_pso_bound_static := bound;
    end if;

    return old_value;
end proc;

# To set a smallest_term_bound or a smallest_term_bound_static.
export SetSmallestTermDegreeBound ::static := proc(_self :: PuiseuxSeriesObject,
                                                   bound::nonnegint,
                                                   {instance :: truefalse := true},
                                                   $)
    local old_value;
    
    if instance then
        old_value := _self:-smallest_term_bound;
        _self:-smallest_term_bound := bound;
    else
        old_value := _self:-smallest_term_bound_static;
        _self:-smallest_term_bound_static := bound;
    end if;

    return old_value;
end proc;

# We check if there is a ray in _self:-rays with weight zero
local HasZeroSumRays::static := proc(_self::PuiseuxSeriesObject, $)
    local r; 

    return orseq(add(r)=0, r in _self:-rays);
end proc;

# Look for first nonzero hom part, and then we select the smallest grevlex 
# term between then.
# The degree_bound represents how much we truncate _self:-pso.
local LookForNonZeroTermInPso::static := proc(_self :: PuiseuxSeriesObject,
                                              bound:: nonnegint,
                                              $)
    local hom_part, i, mons;
                                    
    for i from 0 to bound do
        hom_part := MultivariatePowerSeries:-HomogeneousPart(_self:-pso, i);

        if hom_part <> 0 then
            mons := convert(hom_part, ':-list', ':-`+`');
            mons := map2(subs, _self:-ChangeOfVariables(_self), mons);
            
            return _self:-SmallestGrevlexMonomial(_self, mons, _self:-ord, false), i;
        end if;
    end do;

    error "when inverting a Puiseux series, found no nonzero term of power series degree less than %1. The Puiseux series may be 0, or you may need to increase the bound using the `bound` argument to the Inverse command, or using the SetNonzeroPowerSeriesDegreeBound command", bound;
end proc;

# We look for the smallest grevlex term in _self using d1 and d2
# as bounds of the homogeneus part of _self:-pso
local LookForGrevlexSmallestTerm::static := proc(_self::PuiseuxSeriesObject,
                                                 d1::nonnegint,
                                                 d2::{nonnegint, identical(infinity)},
                                                 list_exp::list(rational),
                                                 haszerosumrays :: truefalse,
                                                 $)
    local hom_part, i, hom_part_as_list, min_exponent_list, sm_from_hom_part,
          hom_part_CV, temp_mon_as_exp;

    local monomials_as_exp := Array(1..0);

    local minimal_ray_degree := map(add, _self:-rays);
    minimal_ray_degree := min(minimal_ray_degree);

    local min_list_exp_degree := add(list_exp);

    # We look for the grevlex smallest exponent vector in lso.
    # We used d1 and d2 as the bounds for the hom part of the 
    # internal pso.
    min_exponent_list := list_exp;


    for i from d1+1 to d2 while haszerosumrays or 
                          min_list_exp_degree >= i * minimal_ray_degree do
        hom_part := MultivariatePowerSeries:-HomogeneousPart(_self:-pso, i);

        if hom_part <> 0 then
            hom_part_CV := subs(_self:-ChangeOfVariables(_self), hom_part);
            hom_part_as_list := convert(hom_part_CV, ':-list', ':-`+`');

            sm_from_hom_part, temp_mon_as_exp:=_self:-SmallestGrevlexMonomial(_self, 
                                                            hom_part_as_list,
                                                            _self:-ord);

            monomials_as_exp ,= op(temp_mon_as_exp);

            if _self:-GrevLexGreater(_self, min_exponent_list, sm_from_hom_part) then
                min_exponent_list := sm_from_hom_part;
                min_list_exp_degree := add(min_exponent_list);
            end if;
          
        end if;
    end do;  
    

    return min_exponent_list, convert(monomials_as_exp, ':-list');                                           
end proc;
