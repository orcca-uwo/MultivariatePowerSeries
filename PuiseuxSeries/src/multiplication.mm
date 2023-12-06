# BinaryMultiply
# to multiply two objects.
# We rely in the BinaryMakeCompatible to look for
# a cone that contains the cone associated to _self and
# the cone associated to other. After the two PuSO have 
# the same associate cone, we can multiply them.
export BinaryMultiply ::static := proc(_self :: PuiseuxSeriesObject,
                                other :: PuiseuxSeriesObject,
                                $)
local newOrdCV, newRays, pso1, pso2, newPso, newOrd, e, v;
    local newSelf := _self;
    local newOther := other;

    if _self:-ord <> other:-ord then
        # If the ord of _self and other is compatible, we extend 
        # other and _self to the Puiseux Series space generated 
        # by this new order.
        local Ord_ext := _self:-MakeOrdCompatible(_self, _self:-ord, other:-ord);
        
        newOther := other:-ExtendPuiseuxSeriesObject(other, Ord_ext);
        newSelf := _self:-ExtendPuiseuxSeriesObject(_self, Ord_ext);
    end if;

    # We modify the exp of newSelf and newOther to make them compatible.
    # We get a new set of rays in which pso1 and pso2 can be operated.
    newOrd, newOrdCV, newRays, pso1, pso2 := PuiseuxSeriesObject:-BinaryMakeCompatible(newSelf, newOther);

    # We create a new PowerSeriesObject
    newPso := (pso1):-BinaryMultiply(pso1, pso2);

    # We multiply the monomials of newSelf and newOther
    e := [seq( v=newSelf:-e[v]+ newOther:-e[v], v in newOrd)];

    return Object(PuiseuxSeriesObject, newPso, newOrd, newOrdCV, newRays, e);
end proc;

# Multiply pso by a monomial given by ray and ord
export MultiplyByMonomial::static := proc(_self :: PuiseuxSeriesObject,
                                pso :: PowerSeriesObject,
                                ray :: list(rational),
                                M :: Matrix,
                                ordCV :: list(name),
                                $)

    if type(ray, ':-list'(0)) then
        return pso;
    end if;                       

    local rayAsVector := convert(ray, ':-Vector');

    local result := LinearAlgebra:-LinearSolve(M, rayAsVector);

    ASSERT(andmap(type, result, ':-nonnegint'), 
                    "the result should be a monomial, but some exponents aren't nonnegative integers");

    local mono := mul(:-`~`[:-`^`](ordCV, ':-` $`', result));
    mono := MultivariatePowerSeries:-PowerSeries(mono);

    return pso:-BinaryMultiply(pso, mono);

end proc;

export NaryMultiply ::static := proc(factors :: seq({PuiseuxSeriesObject,
                                                     PowerSeriesObject,
                                                     UnivariatePolynomialOverPowerSeriesObject,
                                                     algebraic}),
                                $)
    
    local lfactors := [factors];;

    if not type(lfactors, ':-list'(PuiseuxSeriesObject)) then
        local pusos, rest;
        pusos, rest := selectremove(type, lfactors, PuiseuxSeriesObject);

        # We check if there is any rational polynomial.
        if membertype(':-ratpoly', rest) then
            local rat_poly;

            # We remove the rational polynomials from our list.
            rat_poly, rest := selectremove(type, rest, ':-ratpoly');
            # We multiply our rational polynomials.
            local mul_rat_poly := normal(mul(rat_poly));
            # We convert mul_rat_poly to PuSO. 
            local rat_poly_as_puso := Object(PuiseuxSeriesObject, mul_rat_poly);
            # We add rat_poly_as_puso to our list.
            pusos := [op(pusos), rat_poly_as_puso];
        end if;

        local mul_rest := UnivariatePolynomialOverPowerSeriesObject:-NaryMultiply(op(rest));
            
        local rest_as_puso := Object(PuiseuxSeriesObject, 
                                     mul_rest:-ConvertToPowerSeries(mul_rest));

        lfactors := [op(pusos), rest_as_puso];
    end if;

    local N := numelems(lfactors);

    if N = 1 then
        return lfactors[1];
    elif N =2 then
        return lfactors[1]:-BinaryMultiply(lfactors[1], lfactors[2]);
    else 
        local product_factors := lfactors[1];

        for local i from 2 to N do
            product_factors := product_factors:-BinaryMultiply(product_factors, 
                                                               lfactors[i]);
        end do;

        return product_factors;
    end if;

end proc;

# Puiseuxseries:-`*` n-ary operator 
export `*` ::static := proc(factors :: seq({PuiseuxSeriesObject, 
                                            PowerSeriesObject,
                                            UnivariatePolynomialOverPowerSeriesObject,
                                            algebraic}), $)
    return PuiseuxSeriesObject:-NaryMultiply(factors);
end proc;

