# BinaryAdd
# to add two PuSO.
# We rely in the BinaryMakeCompatible to look for
# a cone that contains the cone associated to _self and
# the cone associated to other. After the two PuSO have 
# the same associate cone, we can add them.
export
    BinaryAdd ::static := proc(_self :: PuiseuxSeriesObject,
                                other :: PuiseuxSeriesObject,
                                $)
local newOrdCV, newRays, pso1, pso2, newPso, newOrd, e, v, u_ray;
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

    # We factor the grevlex smallest monomial between _self and other
    if newSelf:-GrevLexGreater(newSelf, [seq(newSelf:-e[v], v in newSelf:-ord)], 
                                    [seq(newOther:-e[v], v in newSelf:-ord)]) then
        return newOther:-BinaryAdd(newOther, newSelf);
    end if;
    e := newSelf:-e;
    u_ray := [seq( newOther:-e[v]- e[v], v in newOther:-ord)];
    
    # We modify the exp of newSelf and newOther to make them compatible.
    # We get a new set of rays in which pso1 and pso2 can be operated.
    newOrd, newOrdCV, newRays, pso1, pso2 := PuiseuxSeriesObject:-BinaryMakeCompatible(newSelf, newOther, [u_ray]);
    local M := map(convert, newRays, ':-Vector');
    M := Matrix(M);

    # We make a change of variable in the monomial given by u_ray and multiply this by pso2
    pso2 := newSelf:-MultiplyByMonomial(newSelf, pso2, u_ray, M, newOrdCV);
    
    # We create a new PowerSeriesObject 
    newPso := (pso1):-BinaryAdd(pso1, pso2);

    return Object(PuiseuxSeriesObject, newPso, newOrd, newOrdCV, newRays, e);
end proc;

export NaryAdd ::static := proc(terms :: seq({PuiseuxSeriesObject,
                                              PowerSeriesObject,
                                              UnivariatePolynomialOverPowerSeriesObject,
                                              algebraic}),
                                $)
    local lterms := [terms];

    if not type(lterms, ':-list'(PuiseuxSeriesObject)) then
        local pusos, rest;
        pusos, rest := selectremove(type, lterms, PuiseuxSeriesObject);

        # We check if there is any rational polynomial.
        if membertype(':-ratpoly', rest) then
            local rat_poly;

            # We remove the rational polynomials from our list.
            rat_poly, rest := selectremove(type, rest, ':-ratpoly');
            # We add our rational polynomials.
            local add_rat_poly := normal(add(rat_poly));
            # We convert add_rat_poly to PuSO. 
            local rat_poly_as_puso := Object(PuiseuxSeriesObject, add_rat_poly);
            # We add rat_poly_as_puso to our list.
            pusos := [op(pusos), rat_poly_as_puso];
        end if;

        local sum_rest := UnivariatePolynomialOverPowerSeriesObject:-NaryAdd(op(rest));
            
        local rest_as_puso := Object(PuiseuxSeriesObject, 
                                     sum_rest:-ConvertToPowerSeries(sum_rest));

        lterms := [op(pusos), rest_as_puso];
    end if;

    local N := numelems(lterms);

    if N = 1 then
        return lterms[1];
    elif N =2 then
        return lterms[1]:-BinaryAdd(lterms[1], lterms[2]);
    else 
        local sum_terms := lterms[1];

        for local i from 2 to N do
            sum_terms := sum_terms:-BinaryAdd(sum_terms, lterms[i]);
        end do;

        return sum_terms;
    end if;

end proc;

#Puiseuxseries:-`+` n-ary operator  
export `+` ::static := proc(terms :: seq({PuiseuxSeriesObject,
                                          PowerSeriesObject,
                                          UnivariatePolynomialOverPowerSeriesObject,
                                          algebraic}), $)
    return PuiseuxSeriesObject:-NaryAdd(terms);
end proc;

