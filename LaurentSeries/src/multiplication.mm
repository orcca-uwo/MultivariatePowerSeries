# BinaryMultiply
# to multiply two objects
export BinaryMultiply ::static := proc(_self :: LaurentSeriesObject,
                                other :: LaurentSeriesObject,
                                $)
local newOrdCV, newRays, pso1, pso2, newPso, newOrd, e, v;

    # We modify the exp of _self and other to make them compatible.
    # We get a new set of rays in which pso1 and pso2 can be operated.
    newOrd, newOrdCV, newRays, pso1, pso2 := LaurentSeriesObject:-BinaryMakeCompatible(_self, other);

    # We create a new PowerSeriesObject
    newPso := (pso1):-BinaryMultiply(pso1, pso2);

    # We multiply the monomials of _self and other
    e := [seq( v=_self:-e[v]+ other:-e[v], v in newOrd)];

    return Object(LaurentSeriesObject, newPso, newOrd, newOrdCV, newRays, e);
end proc;

# Multiply pso by a monomial given by ray and ord
export MultiplyByMonomial::static := proc(_self :: LaurentSeriesObject,
                                pso :: PowerSeriesObject,
                                ray :: list(integer),
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

    local mono := mul(ordCV^~result);
    mono := MultivariatePowerSeries:-PowerSeries(mono);

    return pso:-BinaryMultiply(pso, mono);

end proc;


