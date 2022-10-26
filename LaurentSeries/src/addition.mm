# BinaryAdd
# to multiply two objects
export
    BinaryAdd ::static := proc(_self :: LaurentSeriesObject,
                                other :: LaurentSeriesObject,
                                $)
local newOrdCV, newRays, pso1, pso2, newPso, newOrd, e, v, u_ray;
    # We factor the grevlex greatest monomial between _self and other
    if _self:-ord <> other:-ord then
        error "the order %1 must be equal to the order %2", _self:-ord, other:-ord;
    end if;

    if _self:-GrevLexGreater(_self, [seq(_self:-e[v], v in _self:-ord)], 
                                    [seq(other:-e[v], v in _self:-ord)]) then
        return other:-BinaryAdd(other, _self);
    end if;

    e := _self:-e;
    u_ray := [seq( other:-e[v]- e[v], v in other:-ord)];
    
    # We modify the exp of _self and other to make them compatible.
    # We get a new set of rays in which pso1 and pso2 can be operated.
    newOrd, newOrdCV, newRays, pso1, pso2 := LaurentSeriesObject:-BinaryMakeCompatible(_self, other, [u_ray]);

    local M := map(convert, newRays, ':-Vector');
    M := Matrix(M);

    # We make a change of variable in the monomial given by u_ray and multiply this by pso2
    pso2 := _self:-MultiplyByMonomial(_self, pso2, u_ray, M, newOrdCV);;

    # We create a new PowerSeriesObject 
    newPso := (pso1):-BinaryAdd(pso1, pso2);


    return Object(LaurentSeriesObject, newPso, newOrd, newOrdCV, newRays, e);
end proc;