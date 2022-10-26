# BinaryMultiply
# to multiply two objects
export
    BinaryMultiply ::static := proc(_self :: {OldPuiseuxSeriesObject, PowerSeriesObject
                                        , UnivariatePolynomialOverPowerSeriesObject},
                                other :: {OldPuiseuxSeriesObject, PowerSeriesObject
                                        , UnivariatePolynomialOverPowerSeriesObject},
                                $)
local new_pso, new_expDivs, A, B;
    new_expDivs := table();

    # We modify the exp of _self and other to make them compatible
    new_expDivs, A, B := OldPuiseuxSeriesObject:-BinaryMakeCompatible(_self, other);

    # We create a new PowerSeriesObject
    new_pso := (A):-BinaryMultiply(A, B);


    return Object(OldPuiseuxSeriesObject, new_pso, new_expDivs);
end proc;