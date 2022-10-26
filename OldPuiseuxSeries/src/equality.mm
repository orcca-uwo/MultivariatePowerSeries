# ApproximatelyEqual
# to check approximately the equality of two OldPuiseuxSeriesObject
# deg is the maximum precision; base on that _self and other will be computed and compared 
# if deg = -2, then this algorithm will compare the first max(_self:-deg, other:-deg) homogeneous parts,
# elif deg = -1, returns always true, else, compare _self and other up to deg.
# WE NEED TO CHANGE THE DEGREE FOR SOMETHING ELSE IN THE FUTURE
export ApproximatelyEqual ::static := proc(_self :: OldPuiseuxSeriesObject, 
                            other :: OldPuiseuxSeriesObject,
                            deg :: integer := -2,
                            $)
    local A, B, new_expDivs; 

    # We modify the exp of _self and other to make them compatible
    new_expDivs, A, B := _self:-BinaryMakeCompatible(_self, other);

    return A:-ApproximatelyEqual(A, B, deg);
end proc;