# ApproximatelyZero
# to check approximately (up to a specific precision) if _self is zero
# if deg = -2, then this algorithm will compare ONLY generated hpolys up to _self:-deg with zero polynomial,
# elif deg = -1, returns always true, else, compare _self:-hpoly[i] and 0 for i = 0, ..., deg.
export 
    ApproximatelyZero ::static := proc(_self :: LaurentSeriesObject,
                            deg :: integer := -2,
                            {mode :: identical(:-powerseries, :-absolute) := ':-powerseries'},
                            $)
        # if the algexpr equals 0
        local hpoly := _self:-Truncate(_self, deg, _options['mode']);

        return evalb(normal(hpoly1) = 0);

    end proc;

# ApproximatelyEqual
# to check approximately the equality of two LaurentSeriesObject
# deg is the maximum precision; base on that _self and other will be computed and compared 
# if deg = -2, then this algorithm will compare the first max(_self:-deg, other:-deg) homogeneous parts,
# elif deg = -1, returns always true, else, compare _self and other up to deg.
export 
    ApproximatelyEqual ::static := proc(_self :: LaurentSeriesObject, 
                            other :: LaurentSeriesObject,
                            deg :: integer := -2,
                            {mode :: identical(:-powerseries, :-absolute) := ':-powerseries'},
                            $)
        # if the algexprs equal 
        local hpoly1 := _self:-Truncate(_self, deg, _options['mode']);
        local hpoly2 := other:-Truncate(other, deg, _options['mode']);

        return evalb(normal(hpoly1 - hpoly2) = 0);
    end proc;