export
TaylorShift :: static := proc(_self :: PowerSeriesObject,
                              substitutions :: {thistype, set, list}(name = COEFFICIENT_TYPE),
                              {method :: identical(diff, series) := series},
                              {force :: truefalse := false},
                              $)
local my_substitutions := (if type(substitutions, ':-list') then
                               substitutions;
                           elif type(substitutions, ':-set') then
                               [op(substitutions)];
                           else
                               [substitutions];
                           end if);

local variables_substituted := convert(map(lhs, my_substitutions), ':-set');
    if numelems(variables_substituted) <> numelems(my_substitutions) then
        error "these variables occurred multiple times on the left hand side: %1",
        select((v, lhss) -> numboccur(v, lhss) > 1, variables_substituted, map(lhs, my_substitutions));
    end if;
    
    my_substitutions := remove(type, my_substitutions, ':-anything' = 0);
    my_substitutions := select((s, sv) -> member(lhs(s), sv), my_substitutions, _self:-vars);

    if numelems(my_substitutions) = 0 then
        return _self;
    elif _self:-algexpr = undefined then
        error "cannot compute the Taylor shift of a power series of which the analytic expression is "
        "not known";
    end if;

    my_substitutions := select((s, ae) -> depends(ae, lhs(s)), my_substitutions, _self:-algexpr);
    if numelems(my_substitutions) = 0 then
        # In case there were only variables that are present but that the algexpr doesn't depend on.
        return _self;
    end if;
    
local my_substitutions_algexpr := map(s -> lhs(s) = lhs(s) + rhs(s), my_substitutions);
    try
        return FromAlgebraicExpression(subs(my_substitutions_algexpr, _self:-algexpr),
                                       _options['method', 'force']);
    catch "could not form power series: the expression %1 appears to have a pole at zero",
        "not invertible":
        error "invalid Taylor shift: tried to shift to a pole of the analytic expression";
    end try;
end proc;
