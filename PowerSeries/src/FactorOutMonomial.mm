export
FactorOutMonomial :: static := proc(_self :: PowerSeriesObject,
                                    m :: {{thistype, :-`*`}({name, name^posint}), 1},
                                    $)
    if m = 1 then
        return _self;
    end if;

local deg := degree(m);
local nondivisible := remove(`=`, convert(
    UP_TO_DEGREE_ARRAY(_self, min(deg - 1, _self:-deg)), ':-list'), 0);
    if numelems(nondivisible) > 0 then
        error "found terms not divisible by %1 when trying to factor it out: %2",
        m, nondivisible[1];
    end if;

    _self:-ensure_degree(_self, deg);
local hpoly := map(expr -> expand(expr / m), ArrayTools:-Alias(_self:-hpoly, deg,
                                                               [0 .. _self:-deg - deg]));

local algexpr, vars;
    if _self:-algexpr = undefined then
        algexpr := _self:-algexpr;
        vars := _self:-vars;
    else
        algexpr := _self:-algexpr / m;
        vars := select[2](has, algexpr, _self:-vars);
    end if;
                      
    return Object(PowerSeriesObject, hpoly, _self:-deg - deg, factor_out_monomial_gen, vars,
                  ["parent" = _self, "monomial" = m, "degree" = deg], algexpr);
end proc;

local
factor_out_monomial_gen :: static := proc(_self :: PowerSeriesObject,
                                          d :: nonnegint,
                                          $)
local ancestors := _self:-ancestors;

local hom_part := HOMOGENEOUS_PART(ancestors:-parent, d + ancestors:-degree);
local result := expand(hom_part / ancestors:-monomial);
    if not type(result, ':-polynom'(COEFFICIENT_TYPE, _self:-vars)) then
        error "in attempt to factor out monomial %1, found terms %2, which are not divisible by it",
        ancestors:-monomial, hom_part;
    end if;

    return result;
end proc;
