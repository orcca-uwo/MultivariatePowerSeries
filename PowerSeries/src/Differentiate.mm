# Differentiate generator.
local diff_gen ::static := proc(_self :: PowerSeriesObject,  
                            d :: nonnegint,
                            $)  
    local homogeneous_part := HOMOGENEOUS_PART(_self:-ancestors:-A, d+1);
    local differentiated_term := diff(homogeneous_part, _self:-ancestors:-v);  

    ASSERT(degree(differentiated_term,_self:-ancestors:-v) = d or homogeneous_part = 0);

    return differentiated_term;
end proc;


export Differentiate :: static := proc(_self :: PowerSeriesObject,v :: name, $)
    if numelems(_self:-Variables(_self)) > 1 then
        error "expected univariate Power Series, but recieved %1",_self;
    end if;

    if evalb(v in _self:-Variables(_self)) = false then
        return MultivariatePowerSeries:-PowerSeries(0);
    end if;

    local first_term := diff(HOMOGENEOUS_PART(_self,1),v);
    local hpoly := Array(0..0,[first_term]);
    local deg := 0;

    return Object(PowerSeriesObject, hpoly, deg, diff_gen, _self:-Variables(_self),
          ["A" = _self,"v"=v], ifelse(membertype(undefined, [_self:-algexpr]), undefined,
                               diff(_self:-algexpr,v)));
end proc;