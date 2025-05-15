PowerSeriesWrapper := module()
option object;
export
    # pso: the power series to be wrapped.
    pso;

# ModuleCopy: is called when a new object of the same class is 
# created via a call to Object
local ModuleCopy ::static := proc(new :: PowerSeriesWrapper,
                                    old :: PowerSeriesWrapper,
                                    pso, $)
    if _npassed > 2 then
        new:-pso := pso;
    elif _npassed= 2 then
        new:-pso := old:-pso;
    else
        error "you cannot copy the original OdeObject object";
    end if;

end proc;

export `diff`::static := proc(_self::PowerSeriesWrapper, v::name, $)
    'Object'('PowerSeriesWrapper', _self:-pso:-Differentiate(_self:-pso, v));
end proc;

export Truncate::static := proc(_self::PowerSeriesWrapper, n::integer, $)
    return _self:-pso:-Truncate(_self:-pso, n);
end proc;

# Module Deconstruct
local ModuleDeconstruct :: static := proc(_self :: PowerSeriesWrapper, $)
    return 'Object'('PowerSeriesWrapper', _self:-pso);
end proc;
end module;
