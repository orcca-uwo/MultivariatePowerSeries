#Grevlex computations.

# Check if ray1 is grevlex greater than ray2.
export GrevLexGreater::static := proc(_self::PuiseuxSeriesObject, 
                                      ray1::{list(rational), Vector}, 
                                      ray2::{list(rational), Vector}, $)
    local s1, s2, i;

    if numelems(ray1) <> numelems(ray2) then
        error("the dimension of ray %1 must agree with the dimension of ray %2", ray1, ray2);
    end if;

    s1 := add(ray1);
    s2 := add(ray2);

    if s1 < s2 then
        return false;
    elif s1 > s2 then
        return true;
    else
        for i from numelems(ray1) to 1 by -1 do
            if ray1[i] < ray2[i] then
                return true;
            elif ray1[i] > ray2[i] then
                return false;
            end if;
        end do;
    end if;

    return false;
end proc; 

# Check if a list of exponents is grevlex positive.
export Positive :: static:= proc(_self :: PuiseuxSeriesObject, 
                                 ray::{list(rational), Vector}, $)
    local i;
    local zero := [seq(0, 1..numelems(ray))];

    return _self:-GrevLexGreater(_self, ray, zero);
end proc;

# Look for the grevlex smallest element in a list.
export Smallest::static := proc(_self :: PuiseuxSeriesObject, 
                                rays::list(list(rational)), $)
    local ray;
    local sm := rays[1];

    for ray in rays do
        if _self:-GrevLexGreater(_self, sm, ray) then
            sm := ray;
        end if;
    end do;

    return sm;
end proc;

# Look for the grevlex smallest grevlex monomial in a list of monomials.
local SmallestGrevlexMonomial::static := proc(_self :: PuiseuxSeriesObject, 
                                               mons::list(ratpoly), 
                                               v_order::list(name), 
                                               sw::boolean := true,
                                               $)

    local as_products := map(convert, mons, ':-list', ':-`*`');
    local as_power_lists := map2(map, convert, as_products, ':-list', ':-`^`');
    local as_power_eqns := map2(map, lst -> lst[1] = lst[2], as_power_lists);
    local exponents_except_0 := map(subs, as_power_eqns, v_order);
    local rays := subs(v_order =~ 0, exponents_except_0);
    
    if sw then                  
        return _self:-Smallest(_self, rays), rays;
    else
        return _self:-Smallest(_self, rays);
    end if;
end proc;

# Look for integer a integer vector that's grevlex-less than ray.
local LookForGreatestGrevlexLess::static := proc(_self::PuiseuxSeriesObject, 
                                                ray::list(rational), $)
    local i;

    # Weight of ray
    local s := add(ray);

    local newRay := Array(ray);

    for i from numelems(ray) to 2 by -1 do
        newRay[i] := ceil(ray[i]);
    end do;

    newRay[1] := s - add(newRay[2..]);

    return convert(newRay, ':-list');
end proc; 


local GrevlexFunctional::static := proc(i :: posint,
                                        v :: {list, Vector}, $)
  local j;
  return add(v[j], j = 1..(numelems(v) + 1 - i));
end proc;