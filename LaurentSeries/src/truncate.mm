# Truncate 
# to truncate a multivariate laurent series
export
Truncate :: static := proc(_self :: LaurentSeriesObject,
                           deg :: integer := -2,
                           {mode :: identical(:-powerseries, :-absolute) := ':-powerseries'},
                           $)
  local x;  

  if mode = ':-powerseries' then

    return _self:-GetMonomial(_self) * subs(_self:-ChangeOfVariables(_self),
  _self:-pso:-Truncate(_self:-pso, deg));
  else
    ASSERT(mode = ':-absolute');
    local bound := _self:-AbsoluteDegreeBound(_self, deg);
    local result := _self:-pso:-Truncate(_self:-pso, bound);
    result := subs(_self:-ChangeOfVariables(_self), result);
    
    # Create a list of the monomials.
    local terms := convert(result, ':-list', ':-`+`');
    
    # "Multiply" by the monomial given by e
    terms := map(':-`*`', terms, _self:-GetMonomial(_self));

    # Write each element in 'terms' as a list of the factors.
    local terms_as_products := map(convert, terms, ':-list', ':-`*`');
    
    # Write each factor in the sub-lists of 'terms_as_products' as a 
    #list of the base and exponent.
    local terms_as_power_products := map2(map, convert, terms_as_products, ':-list', ':-`^`');
    
    # Keep in 'terms_as_power_products' only the sub-sub-lists that correspond to variables 
    #(i.e. remove the coefficient terms).
    terms_as_power_products := map2( select, type, terms_as_power_products, [':-name',':-rational'] );

    # Determine the total absolute degrees of each monomial in 'terms'.
    local absolute_total_degrees := map( add, map[3]( map[2], abs@op, 2, terms_as_power_products ) );

    # Determine the indices of the terms in 'terms_as_power_products' that are not larger than 'deg'.
    local include := select(i -> absolute_total_degrees[i] <= deg, [seq(1 .. numelems(terms))]);

    return add(terms[include]);    
  end if;
end proc;

local
AbsoluteDegreeBound :: static := proc(_self :: LaurentSeriesObject,
                                      deg :: integer,
                                      $)
  local signs;
  local n_rays := numelems(_self:-rays);
  local n_ord := numelems(_self:-ord);
  if n_rays = 0 then
    # No rays, so it doesn't matter where we truncate.
    return 0;
  end if;
  local a, i, j;
  local v := [seq(_self:-e[_self:-ord[i]] + add(a[j] * _self:-rays[j][i], j = 1 .. n_rays), i = 1 .. n_ord)];
  v := [op(v)];
  local inequalities := {
    seq(add(signs[i] * v[i], i = 1 .. n_rays) <= deg, signs in Iterator:-CartesianProduct([-1, 1] $ n_rays)),
    seq(a[i] >= 0, i = 1 .. n_rays)};

  try
    local solution := Optimization:-LPSolve(add(a[i], i=1..n_rays), inequalities, ':-maximize', ':-assume' = ':-integer', ':-depthlimit' = 1000);
    return solution[1];
  catch "no feasible point found for LP subproblem":
    # Since there's no feasible point, there are no terms within the requested region.
    return 0;
  end try;
end proc;