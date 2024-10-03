## ODE generators.

# Linear ode with constant coefficients generator.
local cons_coeff_ode_gen ::static := proc(_self :: PowerSeriesObject, 
                                        d :: nonnegint,
                                        $)
    # We get the OdeObject "ancestor" of _self.
    local edo := _self:-ancestors:-ode;
    local ind_var := edo:-GetIndependentVariable(edo);
    local init_cond_sym := edo:-GetInitialConditions(edo);
    local init_cond := edo:-ConvertConstantSymbolsToD(edo, init_cond_sym);
    local ord := edo:-GetOrder(edo);
    local coefs_of_ODE := edo:-GetCoeffs(edo);
    local i;   

    if d < ord then  
        # Return the Taylor coefficient up to our initial conditions.
        return init_cond[d]/d!*ind_var^d;
    else
        # Use our recursion formula to compute terms past the initial conditions.
        # The terms computed up to now.
        local hpoly := _self:-hpoly;
        local e := d-ord;
        local sequence := [seq(coeff(hpoly[e+i], ind_var, e+i)*coefs_of_ODE[i]*(e+i)!/(e)!, 
                                i=0 .. ord-1)]; 

        local s := add(coeff(hpoly[e+i], ind_var, e+i)*coefs_of_ODE[i]*(e+i)!/(e)!, 
                                i=0 .. ord-1);

        return -s*(e)!/(coefs_of_ODE[ord]*d!)*ind_var^d;
    end if; 
end proc;

# Linear ODE with polynoms coefficients generator
local poly_coeff_ode_gen ::static := proc(_self :: PowerSeriesObject, 
                                        d :: nonnegint,
                                        $)
    # We get the OdeObject "ancestor" of _self.
    # As well as all of its attributes
    local edo := _self:-ancestors:-ode;
    local ind_var := edo:-GetIndependentVariable(edo);
    local dep_var := edo:-GetDependentVariable(edo);
    local init_cond_sym := edo:-GetInitialConditions(edo);
    local init_cond := edo:-ConvertConstantSymbolsToD(edo, init_cond_sym);
    local coeffs := edo:-GetCoeffs(edo);
    local ord := edo:-GetOrder(edo);
    

    local i,k;
    local constant_coeff := coeffs[ord][0];
    local max_degree := (ArrayTools:-Size(coeffs[0]))[2]-1;
    local hpoly := _self:-hpoly;

    # Then we assign all initial conditions to the first a[ord-1] terms.
    if d < ord then
        return init_cond[d]/d!*ind_var^d;
    # There is a special condition when min = 0, in this case we ignore the
    # double summation and use single sum we have, similar to case of
    # constant coefficients
    elif min(max_degree,d-ord)= 0 then
        local sequence := add(coeff(hpoly[d-ord+i],ind_var,d-ord+i)*coeffs[i][0]*(d-ord+i)!/(d-ord)! 
                                ,i=0 .. ord-1);
        return -(d-ord)!/(d!*constant_coeff)*sequence*ind_var^(d);
    # Then finally we use the complete algorithm to compute the other terms.
    else
        # This is the single summation.
        local s1 := add(coeff(hpoly[d-ord+i],ind_var,d-ord+i)*coeffs[i][0]*(d-ord+i)!/(d-ord)!, 
                                i=0 .. ord-1);
        # This is the double summation.
        local s2 := add(add(coeff(hpoly[d-ord-k+i],ind_var,d-ord-k+i)*coeffs[i][k]*(d-ord-k+i)!/(d-ord-k)!,i=0..ord)
                            ,k=1.. min(max_degree,d-ord));
        # Final formula all together
        local total := -(d-ord)!/(d!*constant_coeff)*(s1+s2)*ind_var^(d);
        
        return total;
    end if; 
end proc;

## Power series from ode constructors.
# Function for ODE solution with constant coefficients.
local from_linear_coefficient_ode::static := proc(ind_var::name, 
                                                    dep_var::name, ode::set(equation), 
                                                    {exp_pnt :: {integer, name} := 0}, $)
    # NOTE: if exp_pnt is different than 0, OdeManipulator applies 
    # the change pf variable [ind_var=ind_var-exp_pnt] to ode.
    local ode_obj := OdeManipulator:-OdeObject(ind_var, dep_var, ode, _options['exp_pnt']);

    local deg := 0;
    local init_cond := ode_obj:-GetInitialConditions(ode_obj);
    local var := ode_obj:-GetIndependentVariable(ode_obj);

    local a0_sym := init_cond[0];
    local a0 := ode_obj:-ConvertConstantSymbolsToD(ode_obj, a0_sym);
    local hpoly := Array(0..0, [a0]);
    
    local ode_ps := Object(PowerSeriesObject, hpoly, deg, cons_coeff_ode_gen, {var},
          ["ode" = ode_obj], undefined);

    return ode_ps;
end proc:

# Function for ODE solution with polynomial coefficients.
local from_polynom_coefficient_ode :: static := proc(ind_var::name, 
                                                      dep_var::name, ode::set(equation), 
                                                      {exp_pnt :: {integer, name} := 0}, $)
  local ode_obj := OdeManipulator:-OdeObject(ind_var,dep_var,ode, _options['exp_pnt']);
  
  local deg := 0;
  local init_cond := ode_obj:-GetInitialConditions(ode_obj);
  local var := ode_obj:-GetIndependentVariable(ode_obj);
  local ord := ode_obj:-GetOrder(ode_obj);
  local constant_coeff := ode_obj:-GetCoeffs(ode_obj)[ord][0];

  local a0_sym := init_cond[0];
  local a0 := ode_obj:-ConvertConstantSymbolsToD(ode_obj, a0_sym);
  local hpoly := Array(0..0, [a0]);

  # For our algorithm, the leading polynomial must be invertible
  # So we throw an error if that's not the case
  if constant_coeff = 0 then
    error "leading polynomial not invertible";
  end if;
    
    return Object(PowerSeriesObject, hpoly, deg, poly_coeff_ode_gen, {var},
          ["ode" = ode_obj], undefined);

end proc:
