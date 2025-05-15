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
    local oderhs := edo:-GetRHS(edo);
    local i;   
    local hpoly := _self:-hpoly;
    local e := d-ord;

    local oderhs_coeff := 0;
    if ord <= d and d <= degree(oderhs, ind_var) + ord then
      oderhs_coeff := coeff(oderhs, ind_var, d-ord);
    end if;  

    if d < ord then  
        # Return the Taylor coefficient up to our initial conditions.
        return init_cond[d]/d!*ind_var^d;
    else
        # Use our recursion formula to compute terms past the initial conditions.
        # The terms computed up to now. 

        local s := add(coeff(hpoly[e+i], ind_var, e+i)*coefs_of_ODE[i]*(e+i)!/(e)!, 
                                i=0 .. ord-1);

        return (oderhs_coeff-s)*e!/(coefs_of_ODE[ord]*d!)*ind_var^d;
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
    local oderhs := edo:-GetRHS(edo);
    

    local i,k;
    local constant_coeff := coeffs[ord][0];
    local max_degree := (ArrayTools:-Size(coeffs[0]))[2]-1;
    local hpoly := _self:-hpoly;

    local oderhs_coeff := 0;
    if ord <= d and d <= degree(oderhs, ind_var) + ord then
      oderhs_coeff := coeff(oderhs, ind_var, d-ord);
    end if; 
    # Then we assign all initial conditions to the first a[ord-1] terms.
    if d < ord then
        return init_cond[d]/d!*ind_var^d;
    else
        # This is the single summation.
        local s1 := add(coeff(hpoly[d-ord+i],ind_var,d-ord+i)*coeffs[i][0]*(d-ord+i)!/(d-ord)!, 
                                i=0 .. ord-1);
        local s2 := 0;

        # There is a special condition when min = 0, in this case we ignore the
        # double summation and use single sum we have, similar to case of
        # constant coefficients
        if min(max_degree,d-ord) <> 0 then
            # This is the double summation.
            s2 := add(add(coeff(hpoly[d-ord-k+i],ind_var,d-ord-k+i)*coeffs[i][k]*(d-ord-k+i)!/(d-ord-k)!,
                            i=0..ord),k=1.. min(max_degree,d-ord));
        end if;
        # Final formula all together
        local total := (d-ord)!/(d!*constant_coeff)*(oderhs_coeff-s1-s2)*ind_var^(d);
        
        return total;
    end if; 
end proc;

# NOTE: This code is not in p4 of maple because it is not
# working. I will comment it.
(*
local solution_around_exp_pnt_gen := proc(_self :: PowerSeriesObject,
                                          d :: nonnegint)
    local edo := _self:-ancestors:-ode;
    local ind_var := edo:-GetIndependentVariable(edo);
    local dep_var := edo:-GetDependentVariable(edo);
    local init_cond_sym := edo:-GetInitialConditions(edo);
    local init_cond := edo:-ConvertConstantSymbolsToD(edo, init_cond_sym);
    local coefficients := edo:-GetCoeffs(edo);
    local ord := edo:-GetOrder(edo);
    local oderhs := edo:-GetRHS(edo);
    local exp_pnt := edo:-GetExpansionPoint(edo);

    ASSERT(numelems(_self:-vars)=1, "incorrect number of variables");
    # Note: this is a frozen variable obtained with the freeze command.
    # Note: all the computations must be done using this variable that
    # represents ind_var-exp_pnt. The power series MUST contain monomials
    # in var! It will take care of printing in the proper format to the 
    # user, ie, before printing, the PSO substitutes var=thaw(var) if
    # necessary.
    local var := _self:-vars[1];

    local ad_at_x0 := eval(coefficients[ord],[ind_var = exp_pnt]);

    ASSERT(ad_at_x0 <> 0);

    local equation := _self:-ancestors:-diff_eq;
    local Y := _self:-ancestors:-Y;

    if d < ord then 
        Y ,= init_cond[d]/d!;
        ASSERT(numelems(Y)=d+1);
        print(Y[d]);
        return Y[d]*(var)^d;
    else
        #check id we need to diff d-ord times
        
        local j;
        local R := equation - coefficients[ord]*(D@@ord+1)(ind_var)(dep_var);
        print("HI1", R);
        print(numelems(Y),d);
        local u := eval(R, [seq((D@@j)(dep_var)(ind_var)=Y[j],j=0..d-1)]);
        print("HI2", [seq((D@@j)(dep_var)(ind_var)=Y[j],j=0..d-1)]);
        u := normal(u/ad_at_x0);
        Y ,= u;
        ASSERT(numelems(Y)=d+1);
        print("HI3");
        equation := diff(equation,ind_var);
        return (u/d!)*(ind_var-exp_pnt)^d;
    end if;
end proc;
*)

## Power series from ode constructors.
# Function for ODE solution with constant coefficients.
local from_linear_coefficient_ode::static := proc(ind_var::name, 
                                                    dep_var::name, ode::set(equation), 
                                                    {exp_pnt :: {integer, name} := 0}, 
                                                    ode_obj_in::OdeManipulator:-OdeObject:=NULL, $)
    # NOTE: if exp_pnt is different than 0, OdeManipulator applies 
    # the change of variable [ind_var=ind_var-exp_pnt] to ode.
    local ode_obj;

    if ode_obj_in=NULL then
        ode_obj := OdeManipulator:-OdeObject(ind_var, dep_var, ode, _options['exp_pnt']);
    else
        ode_obj := ode_obj_in;
    end if;

    ASSERT(ode_obj:-GetOdeType(ode_obj)="lin_cons_coeff");

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
                                                      {exp_pnt :: {integer, name} := 0}, 
                                                      ode_obj_in::OdeManipulator:-OdeObject:=NULL, $)

    local ode_obj;

    if ode_obj_in=NULL then
        ode_obj := OdeManipulator:-OdeObject(ind_var,dep_var,ode, _options['exp_pnt']);
    else
        ode_obj := ode_obj_in;
    end if;

    ASSERT(ode_obj:-GetOdeType(ode_obj)="lin_pol_coeff");

    local deg := 0;
    local init_cond := ode_obj:-GetInitialConditions(ode_obj);
    local var := ode_obj:-GetIndependentVariable(ode_obj);
    local ord := ode_obj:-GetOrder(ode_obj);
    local constant_coeff := ode_obj:-GetCoeffs(ode_obj)[ord][0];

    local a0_sym := init_cond[0];
    local a0 := ode_obj:-ConvertConstantSymbolsToD(ode_obj, a0_sym);
    local hpoly := Array(0..0, [a0]);

    # For our algorithm, the leading polynomial must be invertible
    # So we throw an error if that's not the case.
    if constant_coeff = 0 then
        error "leading polynomial not invertible";
    end if;

    return Object(PowerSeriesObject, hpoly, deg, poly_coeff_ode_gen, {var},
          ["ode" = ode_obj], undefined);

end proc:

# NOTE: This code is not in p4 of maple because it is not
# working. I will comment it.
(*
local solution_around_exp_pnt := proc(ind_var::name, 
                                      dep_var::name, ode::set(equation), 
                                      {exp_pnt :: {integer, name} := 0}, 
                                      ode_obj_in::OdeManipulator:-OdeObject:=NULL, $)

    
   
    local ode_obj;

    if ode_obj_in=NULL then
        ode_obj := OdeManipulator:-OdeObject(ind_var, dep_var, ode, _options['exp_pnt']);
    else
        ode_obj := ode_obj_in;
    end if;

    ASSERT(ode_obj:-GetOdeType(ode_obj)="lin_func_coeff");

    local deg := 0;

    local init_cond := ode_obj:-GetInitialConditions(ode_obj);
    local ord := ode_obj:-GetOrder(ode_obj);
    
    local diff_eq := ode_obj:-GetOdeEquation(ode_obj); #TODO: DOuble check rhs is there w/ example

    local a0_sym := init_cond[0];
    local a0 := ode_obj:-ConvertConstantSymbolsToD(ode_obj, a0_sym);
    local Y := Array(0..0,[a0]);
    local hpoly := Array(0..0, [a0]);

    return Object(PowerSeriesObject, hpoly, deg, solution_around_exp_pnt_gen, {freeze(ind_var-exp_pnt)},
          ["ode" = ode_obj, "Y" = Y, "diff_eq" = diff_eq], undefined);
end proc:
*)

## Dummy eval.
local my_eval::static := proc()
    if _npassed=1 then
        return eval(_passed);
    end if; 

    if _npassed<1 or _npassed>2 then
        error "my_eval requires 2 parameters but received [%1]", _passed;
    end if;

    local pso_dummy := _passed[1]; 

    local old := kernelopts( 'opaquemodules' = 'false' );
    if not type(pso_dummy, 'MultivariatePowerSeries:-OdeManipulator:-PowerSeriesWrapper') then
        return eval(pso_dummy, _passed[2]);
    end if;
    kernelopts( 'opaquemodules' = old );

    local q := _passed[2];
    if not type(q, equation) then
        error "%1 should be of type equation", q;
    end;

    if rhs(q)<>0 or (not pso_dummy:-pso:-IsUnit(pso_dummy:-pso)) then 
        return FAIL; 
    end if;

    return pso_dummy:-pso:-HomogeneousPart(pso_dummy:-pso, 0);
end proc;


# Wrapper for all the ODE solution algorithms.
export ode_sol_wrapper::static := proc(ind_var::name, 
                                        dep_var::name, ode::set(equation), 
                                        {exp_pnt :: {integer, name} := 0}, $)
    
    local ode_obj := OdeManipulator:-OdeObject(ind_var, dep_var, ode, _options['exp_pnt']);
    
    # NOTE: "lin_cons_coeff" is the type of linear ODEs with constant coefficient 
    #       (for exp_pnt=0).
    # NOTE: "lin_pol_coeff" is the type of linear ODEs with polynomial coefficient.
    #       (for exp_pnt=0).
    # NOTE: "lin_func_coeff" is the type of linear ODEs with function coefficient.
    if ode_obj:-GetOdeType(ode_obj)="lin_cons_coeff" then
        return PowerSeriesObject:-from_linear_coefficient_ode(_passed, ode_obj);
    elif ode_obj:-GetOdeType(ode_obj)="lin_pol_coeff" then
        return PowerSeriesObject:-from_polynom_coefficient_ode(_passed, ode_obj);
    else 
        error "Ode not supported";
    end if;

end proc;            