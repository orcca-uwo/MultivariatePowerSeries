#Gets.

# To get the initial conditions of our ODE.
export ConvertConstantSymbolsToD::static := proc(_self :: OdeObject, exp, $)
    if numelems(_self:-init_cond_eq)=0 then 
        return exp;
    end if; 
    
    local exp_new := subs(seq( _self:-init_cond_eq ), exp);
    
    return exp_new; 
end proc;

# To get the initial conditions of our ODE.
export GetInitialConditions::static := proc(_self :: OdeObject, $)
    return _self:-init_cond; 
end proc;

# To get the independent variable of our ODE.
export GetIndependentVariable::static := proc(_self :: OdeObject, $)
    return _self:-ind_var; 
end proc;

# To get the dependent variable of our ODE.
export GetDependentVariable :: static := proc(_self :: OdeObject,$)
    return _self:-dep_var;
end proc;

# To get the order of our ODE.
export GetOrder :: static := proc(_self :: OdeObject,$)
    return _self:-ord;
end proc;

# To get the ode.
export GetOde :: static := proc(_self :: OdeObject,$)
    return _self:-equ;
end proc;

# To get the coefficients of our ODE.
export GetCoeffs :: static := proc(_self :: OdeObject,$)
    return _self:-coefs;
end proc;

# To get the expansion point of our ODE.
export GetExpansionPoint :: static := proc(_self :: OdeObject,$)
    return _self:-exp_pnt;
end proc;





