#OdeObject

# TODO: expansion around a point different than 0.
# TODO: for ODEs with PowerSeries coefficients, the printing 
#       of the object is ugly. This should be improved.

OdeObject := module()
option object;

    # OdeObject Module Attributes
    local
        # ode: the ode.
        ode::algebraic, 
        # equ : the set of the ode and the initial conditions (used to print the object).
        equ::string,
        # coefs : An Array which stores all coefficients of the ODE 
        coefs::Array,
        # ord : a nonnegint which is the order of the ODE
        ord::nonnegint,
        # ind_var : the independent variable of the ODE. 
        #           There may be ways of finding this, but for the sake of simplicity
        #           this is easier.
        ind_var::name,
        #dep_var : the dependent variable 
        #         (a.k.a the name of the function we are trying to find)
        dep_var::name,
        # init_cond : the array of initial conditions.
        init_cond::Array,
        # init_cond_eq : the array of initial conditions equations.
        init_cond_eq::Array,
        # exp_pnt : the expansion is done around this point.
        exp_pnt::simple_algebraic,
        # ode_rhs : We assume the rhs of our ode is a polynomial of the 
        #         independent variable with no dependence on the dependent variable
        ode_rhs :: polynom,
        # ode_type : flag used to identify the type of ode.
        ode_type :: identical("lin_cons_coeff", "lin_pol_coeff", "lin_func_coeff", "lin_pso_coeff"),
        # dstyle : a display style for a specific odes  
        dstyle :: PSDS_TYPE; 



# ModuleApply: to make OdeObject an object factory
local ModuleApply ::static := proc()
    Object(OdeObject, _passed);
end proc;

local UnifyInput::static := proc(sys :: { list, set, thistype }( { algebraic, equation } ), $)
  if sys :: ':-list' then
    return thisproc( convert( sys, ':-set' ));
  elif not sys :: { ':-set', ':-list' } then
    return thisproc( [ sys ] );
  elif membertype( ':-algebraic', sys ) then
    return thisproc( map( a -> ifelse( a :: ':-equation', a, a = 0 ), sys )); 
  end if;

  return sys;  
end proc;

local CheckIC::static := proc(ic::equation, ind_var::And(name,Not(constant)), 
                                dep_var::And(name,Not(constant)), exp_pnt::{integer, name} ,$)
    
    local rhs_type := ':-And'( ':-simple_algebraic', ':-freeof'( dep_var ) );
    local lhs_type := { typefunc( 'identical'(exp_pnt), 'identical'(dep_var) ), typefunc( 'identical'(exp_pnt), typefunc( 'identical'(dep_var), { 'identical(D) @@ posint', 'identical(D)' } ) ) };

    local ic_rhs := rhs(ic);
    if not type(ic_rhs, rhs_type) then
        return false;
    end if;

    local ic_lhs := lhs(ic);
    if not type(ic_lhs, ':-dependent'(dep_var)) then
        return false;
    end if;

    return true;
end proc;

# To check ode is linear.
local ClassifyODE::static := proc( 
                            de :: { thistype, equation }( simple_algebraic ),
                            ind_var :: And( name, Not( constant ) ),
                            dep_var :: And( name, Not( constant ) ),
                            de_order::nonnegint, $)

    if de :: ':-equation' then
        return thisproc( ( lhs - rhs )( de ), _passed[ 2 .. ] );
    end if; 
    ASSERT( de :: ':-algebraic' );

    # Type for the dependent variable function and its derivatives, write in terms of 'D'.
    local T1 := ':-typefunc'( ':-identical'( ind_var ), ':-identical'( dep_var ) );
    local T2 := ':-typefunc'( ':-identical'( ind_var ), ':-typefunc'( ':-identical'( dep_var ), ':-identical'( ':-D' ) ) );
    local T3 := ':-typefunc'( ':-identical'( ind_var ), ':-typefunc'( ':-identical'( dep_var ), ':-identical'( ':-D' ) @@ ':-posint' ) );
    local T := { T1, T2, T3 };

    # Determine if the ODE is linear. 
    local f := convert( de, ':-D' );  
    local i :: ':-nothing', z :: ':-nothing';
    local B := [ seq( z[ i ], i = 0 .. de_order ) ];
    local g := subs( seq( `@@`( D, i )( dep_var )( ind_var ) = z[ i ], i = de_order .. 0, -1 ), f );

    if not g :: ':-linear'( B ) then
        return ':-false';
    end if;
    
    # Determine if the coefficients are polynomial in the independent variable.
    local C := PolynomialTools:-CoefficientList( g, [ seq( B ) ] );

    if andmap( type, C, ':-polynom'( ':-anything', ind_var ) ) then
        return "lin_pol_coeff";
    end if;

    local C_unapply := map(unapply, C, ind_var);
    if andmap( type, C_unapply, ':-mathfunc'( ':-anything', ind_var ) ) then
        return "lin_func_coeff";
    else 
        return false;
    end if;
    
end proc;

local CheckLinearPolynomialODE :: static := proc(ind_var :: And(name,Not(constant)), dep_var :: And(name,Not(constant)), exp_pnt::{integer, name},
                                    ics, dec_ord, dec_D_in, derivatives, de, 
                                    ode_type_in::identical("lin_pol_coeff", "lin_func_coeff"), $)
    local dec_D := dec_D_in;
    local ics_array, ic, ics_eq;

    if numelems(ics)< dec_ord then
        # We create dec_ord initial conditions.
        local ics_dummy := [seq((':-D'@@ ic)(dep_var)(exp_pnt)=`tools/gensym`(cat(dep_var,ic)), ic in 0..dec_ord-1)];
        local ics_to_add := remove( a -> membertype( 'identical'(lhs( a )) = 'anything', ics ), ics_dummy);
        local ics_to_add_set := convert(ics_to_add, ':-set');
        local new_ics := convert(`union`(ics, ics_to_add_set), ':-list');
        ics_array := Array(0..dec_ord-1, map(rhs, new_ics));
        dec_D := `union`(dec_D, ics_to_add_set);
        ics_eq := Array(map(rhs=lhs,ics_dummy));
    else
        # The rhs of the initial conditions as array.
        
        ics_array := Array(0..(numelems(ics)-1), [seq(rhs(ic), ic in ics)]);
        ics_eq := Array(0..-1);
    end if;

    local oderhs_dummy := -remove(has,de,convert(derivatives,':-set'));
    if type(oderhs_dummy,polynom(anything,ind_var)) = false then
        error "expecting the right-hand side of the ODE to be a polynomial of the independent variable, but received", oderhs_dummy
    end if;

    local derivative;
    local coefficients := Array(0..numelems(derivatives)-1,[seq(coeff(de, derivative), derivative in derivatives)]);
        
    local ode_type := "lin_cons_coeff";
    # If statement to determine whether we are in the case of polynomials 
    # or constant coefficients
    if 'exp_pnt'=0 and ode_type_in="lin_pol_coeff" and  depends(coefficients, ind_var) then
        local i; 
        # computes the coefficients of the polynomials
        local coeff_poly := map(PolynomialTools:-CoefficientList,coefficients,ind_var);
        # finds max degree of the polynomial coefficients
        local max_degree := max(map(degree,coefficients,ind_var));
        # makes all coefficients of the polynomials of the same size
        # (a.k.a inserts 0s in place of the poly in degrees higher than the initial poly)
        coefficients := map2(Array,0..max_degree,coeff_poly);
        ode_type := "lin_pol_coeff";
    elif ode_type_in ="lin_func_coeff" then 
        ode_type := "lin_func_coeff";
    end if;

    return [dec_D, dec_ord, coefficients, 
            ics_array, ics_eq, oderhs_dummy, ode_type, de];
end proc:


# We decompose the ode wrt D and return all we need for constructors
# CheckSystem.
local CheckSystem::static := proc(ind_var :: And(name,Not(constant)), dep_var :: And(name,Not(constant)), exp_pnt:: simple_algebraic, 
                                    sys :: { list, set, thistype }( { algebraic, equation } ), $)
    local ode := OdeObject:-UnifyInput(sys);

    # Rewrite ODE in terms of diff operator D.
    local dec_D := convert(ode, ':-D');
    local T := { typefunc( 'identical'(ind_var), 'identical'(dep_var) ), 
                typefunc( 'identical'(ind_var), typefunc( 'identical'(dep_var), { 'identical(D) @@ posint', 'identical(D)' } ) )};
    local funs := indets(dec_D, T);
 
    local op0s := map2(op, 0, funs);
    local Dfs := select(type, op0s, typefunc(identical(':-D')));                 
    local Dnfs := select(type, op0s, typefunc(identical(':-D') @@ posint));
    
    # Determines the order of the ODE.
    local dec_ord := max(map2(op, [0, 2], Dnfs) union map(1, Dfs));

    if dec_ord = 0 then
        error "expecting one or more derivatives of %1 to occur in the ordinary differential equation", dep_var(ind_var);
    end if;

    local n;
    # Makes a list of the dependent variable and its derivatives.
    local derivatives;

    derivatives := [seq((':-D' @@ n)(dep_var)(ind_var), n = 0 .. dec_ord)];

    # We divide the input in the ode and the initial conditions.
    local (des, ics) := selectremove(has, dec_D, convert(derivatives, ':-set'));

    if numelems(des) <> 1 then
        error "currently the command only supports one ODE but received: %1", des;
    end if; 

    # We check the initial conditions.
    local ic;
    if not andseq(CheckIC(ic, ind_var, dep_var, exp_pnt), ic in ics) then 
        error "the format of the initial conditions %1 is incorrect", ics;
    end if;

    # The differential equation.
    local de := op(des); 

    if type(de, equation) = true then
        de := (lhs-rhs)(de);
    end if;

    local chk := ClassifyODE(de, ind_var, dep_var, dec_ord);        
    if chk=false then 
        error "expecting one linear ODE but received: %1", de;
    else
        return CheckLinearPolynomialODE(ind_var, dep_var, exp_pnt, 
                                    ics, dec_ord, dec_D, derivatives, de, chk);
    end if; 
   
end proc:

local pow_ser_ode_decom::static := proc(f :: uneval, p :: name, $ )
    local w, i, t;
    local T := [ `+`, `*`, `^` ];
    local g := f;

    for i, t in T do
        g := subsindets( g, t, z -> w[i]( op(z) ) );
    end do;

    g := eval( g );
    local old := kernelopts( 'opaquemodules' = 'false' );
    local Q := convert( indets( g, 'MultivariatePowerSeries:-PowerSeriesObject' ), 'list' );
    kernelopts( 'opaquemodules' = old );

    local j :: 'nothing';
    local R := [ seq( p[j] = Q[j], j = 1 .. numelems( Q ) ) ];
    local R2 := [ seq( Q[j] = Object(PowerSeriesWrapper, Q[j]), j = 1 .. numelems( Q ) ) ];
    g := subs( seq( R2 ), g );
    g := eval(g, [seq( w[j] = T[j], j = 1 .. numelems( T ) )] );

    return [g, R];
end proc;

# To handle ODEs with PSO coefficients.
export DummyConstructor::static := proc(ind_var :: And(name,Not(constant)),
                                    dep_var :: And(name,Not(constant)),
                                    ode::uneval, 
                                    {exp_pnt :: {integer, name} := 0}, 
                                    dstyle :: PSDS_TYPE := [], $)

    local p;
    local sys := OdeObject:-pow_ser_ode_decom(ode, p);

    if sys[2]<>[] then
        return Object(OdeObject, ind_var, dep_var, sys[1], _options['exp_pnt'], sys[2], dstyle);
    else 
        return Object(OdeObject, ind_var, dep_var, ode, _options['exp_pnt'], [], dstyle);
    end if;
end proc;

# ModuleCopy: is called when a new object of the same class is created via a call to Object
local ModuleCopy ::static := proc(new :: OdeObject,
                                    old :: OdeObject,
                                    ind_var :: And(name,Not(constant)),
                                    dep_var :: And(name,Not(constant)),
                                    ode::uneval, 
                                    {exp_pnt :: {integer, name} := 0}, 
                                    pso_equ::list(equation):=[],
                                    dstyle :: PSDS_TYPE := [], $)
    if _npassed > 2 then
        if ind_var=dep_var then 
            error "expecting independent and dependent variable names to be different, but received %1 , %2", ind_var, dep_var;
        end if;

        local result := OdeObject:-CheckSystem(ind_var, dep_var, exp_pnt, ode);

        if pso_equ<>[] then
            new:-ode_type := "lin_pso_coeff";
        else 
            ASSERT(result[7]::identical("lin_cons_coeff", "lin_pol_coeff", "lin_func_coeff", "lin_pso_coeff"));
            new:-ode_type := result[7];
        end if;

        # Note: in the case "lin_pso_coeff" ode does not have pso coefficients.
        new:-ode := result[8];
        new:-equ := convert(result[1], ':-string');
        new:-ord := result[2];
        new:-coefs := result[3];
        new:-init_cond := result[4];
        new:-init_cond_eq := result[5];
        new:-ode_rhs := result[6];

        new:-ind_var := ind_var;
        new:-dep_var := dep_var;
        new:-exp_pnt := exp_pnt;
        new:-dstyle := dstyle;
    elif _npassed= 2 then
        new:-ode := old:-ode;
        new:-equ := old:-equ;
        new:-ord := old:-ord;
        new:-coefs := copy(old:-coefs);
        new:-init_cond := copy(old:-init_cond);
        new:-init_cond_eq := copy(old:-init_cond_eq);
        
        new:-ind_var := old:-ind_var;
        new:-dep_var := old:-dep_var;
        new:-exp_pnt := old:-exp_pnt;
        new:-ode_rhs := old:-ode_rhs;
        new:-ode_type := old:-ode_type;
        new:-dstyle := old:-dstyle;
    else
        error "you cannot copy the original OdeObject object";
    end if;

end proc;
    
local defaultDisplayStyle ::static := [];

# SetDefaultDisplayStyle
# To set the local global variable defaultDisplayStyle 
export SetDefaultDisplayStyle ::static := proc(s :: PSDS_TYPE, $)
    defaultDisplayStyle := s;
end proc;

# GetDefaultDisplayStyle
export GetDefaultDisplayStyle ::static := proc($)
    return defaultDisplayStyle;
end proc;

# SetDisplayStyle
# to set the dstyle of the input object
export SetDisplayStyle ::static := proc(obj :: OdeObject, s :: PSDS_TYPE, $)
    obj:-dstyle := s; 
end proc;

# Display
# to display the input object, 
# user_dstyle : will overwrite the default and _self:-dstyle style 
# output : (doesn't need to be documented) used for internal purposes
export Display ::static := proc(_self :: OdeObject, 
                                 user_dstyle :: list := [],
                                 output :: identical("full", "typeset", "fullstring", "string")
                                 := "full",
                                 {updateterms :: name := NULL},
                                 $)
local mystyle;
    if user_dstyle <> [] then
        mystyle := user_dstyle;
    elif _self:-dstyle <> [] then
        mystyle := _self:-dstyle;
    else
        mystyle := OdeObject:-GetDefaultDisplayStyle();
    end if;

    return _self:-equ;
end proc;

# ModulePrint: print _self :: OdeObject 
local ModulePrint :: static := proc(_self :: OdeObject, $)
    if IsWorksheetInterface() then
        return Display(_self);
    else
        return convert(Display(_self), ':-name');
    end if;
end proc;

# Module Deconstruct
local ModuleDeconstruct :: static := proc(_self :: OdeObject, $)
    return 'Object'('OdeObject', _self:-ode, _self:-equ, _self:-ord, 
                        _self:-ind_var, _self:-dep_var, _self:-init_cond,
                        _self:-ode_rhs, _self:-ode_type);
end proc;


$include "MultivariatePowerSeries/OdeManipulator/OdeObject/src/getters.mm"
end module:
