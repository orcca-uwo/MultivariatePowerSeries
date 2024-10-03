#test
##############
## ODEtests ##
##############
with(TestTools):
# TODO: FIX TEST ENUMERATION

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ODEO := MultivariatePowerSeries:-OdeManipulator:-OdeObject:
kernelopts(opaquemodules=true):

local MyIsEqual := proc(A::Array(Array), B::Array(Array), $)
    local n := numelems(A);

    if [rtable_dims(A)]<>[rtable_dims(B)] then
        return false;
    end if; 

    return andseq(ArrayTools:-IsEqual(A[i], B[i]), i=0..n-1);
end proc:

#Test 1: Determining if it prints to correct ODE	
Try[testnoerror]("test 1", ODEO(t, x,{diff(x(t),t)=x(t), x(0)=1}), 'assign'='eOde');
Try[testnoerror]("test 2", ODEO(t, x,[diff(x(t),t)=x(t), x(0)=1]));
Try[testnoerror]("test 3",eOde,{x(0) = 1, D(x)(t) = x(t)});
my_coeff :=  Array(0..0,[1]):
Try("test 4", ArrayTools:-IsEqual(eOde:-GetInitialConditions(eOde), my_coeff), true); ## write for all ODEs.
my_coeff :=  Array(0..1,[-1, 1]):
Try("test 5", ArrayTools:-IsEqual(eOde:-GetCoeffs(eOde), my_coeff), true);
Try("test 6", eOde:-GetOrder(eOde), 1);


#Test 2: Determining if it prints the correct ODE
Try[testnoerror]("test 7", ODEO(t,x,{diff(x(t),t,t)+x(t)=0,x(0)=0,D(x)(0)=1}), 'assign'='SinOde');
Try[testnoerror]("test 8", ODEO(t,x,[diff(x(t),t,t)+x(t)=0,x(0)=0,D(x)(0)=1]), 'assign'='SinOde');
my_coeff :=  Array(0..1,[0,1]):
Try("test 9", ArrayTools:-IsEqual(SinOde:-GetInitialConditions(SinOde), my_coeff), true);
my_coeff := Array(0..2,[1,0,1]):
Try("test 10",ArrayTools:-IsEqual(SinOde:-GetCoeffs(SinOde), my_coeff),true);
Try("test 11",SinOde:-GetOrder(SinOde),2);

#Test 3: Determining if it prints the correct ODE
Try[testnoerror]("test 12", ODEO(t,x,{diff(x(t),t,t)+x(t)=0, x(2*pi)=1,D(x)(0)=0}), 'assign'='CosOde');
Try[testnoerror]("test 13", ODEO(t,x,[diff(x(t),t,t)+x(t)=0, x(2*pi)=1,D(x)(0)=0]), 'assign'='CosOde');
my_coeff :=  Array(0..1,[1, 0]):
Try("test 14",ArrayTools:-IsEqual(CosOde:-GetInitialConditions(CosOde), my_coeff), true);
my_coeff :=  Array(0..2,[1,0,1]):
Try("test 15",ArrayTools:-IsEqual(CosOde:-GetCoeffs(CosOde), my_coeff), true);
Try("test 16",CosOde:-GetOrder(CosOde),2);

#Test 4: Determining if it prints the correct ODE
Try[testnoerror]("test 17", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0, x(0)=0, D(x)(0)=1, (D@@2)(x)(0)=1}), 'assign'='Ode4');
Try[testnoerror]("test 18", ODEO(t,x,[diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0, x(0)=0, D(x)(0)=1, (D@@2)(x)(0)=1]), 'assign'='Ode4');
my_coeff :=  Array(0..2,[0, 1, 1]):
Try("test 19", ArrayTools:-IsEqual(Ode4:-GetInitialConditions(Ode4), my_coeff), true);
my_coeff := Array(0..3,[0, -20, -1, 1]):
Try("test 20",ArrayTools:-IsEqual(Ode4:-GetCoeffs(Ode4), my_coeff), true);
Try("test 21",Ode4:-GetOrder(Ode4),3);

#Test 5: Determining if it prints the correct ODE
Try[testnoerror]("test 22", ODEO(t, x,{diff(x(t),t)-2*t*x(t)=0, x(0)=1}), 'assign'='ex2Ode');
Try[testnoerror]("test 23", ODEO(t, x,[diff(x(t),t)-2*t*x(t)=0, x(0)=1]), 'assign'='ex2Ode');
my_coeff :=  Array(0..0,[1]):
Try("test 24", ArrayTools:-IsEqual(ex2Ode:-GetInitialConditions(ex2Ode), my_coeff), true); ## write for all ODEs.

my_coeff :=  Array(0.. 1, [Array(0.. 1, [0, -2]), Array(0.. 1, [1, 0])]):
Try("test 25", MyIsEqual(ex2Ode:-GetCoeffs(ex2Ode), my_coeff), true);
Try("test 26", ex2Ode:-GetOrder(ex2Ode), 1);

Try[testnoerror]("test 27", ODEO(t, x,{diff(x(t),t,t)-t*diff(x(t),t)+4*x(t)=0, x(0)=1,D(x)(0)=0}), 'assign'='polyOde');
Try[testnoerror]("test 28", ODEO(t, x,[diff(x(t),t,t)-t*diff(x(t),t)+4*x(t)=0, x(0)=1,D(x)(0)=0]), 'assign'='polyOde');
my_coeff :=  Array(0..1,[1, 0]):
Try("test 29", ArrayTools:-IsEqual(polyOde:-GetInitialConditions(polyOde), my_coeff), true); ## write for all ODEs.

my_coeff :=  Array(0..2,[Array(0.. 1, [4, 0]), Array(0.. 1, [0, -1]), Array(0.. 1, [1, 0])]):
Try("test 30", MyIsEqual(polyOde:-GetCoeffs(polyOde), my_coeff), true);

Try("test 31", polyOde:-GetOrder(polyOde), 2);

# Test errors.
Try[testerror]("test 32", ODEO(x,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0, D(x)(0)=1, x(0)=0}), 'assign'='my_edo');
Try[testerror]("test 33", ODEO(x,x,{D(x)(0)=1, x(0)=0}), 'assign'='my_edo');

#Test for ode without initial conditions.
Try[testnoerror]("test 34", ODEO(y, x,{diff(x(y),y)=2*y*x(y)}), 'assign'='my_edo');
my_coeff :=  Array(0..0,[x(0)]):
P:=my_edo:-GetInitialConditions(my_edo):
Try("test 35", ArrayTools:-IsEqual(my_edo:-ConvertConstantSymbolsToD(my_edo, P), my_coeff), true);

Try[testnoerror]("test 36", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0}), 'assign'='my_edo');
my_coeff :=  Array(0..2,[x(0), (D)(x)(0), (D@@2)(x)(0)]):
P:=my_edo:-GetInitialConditions(my_edo):
Try("test 37", ArrayTools:-IsEqual(my_edo:-ConvertConstantSymbolsToD(my_edo, P), my_coeff), true);

Try[testnoerror]("test 38", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0, D(x)(0)=1, x(0)=0}), 'assign'='my_edo');
my_coeff :=  Array(0..2,[0, 1, (D@@2)(x)(0)]):
P:=my_edo:-GetInitialConditions(my_edo):
Try("test 39", ArrayTools:-IsEqual(my_edo:-ConvertConstantSymbolsToD(my_edo, P), my_coeff), true);

Try[testnoerror]("test 40", ODEO(x,y,{diff(y(x),x,x,x)-diff(y(x),x,x)-20*diff(y(x),x)=0, D(y)(0)=1, (D@@2)(y)(0)=1}), 'assign'='my_edo');
my_coeff :=  Array(0..2,[y(0), 1, 1]):
P:=my_edo:-GetInitialConditions(my_edo):
Try("test 41", ArrayTools:-IsEqual(my_edo:-ConvertConstantSymbolsToD(my_edo, P), my_coeff), true);

#Test change order initial conditions.
Try[testnoerror]("test 42", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0, D(x)(0)=1, (D@@2)(x)(0)=1, x(0)=0}), 'assign'='Ode4');
my_coeff :=  Array(0..2,[0, 1, 1]):
Try("test 43", ArrayTools:-IsEqual(Ode4:-GetInitialConditions(Ode4), my_coeff), true);

#Test expansion point.
Try[testnoerror]("test 44", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0}, exp_pnt=1), 'assign'='my_edo');
Try("test 45", my_edo:-GetExpansionPoint(my_edo), 1);
Try[testnoerror]("test 46", ODEO(t,x,{diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0}, exp_pnt=a), 'assign'='my_edo');
Try("test 47", my_edo:-GetExpansionPoint(my_edo), a);

#More tests.
de := diff(x(t),t$2) + omega^2 * x(t) = 0;
Try[testnoerror]("test 48", ODEO(t,x,{de}), 'assign'='my_de');
my_coeff :=  Array(0..1,[x(0), (D)(x)(0)]):
P:=my_de:-GetInitialConditions(my_de):
Try("test 49", ArrayTools:-IsEqual(my_de:-ConvertConstantSymbolsToD(my_de, P), my_coeff), true);
my_coeff := Array(0..2,[omega^2,0,1]):
Try("test 50",ArrayTools:-IsEqual(my_de:-GetCoeffs(my_de), my_coeff),true);
Try("test 51",my_de:-GetOrder(my_de),2);

#Test a non-linear ode gives an error.
non_linear_de := diff(x(t),t$2) + omega^2 * x(t) = diff(y(t),t$3);
Try[testerror]("test 52", ODEO(t,x,{non_linear_de}));

non_linear_de := diff(x(t),t$2) + x(t)^2 = diff(y(t),t$3);
Try[testerror]("test 53", ODEO(t,x,{non_linear_de}));

de_not_supported_yet := diff(x(t),t)+sin(x(t))=0;
Try[testerror]("test 54", ODEO(t,x,{de_not_supported_yet}));

#end test
