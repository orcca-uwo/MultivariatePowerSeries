#test
##############
## ODEtests ##
##############
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
ODEO := MultivariatePowerSeries:-OdeManipulator:-OdeObject:
DummyConstructor := MultivariatePowerSeries:-OdeManipulator:-OdeObject:-DummyConstructor;
kernelopts(opaquemodules=true):

local MyIsEqual := proc(A::Array(Array), B::Array(Array), $)
    local n := numelems(A);

    if [rtable_dims(A)]<>[rtable_dims(B)] then
        return false;
    end if; 

    return andseq(ArrayTools:-IsEqual(A[i], B[i]), i=0..n-1);
end proc:

# Tests: ode with function coefficient.
Try[testnoerror]("test func coeff 1", ODEO(t, x,{diff(x(t),t)+sin(t)*x(t)=0, x(0)=1}), 'assign'='my_ode');
Try[testnoerror]("test func coeff 2", ODEO(t,x,{sin(t)*diff(x(t),t,t)-2*sin(t)*diff(x(t),t)-cos(t)*x(t)= t^2}),'assign'='my_ode2');
my_coeffs := Array(0..2,[-cos(t),-2*sin(t),sin(t)]);
Try("test 3", ArrayTools:-IsEqual(my_ode2:-GetCoeffs(my_ode2),my_coeffs),true);

#Testing if order matters
Try[testnoerror]("test func coeff 3", ODEO(t,x,{-2*sin(t)*diff(x(t),t)= -cos(t)*diff(x(t),t,t)+ cos(t)*x(t)+ t^2}),'assign'='my_ode2');
my_coeffs := Array(0..2,[-cos(t),-2*sin(t),cos(t)]);
Try("test 4", ArrayTools:-IsEqual(my_ode2:-GetCoeffs(my_ode2),my_coeffs),true);

#Testing if we can reassign coefficients
a(t) := sin(t);
Try[testnoerror]("test func coeff 4", ODEO(t,x,{a(t)*diff(x(t),t,t)+ b(t)*diff(x(t),t)+ c(t)*x(t)= 0}), 'assign'='my_ode3');
my_coeffs := Array(0..2,[c(t),b(t),a(t)]);
Try("test 5", ArrayTools:-IsEqual(my_ode3:-GetCoeffs(my_ode3),my_coeffs),true);

# Testing if we can apply derivative rules to get linear ODE
Try[testnoerror]("test func coeff 5", ODEO(t,x,{diff(sin(t)*x(t),t)=0}),'assign'= 'my_ode4');
my_coeffs := Array(0..1,[cos(t),sin(t)]);
Try("test 6", ArrayTools:-IsEqual(my_ode4:-GetCoeffs(my_ode4),my_coeffs),true);

# Tests: ode with pso coefficient.
Try[testnoerror]("test pso coeff 1", DummyConstructor(t, x,{PowerSeries(1)*diff(x(t),t)+PowerSeries(t)*x(t)=0, x(0)=1}), 'assign'='my_ode');
Try[testnoerror]("test pso coeff 2", DummyConstructor(t, x,{diff(x(t),t)+PowerSeries(t)*x(t)=0, x(0)=1}), 'assign'='my_ode');

#end test
