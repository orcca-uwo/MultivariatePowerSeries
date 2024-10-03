#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
from_polynom_coefficient_ode := MultivariatePowerSeries:-PowerSeriesObject:-from_polynom_coefficient_ode:
kernelopts(opaquemodules=true):
#TODO: change order
ord := 10:

#Test on ode1
my_ode:= {diff(x(t),t)- 2*t*x(t)= 0, x(0)=1}:

Try[testnoerror]("test 1", from_polynom_coefficient_ode(t, x, my_ode), 'assign'='ps1');
sol1 := dsolve(my_ode, x(t), 'type=series', 'order' = ord):
sol_pol := convert( sol1, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 2", Truncate(ps1, ord-1), sol_pol);

#Test on ode2
my_ode := {diff(x(t),t,t)-t^2*x(t)=0}:

Try[testnoerror]("test 3", from_polynom_coefficient_ode(t, x, my_ode), 'assign'='ps2');
sol2 := dsolve(my_ode, x(t), 'type=series', 'order' = ord):
sol_pol := convert( sol2, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 4", Truncate(ps2, ord-1), sol_pol);

#Test on ode3
my_ode := {diff(x(t),t,t)+ 3*t*diff(x(t),t)- x(t)= 0,x(0)= 2, D(x)(0)= 0}:

Try[testnoerror]("test 5", from_polynom_coefficient_ode(t, x, my_ode), 'assign'='ps3');
sol3 := dsolve(my_ode, x(t), 'type=series', 'order' = ord):
sol_pol := convert( sol3, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 6", Truncate(ps3, ord-1), sol_pol);

#Test on ode4
my_ode := {diff(x(t),t,t)- t*diff(x(t),t)- x(t)= 0}:

Try[testnoerror]("test 7", from_polynom_coefficient_ode(t, x, my_ode), 'assign'='ps4');
sol4 := dsolve(my_ode, x(t), 'type=series', 'order' = ord):
sol_pol := convert( sol4, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 8", Truncate(ps4, ord-1), sol_pol);


#end test