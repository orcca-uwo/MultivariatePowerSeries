#test
with(TestTools):

Try[testnoerror]("test 0", with(MultivariatePowerSeries));

kernelopts(opaquemodules=false):
PSO := MultivariatePowerSeries:-PowerSeriesObject:
from_linear_coefficient_ode := MultivariatePowerSeries:-PowerSeriesObject:-from_linear_coefficient_ode:
kernelopts(opaquemodules=true):
# TODO: Change the orders
ord := 10;
# the number of terms we consider for all tests

# tests for ode1
my_ode := {diff(x(t),t)- x(t)= 0}:

Try[testnoerror]("test 1", from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps1');
sol1 := dsolve(my_ode, x(t), 'type=series', 'order' = ord);
sol_pol := convert( sol1, 'polynom');
sol_pol := rhs(sol_pol):;
Try[verify,normal]("test 2", Truncate(ps1, ord-1), sol_pol);

# tests for ode2
my_ode:= {diff(x(t), t, t)+x(t)=0}:

Try[testnoerror]("test 3", from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps2');
sol2 := dsolve(my_ode, x(t), 'type=series', 'order'= ord):
sol_pol := convert( sol2, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 4", Truncate(ps2, ord-1), sol_pol);

# tests for ode 3
my_ode :={diff(x(t), t, t)+ 6*diff(x(t),t)= -5*x(t)}:

Try[testnoerror]("test 5", from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps3');
sol3 := dsolve(my_ode, x(t), 'type=series','order'= ord):
sol_pol := convert( sol3, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 6", Truncate(ps3, ord-1), sol_pol);

# tests for ode 4
my_ode := {diff(x(t),t,t)+ 6*diff(x(t),t)+9*x(t)=0}:

Try[testnoerror]("test 7",from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps4');
sol4 := dsolve(my_ode, x(t), 'type=series','order'= ord):
sol_pol := convert( sol4, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 8", Truncate(ps4, ord-1), sol_pol);

# tests for ode 5
my_ode := {diff(x(t),t,t,t)-3*diff(x(t),t,t)-diff(x(t),t)+3*x(t)=0}:

Try[testnoerror]("test 9",from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps5');
sol5 := dsolve(my_ode, x(t), 'type=series','order'= ord):
sol_pol := convert( sol5, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 10", Truncate(ps5, ord-1), sol_pol);

# tests for ode 6
my_ode:= {diff(x(t),t,t,t)-diff(x(t),t,t)-20*x(t)=0}:

Try[testnoerror]("test 11",from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps6');
sol6 := dsolve(my_ode, x(t), 'type=series','order'= ord):
sol_pol := convert( sol6, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 12", Truncate(ps6, ord-1), sol_pol);

# tests for ode 7
my_ode:= {diff(x(t),t,t,t,t)-4*diff(x(t),t,t,t)+8*diff(x(t),t,t)-8*diff(x(t),t)+4*x(t)= 0}:

Try[testnoerror]("test 13",from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps7');
sol7 := dsolve(my_ode, x(t), 'type=series','order'= ord):
sol_pol := convert( sol7, 'polynom'):
sol_pol := rhs(sol_pol):
Try[verify,normal]("test 14", Truncate(ps7, ord-1), sol_pol);

# Test ode without enough initial conditions.
my_ode := {diff(x(t),t)- x(t)= 0};
Try[testnoerror]("test 15", from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps');

# We compute the solution and write it in a nice format.
sol := dsolve(my_ode, x(t), 'type=series', 'order' = ord );
sol_pol := convert( sol, 'polynom'); 
sol_pol := rhs(sol_pol);
Try[verify,normal]("test 17", Truncate(ps, ord-1), sol_pol);

# "Model" test.
my_ode := {diff(x(t),t,t,t)-diff(x(t),t,t)-20*diff(x(t),t)=0};
Try[testnoerror]("test 18", from_linear_coefficient_ode(t, x, my_ode), 'assign'='ps');

# We compute the solution and write it in a nice format.
sol := dsolve(my_ode, x(t), 'type=series', 'order' = ord );
sol_pol := convert( sol, 'polynom'); 
sol_pol := rhs(sol_pol);

Try[verify,normal]("test 19", Truncate(ps, ord-1), sol_pol);

#end test