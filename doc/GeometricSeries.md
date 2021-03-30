##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/GeometricSeries",versionnew="{XENOMORPH}") MultivariatePowerSeries[GeometricSeries]
##TITLE(halfline="Create the geometric series over a set of variables")
##    MultivariatePowerSeries[GeometricSeries]
##ALIAS GeometricSeries, MultivariatePowerSeries:-GeometricSeries, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- GeometricSeries('x')
##
##PARAMETERS
##- 'x' : variable, or nonempty list of variables
##
##DESCRIPTION
##- ~GeometricSeries(x)~ creates the geometric power series in 'x'.
##- If 'x' has more than one element then the result is the geometric series
##  in the sum of the variables in 'x'.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- We create the geometric power series in _x_, that is, the power series for _1/(1+x)_.
##> a := GeometricSeries(x);
##<(verification="type") object
##- If we truncate _a_ at homogeneous degree 5, this is what we obtain.
##> Truncate(a, 5);
##< 1 + x + x^2 + x^3 + x^4 + x^5
##- We create the geometric power series in _x + y_, that is, the power series for _1/(1+x+y)_.
##> b := GeometricSeries([x,y]);
##- The homogeneous part of _b_ of degree 4 is the following expression.
##> HomogeneousPart(b, 4);
##< expand((x + y)^4)
##- If we truncate _b_ at homogeneous degree 3, this is what we obtain.
##> Truncate(b, 3);
##< add(expand((x + y)^i), i=0..3)
##
##SEEALSO
##- "PowerSeries"
##- "SumOfAllMonomials"
##- "HomogeneousPart"
##- "Truncate"
##
##XREFMAP
##- PowerSeries : Help:MultivariatePowerSeries[PowerSeries]
##- SumOfAllMonomials : Help:MultivariatePowerSeries[SumOfAllMonomials]
##- HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##- Truncate : Help:MultivariatePowerSeries[Truncate]
