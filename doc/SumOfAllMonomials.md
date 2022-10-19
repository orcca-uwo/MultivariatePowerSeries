##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/SumOfAllMonomials",versionnew="{XENOMORPH}") MultivariatePowerSeries[SumOfAllMonomials]
##TITLE(halfline="get the power series representing the sum of all monomials over a set of variables")
##    MultivariatePowerSeries[SumOfAllMonomials]
##ALIAS SumOfAllMonomials, MultivariatePowerSeries:-SumOfAllMonomials, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- SumOfAllMonomials('x')
##
##PARAMETERS
##- 'x' : list of variables
##
##DESCRIPTION
##- The command ~SumOfAllMonomials(x)~ returns the power series containing each monomial in the variables 'x' with coefficient 1; that is,
##  its homogeneous part of degree 'd' is
##  the sum of all monomials of degree 'd'.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- We create a power series that is the sum of all monomials in _x_, and compute its homogeneous
##  part of degree 2 and its truncation at precision 4.
##> a := SumOfAllMonomials([x]);
##<(verification="type") object
##> HomogeneousPart(a, 2);
##< x^2
##> Truncate(a, 4);
##< 1 + x + x^2 + x^3 + x^4
##- We create a power series that is the sum of all monomials in _x_, _y_, and/or _z_. We compute
##  its homogeneous part of degree 2 and its truncation at precision 3.
##> b := SumOfAllMonomials([x,y,z]);
##<(verification="type") object
##> HomogeneousPart(b, 2);
##< x^2 + x*y + x*z + y^2 + y*z + z^2
##> Truncate(b, 3);
##< add(add(add(x^i * y^j * z^k, k=0..3-i-j), j = 0 .. 3-i), i=0..3)
##
##SEEALSO
##- "HomogeneousPart"
##- "GeometricSeries"
##- "Truncate"
##
##XREFMAP
##- HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##- GeometricSeries : Help:MultivariatePowerSeries[GeometricSeries]
##- Truncate : Help:MultivariatePowerSeries[Truncate]
