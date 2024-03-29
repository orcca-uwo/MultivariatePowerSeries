##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/GetMonomial",versionnew="{ZONKEY}") MultivariatePowerSeries[GetMonomial]
##TITLE(halfline="get the monomial that multiplies a Puiseux series")
##    MultivariatePowerSeries[GetMonomial]
##ALIAS GetMonomial, MultivariatePowerSeries:-GetMonomial, MultivariatePowerSeries
##AUTHOR Matt Calder, Juan Gonzalez Trochez jgonza55@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca, Erik Postma
##
##CALLINGSEQUENCE
##- GetMonomial('s')
##
##PARAMETERS
##- 's' : Puiseux series generated by this package
##
##DESCRIPTION
##- This command returns the monomial _X^e_ that multiplies a Puiseux series.
##INCLUDE PuiseuxSeriesIntro.mi
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##-(lead=indent) Create a Puiseux series. 
##> p := PowerSeries(1+u*v); X := [x, y]; U := [u, v]; R := [[1,0], [1,-1]]; E := [x=-5, y=3];
##> s := PuiseuxSeries(p, X, U, R, E);
##<(verification="type") object
##- We get the monomial that multiplies 's'.
##> GetMonomial(s);
##<	y^3/(x^5)
##
##SEEALSO
##- "MultivariatePowerSeries"
##- "PuiseuxSeries"
## 
##INCLUDE PxRef.mi
##
##XREFMAP
##- MultivariatePowerSeries : Help:MultivariatePowerSeries
##- PuiseuxSeries : Help:MultivariatePowerSeries[PuiseuxSeries]
