##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/GrevlexPositive",versionnew="{ZONKEY}") MultivariatePowerSeries[PuiseuxSeries]
##TITLE(halfline="to know if a vector or a list is grevlex positive")
##    MultivariatePowerSeries[GrevlexPositive]
##ALIAS GrevlexPositive, MultivariatePowerSeries:-GrevlexPositive, MultivariatePowerSeries
##AUTHOR Matt Calder, Juan Gonzalez Trochez jgonza55@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca, Erik Postma
##
##CALLINGSEQUENCE
##- GrevlexPositive('v')
##
##PARAMETERS
##- 'v' : a rational list or vector
##
##DESCRIPTION
##- We get true if 'v' is grevlex greater than '0' and false otherwise.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##-(lead=indent) define two lists and a vector. 
##> l1 := [1,2/5,4];
##> l2 := [-1,1];
##> v := <1,2,-3>;
##- We check if 'l1' is grevlex positive.
##> GrevlexPositive(l1);
##- We check if 'l2' is grevlex positive.
##> GrevlexPositive(l2);
##- We check if 'v' is grevlex positive.
##> GrevlexPositive(v);
##
##SEEALSO
##- "MultivariatePowerSeries"
##- "PuiseuxSeries"
## 
##XREFMAP
##- MultivariatePowerSeries : Help:MultivariatePowerSeries
##- PuiseuxSeries : Help:MultivariatePowerSeries[PuiseuxSeries]