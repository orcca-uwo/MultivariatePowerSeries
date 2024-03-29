##ENCODING ISO-8859-1
##PROCEDURE(help,nospec,label="MultivariatePowerSeries/Copy",versionnew="{XENOMORPH}",versionupdated="{ZONKEY}") MultivariatePowerSeries[Copy]
##TITLE(halfline="copy a power series, a Puiseux series, or univariate polynomial over these series")
##    MultivariatePowerSeries[Copy]
##ALIAS Copy, MultivariatePowerSeries:-Copy, MultivariatePowerSeries
##AUTHOR Ali Asadi masadi4@uwo.ca, Alex Brandt abrandt5@uwo.ca, Marc Moreno Maza moreno@csd.uwo.ca
##
##CALLINGSEQUENCE
##- Copy('p')
##- Copy('s')
##- Copy('u')
##
##PARAMETERS
##- 'p' : power series generated by this package
##- 's' : Puiseux series generated by this package
##- 'u' : univariate polynomial over power series or Puiseux series generated by this package
##
##DESCRIPTION
##- ~Copy(p)~, ~Copy(s)~, and ~Copy(u)~ return copies of ~p~, ~s~, and ~u~, respectively. If the original object is subsequently
## modified (for example, by computing extra coefficients, or modifying the "display style"), these changes are not
## reflected in the copy, and vice versa. Note that the original object and its copy may share ancestors, such as
## power series or Puiseux series objects from which they were computed.
##
##INCLUDE assignment_warning.mi
##
##EXAMPLES
##> with(MultivariatePowerSeries):
##- We create two power series, _a_ and _b_.
##> a := Inverse(PowerSeries(1 + x - y)):
##<(verification="type") object
##> b := Inverse(PowerSeries(y^2 - x + 1)):
##<(verification="type") object
##- The power series _c_ keeps a record of _a_ and _b_ as its ancestors.
##> c := a + b;
##<(verification="type") object
##> d := Copy(c);
##<(verification="type") object
##- We can set the "display styles" for _c_ and _d_ independently. We ensure that enough
##  terms are computed to show the difference.
##> SetDisplayStyle(c, ['precision' = 7]);
##> SetDisplayStyle(d, ['precision' = 4]);
##> UpdatePrecision(c, 7): UpdatePrecision(d, 7):
##> c;
##<(verification="type") object
##> d;
##<(verification="type") object
##- The power series _d_ is a copy of _c_, however, they share _a_ and _b_, therefore when more
##  terms of _c_ are computed, the precision of _a_ and _b_ will be updated too. Consequently, the
##  computation of new terms of _d_ must be cheaper as it does not involve computing the
##  coefficients of _a_ and _b_.
##> gc():
##> CodeTools:-Usage(HomogeneousPart(c, 500)):
##> Precision(a);
##< 500
##> Precision(b);
##< 500
##> Precision(c);
##< 500
##> Precision(d);
##< 7
##> gc():
##> CodeTools:-Usage(HomogeneousPart(d, 500)):
##> Precision(d);
##< 500
##
##SEEALSO
##- "PowerSeries"
##- "PuiseuxSeries"
##- "Inverse"
##- "Precision"
##- "HomogeneousPart"
##- "Truncate"
##
##INCLUDE PxRef.mi
##
##XREFMAP
##- PowerSeries : Help:MultivariatePowerSeries[PowerSeries]
##- PuiseuxSeries : Help:MultivariatePowerSeries[PuiseuxSeries]
##- Inverse : Help:MultivariatePowerSeries[Inverse]
##- Precision : Help:MultivariatePowerSeries[Precision]
##- HomogeneousPart : Help:MultivariatePowerSeries[HomogeneousPart]
##- Truncate : Help:MultivariatePowerSeries[Truncate]
##- "display style" : Help:MultivariatePowerSeries[SetDisplayStyle]
##- "display styles" : Help:MultivariatePowerSeries[SetDisplayStyle]
