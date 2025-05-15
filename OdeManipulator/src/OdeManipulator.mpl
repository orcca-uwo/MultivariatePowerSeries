unprotect( 'OdeManipulator' ):

OdeManipulator := module()
    option package;

    export PowerSeriesWrapper; # A wrapper for the pso.
    export OdeObject; # OdeObject.


$include "MultivariatePowerSeries/OdeManipulator/OdeObject/src/OdeObject.mm"
$include "MultivariatePowerSeries/OdeManipulator/OdeObject/src/PowerSeriesWrapper.mm"

end module:

protect( 'OdeManipulator' ):
#with(LibraryTools);
#savelib( 'MultivariatePowerSeries' ):
 
