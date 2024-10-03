unprotect( 'OdeManipulator' ):

OdeManipulator := module()
    option package;

    export OdeObject; # OdeObject


$include "MultivariatePowerSeries/OdeManipulator/OdeObject/src/OdeObject.mm"

end module:

protect( 'OdeManipulator' ):
#savelib( 'MultivariatePowerSeries' ):
 
