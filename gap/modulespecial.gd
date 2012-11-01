# GAP Declarations
# $Id: specialreps.gd,v 1.6 2012/08/01 16:01:10 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareAttribute( "IndecProjectiveModules", IsQuotientOfPathAlgebra ); 
DeclareAttribute( "IndecInjectiveModules", IsAlgebra ); 
DeclareAttribute( "SimpleModules", IsAlgebra ); 
DeclareAttribute( "ZeroModule", IsAlgebra );
DeclareAttribute( "BasisOfProjectives", IsQuotientOfPathAlgebra );  
