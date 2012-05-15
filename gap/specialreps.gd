# GAP Declarations
# $Id: specialreps.gd,v 1.5 2012/05/15 06:58:20 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareAttribute( "IndecProjectiveModules", IsQuotientOfPathAlgebra ); 
DeclareAttribute( "IndecInjectiveModules", IsAlgebra ); 
DeclareOperation( "SimpleModules", [IsAlgebra] ); 
DeclareOperation( "ZeroModule", [IsAlgebra] ); 