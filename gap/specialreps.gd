# GAP Declarations
# $Id: specialreps.gd,v 1.4 2012/02/27 12:26:34 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareOperation( "IndecProjectiveModules", 
                        [IsQuotientOfPathAlgebra, IsList] ); 
DeclareOperation( "IndecInjectiveModules", 
                        [IsAlgebra, IsList] ); 
DeclareOperation( "SimpleModules", 
                        [IsAlgebra] ); 
DeclareOperation( "ZeroModule", 
                        [IsAlgebra] ); 


