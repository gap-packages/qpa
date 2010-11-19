# GAP Declarations
# $Id: specialreps.gd,v 1.3 2010/11/19 13:24:48 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareOperation( "IndecomposableProjectiveRepresentations", 
                        [IsSubalgebraFpPathAlgebra, IsList] ); 
DeclareOperation( "IndecomposableInjectiveRepresentations", 
                        [IsAlgebra, IsList] ); 
DeclareOperation( "VertexSimpleRepresentations", 
                        [IsAlgebra] ); 
DeclareOperation( "ZeroRepresentation", 
                        [IsAlgebra] ); 


