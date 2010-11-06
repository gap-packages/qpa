# GAP Declarations
# $Id: specialreps.gd,v 1.2 2010/11/06 14:48:40 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareOperation( "IndecomposableProjectiveRepresentations", 
                        [IsSubalgebraFpPathAlgebra] ); 
DeclareOperation( "IndecomposableInjectiveRepresentations", 
                        [IsAlgebra] ); 
DeclareOperation( "VertexSimpleRepresentations", 
                        [IsAlgebra] ); 
DeclareOperation( "ZeroRepresentation", 
                        [IsAlgebra] ); 
