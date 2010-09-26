# GAP Declarations
# $Id: specialreps.gd,v 1.1 2010/09/26 19:55:42 sunnyquiver Exp $

# specialreps.gd: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

DeclareOperation( "IndecomposableProjectiveRepresentations", 
                        [IsPathAlgebra, IsList] ); 
DeclareOperation( "IndecomposableInjectiveRepresentations", 
                        [IsPathAlgebra, IsList] ); 
DeclareOperation( "VertexSimpleRepresentations", 
                        [IsPathAlgebra, IsList] ); 
DeclareOperation( "ZeroRepresentation", 
                        [IsPathAlgebra, IsList] ); 