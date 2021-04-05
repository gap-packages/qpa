# GAP Declarations
# $Id: patensor.gd,v 1.3 2012/08/01 16:01:10 sunnyquiver Exp $

# patensor.gd: Tensor product of path algebras.

# Declares the operation TensorProductOfAlgebras, which takes two
# algebras (over the same field) as arguments and returns their tensor
# product (a new algebra).

# The operation SimpleTensor takes as arguments a list of two algebra
# elements -- one from each of the operands of a tensor product -- as
# well as the tensor product algebra, and returns the corresponding
# element in the tensor product.

# The QuiverProduct operation is used for creating the quiver for the
# tensor product of two path algebras.

# Example (after reading patensor.gd and patensor.gi):
#
# gap> q := Quiver( 3, [ [ 1, 2, "a" ], [ 2, 3, "b" ] ] );
# <quiver with 3 vertices and 2 arrows>
# gap> pa := PathAlgebra( Rationals, q );
# <algebra-with-one over Rationals, with 5 generators>
# gap> tp := TensorProductOfAlgebras( pa, pa );
# <algebra-with-one over Rationals, with 21 generators>
# gap> Print( QuiverOfPathAlgebra( tp ), "\n" );
# Quiver( ["v1_v1","v1_v2","v1_v3","v2_v1","v2_v2","v2_v3","v3_v1","v3_v2","v3_v3"],
#         [["v1_v1","v1_v2","v1_a"], ["v1_v2","v1_v3","v1_b"],
#          ["v2_v1","v2_v2","v2_a"], ["v2_v2","v2_v3","v2_b"],
#          ["v3_v1","v3_v2","v3_a"], ["v3_v2","v3_v3","v3_b"],
#          ["v1_v1","v2_v1","a_v1"], ["v1_v2","v2_v2","a_v2"],
#          ["v1_v3","v2_v3","a_v3"], ["v2_v1","v3_v1","b_v1"],
#          ["v2_v2","v3_v2","b_v2"], ["v2_v3","v3_v3","b_v3"]] )
# gap> RelatorsOfFpAlgebra( tp );
# [ (-1)*v1_a*a_v2+(1)*a_v1*v2_a, (-1)*v1_b*a_v3+(1)*a_v2*v2_b,
#   (-1)*v2_a*b_v2+(1)*b_v1*v3_a, (-1)*v2_b*b_v3+(1)*b_v2*v3_b ]
# gap> SimpleTensor( [ pa.a * pa.b, pa.a ], tp );
# [(1)*a_v1*b_v1*v3_a]


DeclareOperation( "QuiverProduct", [ IsQuiver, IsQuiver ] );
DeclareAttribute( "QuiverProductDecomposition", IsQuiver);
DeclareOperation( "IncludeInProductQuiver", [ IsDenseList, IsQuiver and HasQuiverProductDecomposition ] );
DeclareOperation( "ProjectFromProductQuiver", [ IsPosInt, IsPath ] );

DeclareCategory( "IsQuiverProductDecomposition", IsList );

DeclareOperation( "WalkOfPathOrVertex", [ IsPath ] );
DeclareOperation( "ReversePath", [ IsPath ] );

DeclareOperation( "IncludeInPathAlgebra", [ IsPath, IsQuiverAlgebra ] );
DeclareOperation( "VerticesOfPathAlgebra", [ IsQuiverAlgebra ] );
DeclareGlobalFunction( "PathAlgebraElementTerms" );
DeclareOperation( "SimpleTensor", [ IsDenseList, IsQuiverAlgebra ] );

DeclareOperation( "TensorProductOfAlgebras", [ IsQuiverAlgebra, IsQuiverAlgebra ] );
DeclareGlobalFunction( "TensorProductOfPathAlgebras" );
DeclareAttribute( "TensorProductDecomposition", IsAlgebra );
DeclareOperation( "TensorAlgebraInclusion", [  IsQuiverAlgebra, IS_INT ] );

DeclareAttribute( "EnvelopingAlgebra", IsQuiverAlgebra );
DeclareProperty( "IsEnvelopingAlgebra", IsAlgebra );

DeclareAttribute( "AlgebraAsModuleOverEnvelopingAlgebra", IsQuiverAlgebra );
DeclareSynonym( "AlgebraAsModuleOfEnvelopingAlgebra", AlgebraAsModuleOverEnvelopingAlgebra );
DeclareAttribute( "DualOfAlgebraAsModuleOverEnvelopingAlgebra", IsQuiverAlgebra );
DeclareAttribute( "TrivialExtensionOfQuiverAlgebraLevel", IsQuiverAlgebra );
DeclareAttribute( "TrivialExtensionOfQuiverAlgebra", IsQuiverAlgebra );
DeclareAttribute( "EnvelopingAlgebraHomomorphism", IsAlgebraHomomorphism );
DeclareAttribute( "TrivialExtensionOfQuiverAlgebraProjection", IsQuiverAlgebra );