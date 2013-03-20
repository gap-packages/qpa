# GAP Declarations
# This file was generated from
# $Id: algebra.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
KeyDependentOperation("PowerSubalgebra", IsAlgebra, IsPosInt, IsPosInt);
DeclareCategory( "IsOppositeAlgebraElement", IsRingElement );
DeclareCategoryCollections( "IsOppositeAlgebraElement" );
DeclareCategoryFamily( "IsOppositeAlgebraElement" );
DeclareProperty("IsOppositeAlgebra", IsAlgebra);
DeclareAttribute("UnderlyingAlgebra", IsOppositeAlgebra);
DeclareOperation("OppositeAlgebra", [IsAlgebra]);
DeclareAttribute( "RadicalSeriesOfAlgebra", IsAlgebra );