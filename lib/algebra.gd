# GAP Declarations
# This file was generated from
# $Id: algebra.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
KeyDependentOperation( "PowerSubalgebra", IsAlgebra, IsPosInt, IsPosInt);
DeclareCategory( "IsOppositeAlgebraElement", IsRingElement );
DeclareCategoryCollections( "IsOppositeAlgebraElement" );
DeclareCategoryFamily( "IsOppositeAlgebraElement" );
DeclareProperty( "IsOppositeAlgebra", IsAlgebra);
DeclareAttribute( "UnderlyingAlgebra", IsOppositeAlgebra);
DeclareOperation( "OppositeAlgebra", [IsAlgebra]);
DeclareAttribute( "RadicalSeriesOfAlgebra", IsAlgebra );
DeclareProperty( "IsSemisimpleAlgebra", IsAlgebra );
DeclareProperty( "IsHereditaryAlgebra", IsAlgebra );
DeclareProperty( "IsFiniteGlobalDimensionAlgebra", IsAlgebra );
DeclareProperty( "IsGorensteinAlgebra", IsAlgebra );
DeclareProperty( "IsRadicalSquareZeroAlgebra", IsAlgebra );
InstallTrueMethod( IsAlgebra, IsSemisimpleAlgebra );
InstallTrueMethod( IsAlgebra, IsHereditaryAlgebra );
InstallTrueMethod( IsAlgebra, IsFiniteGlobalDimensionAlgebra );
InstallTrueMethod( IsAlgebra, IsGorensteinAlgebra );
InstallTrueMethod( IsAlgebra, IsRadicalSquareZeroAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsSemisimpleAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsHereditaryAlgebra );
DeclareAttribute( "GlobalDimension", IsAlgebra );