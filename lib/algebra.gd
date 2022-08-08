# GAP Declarations
# This file was generated from
# $Id: algebra.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
KeyDependentOperation( "PowerSubalgebra", IsAlgebra, IsPosInt, IsPosInt);
DeclareCategory( "IsOppositeAlgebraElement", IsRingElement );
DeclareCategoryCollections( "IsOppositeAlgebraElement" );
DeclareCategoryFamily( "IsOppositeAlgebraElement" );
DeclareProperty( "IsOppositeAlgebra", IsAlgebra);
InstallTrueMethod( IsAlgebra, IsOppositeAlgebra );
DeclareAttribute( "UnderlyingAlgebra", IsOppositeAlgebra);
DeclareOperation( "OppositeAlgebra", [IsAlgebra]);
DeclareAttribute( "RadicalSeriesOfAlgebra", IsAlgebra );
DeclareProperty( "IsSemisimpleAlgebra", IsAlgebra );
DeclareProperty( "IsHereditaryAlgebra", IsAlgebra );
DeclareProperty( "IsFiniteGlobalDimensionAlgebra", IsAlgebra );
DeclareProperty( "IsGorensteinAlgebra", IsAlgebra );
DeclareProperty( "IsRadicalSquareZeroAlgebra", IsAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsSemisimpleAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsHereditaryAlgebra );
DeclareAttribute( "GlobalDimension", IsAlgebra );
DeclareProperty( "IsBasicAlgebra", IsAlgebra );
DeclareProperty( "IsElementaryAlgebra", IsAlgebra );
InstallTrueMethod( IsBasicAlgebra, IsElementaryAlgebra );
DeclareOperation( "LiftingIdempotent", [IsAlgebraGeneralMapping, IsObject ] ); 
DeclareOperation( "LiftingCompleteSetOfOrthogonalIdempotents", [IsAlgebraGeneralMapping, IsHomogeneousList ] ); 
DeclareOperation( "FindMultiplicativeIdentity", [ IsAlgebra ] );
DeclareOperation( "AlgebraAsQuiverAlgebra", [ IsAlgebra ]);
DeclareOperation( "BlocksOfAlgebra", [ IsAlgebra ]);
DeclareAttribute( "OppositeAlgebraHomomorphism", IsAlgebraHomomorphism );

