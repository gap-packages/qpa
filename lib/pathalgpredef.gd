# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gd,v 1.5 2012/08/01 16:01:10 sunnyquiver Exp $
DeclareOperation( "NakayamaAlgebra", [IsField, IsList] );
DeclareOperation( "CanonicalAlgebra", [IsField, IsList, IsList] );
DeclareOperation( "KroneckerAlgebra", [IsField, IS_INT] );
DeclareProperty("IsNakayamaAlgebra", IsAlgebra );
DeclareProperty("IsCanonicalAlgebra", IsAlgebra );
DeclareProperty("IsKroneckerAlgebra", IsAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsCanonicalAlgebra );
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsKroneckerAlgebra );
InstallTrueMethod( IsHereditaryAlgebra, IsKroneckerAlgebra );
DeclareProperty("IsSpecialBiserialQuiver", IsQuiver);
DeclareProperty("IsSpecialBiserialAlgebra", IsPathAlgebra);
DeclareProperty("IsSpecialBiserialAlgebra", IsQuotientOfPathAlgebra);
DeclareProperty("IsStringAlgebra", IsPathAlgebra);
DeclareProperty("IsStringAlgebra", IsQuotientOfPathAlgebra);
DeclareProperty( "IsPosetAlgebra", IsQuiverAlgebra ); 
DeclareOperation( "PosetAlgebra", [ IsField, IsPoset ] ); 
DeclareAttribute( "PosetOfPosetAlgebra", IsPosetAlgebra );
DeclareOperation( "BrauerConfigurationAlgebra", [IsField, IsList] );
DeclareOperation( "PreprojectiveAlgebra", [ IsPathAlgebraMatModule, IsInt ] );
DeclareOperation( "AdmissibleSequenceGenerator", [ IsPosInt, IsPosInt ] );