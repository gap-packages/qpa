# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gd,v 1.4 2012/06/20 17:00:31 andrzejmroz Exp $
DeclareOperation( "NakayamaAlgebra", [IsList, IsField] );
DeclareOperation( "CanonicalAlgebra", [IsField, IsList, IsList] );
DeclareOperation( "KroneckerAlgebra", [IsField, IS_INT] );
DeclareCategory("IsNakayamaAlgebras", IsAlgebra );
DeclareCategory("IsCanonicalAlgebras", IsAlgebra );
DeclareCategory("IsKroneckerAlgebras", IsAlgebra );

DeclareProperty("IsSpecialBiserialQuiver", IsQuiver);
DeclareProperty("IsSpecialBiserialAlgebra", IsPathAlgebra);
DeclareProperty("IsSpecialBiserialAlgebra", IsQuotientOfPathAlgebra);
DeclareProperty("IsStringAlgebra", IsPathAlgebra);
DeclareProperty("IsStringAlgebra", IsQuotientOfPathAlgebra);
#DeclareProperty("IsGentleAlgebra", IsPathAlgebra);
#DeclareProperty("IsGentleAlgebra", IsQuotientOfPathAlgebra);