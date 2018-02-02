# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)

DeclareProperty( "IsIdealInPathAlgebra", IsFLMLOR );
DeclareProperty( "IsAdmissibleIdeal",  IsIdealInPathAlgebra );
DeclareProperty( "IsMonomialIdeal", IsIdealInPathAlgebra );
DeclareOperation( "ProductOfIdeals", [ IsIdealInPathAlgebra, IsIdealInPathAlgebra ]);
DeclareProperty( "IsAdmissibleQuotientOfPathAlgebra", IsQuotientOfPathAlgebra );
InstallTrueMethod( IsQuotientOfPathAlgebra, IsAdmissibleQuotientOfPathAlgebra );
DeclareProperty( "IsGentleAlgebra",  IsQuiverAlgebra );
InstallTrueMethod( IsGorensteinAlgebra, IsGentleAlgebra ); 
InstallTrueMethod( IsSpecialBiserialAlgebra, IsGentleAlgebra ); 
InstallTrueMethod( IsFiniteGlobalDimensionAlgebra, IsPosetAlgebra );
InstallTrueMethod( IsElementaryAlgebra, IsAdmissibleQuotientOfPathAlgebra );
DeclareOperation( "MinimalGeneratingSetOfIdeal", [ IsAdmissibleIdeal ] );