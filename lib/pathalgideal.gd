# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)

DeclareProperty( "IsIdealInPathAlgebra", IsFLMLOR );
DeclareProperty( "IsAdmissibleIdeal",  IsIdealInPathAlgebra );
DeclareProperty( "IsMonomialIdeal", IsIdealInPathAlgebra );
DeclareOperation( "ProductOfIdeals", [ IsIdealInPathAlgebra, IsIdealInPathAlgebra ]);
