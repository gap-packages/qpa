# GAP Implementation
# This contains various tools for determining a representation type of a (quotient of) path algebra
# (created A. Mroz, 08.01.2013)

DeclareAttribute( "IsFiniteTypeAlgebra", IsPathAlgebra );
DeclareAttribute( "IsFiniteTypeAlgebra", IsQuotientOfPathAlgebra );
DeclareOperation(  "BongartzTest",  [ IsQuiverAlgebra, IS_INT ] );
