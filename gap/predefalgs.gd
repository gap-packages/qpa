# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gd,v 1.3 2012/06/10 08:19:10 sunnyquiver Exp $
DeclareOperation( "NakayamaAlgebra", [IsList, IsField] );
DeclareOperation( "CanonicalAlgebra", [IsField, IsList, IsList] );
DeclareOperation( "KroneckerAlgebra", [IsField, IS_INT] );
DeclareCategory("IsNakayamaAlgebras", IsAlgebra );
DeclareCategory("IsCanonicalAlgebras", IsAlgebra );
DeclareCategory("IsKroneckerAlgebras", IsAlgebra );