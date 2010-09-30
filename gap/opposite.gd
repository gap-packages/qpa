# GAP Declarations
# $Id: opposite.gd,v 1.1 2010/09/30 14:02:22 oysteini Exp $

# opposite.gd: Operations for producing the opposites of quivers and
# path algebras.

DeclareAttribute( "OppositeQuiver", IsQuiver );
DeclareOperation( "OppositePath", [ IsPath ] );
DeclareAttribute( "OppositeQuiverNameMap", IsQuiver );

DeclareAttribute( "OppositePathAlgebra", IsAlgebra );
DeclareGlobalFunction( "OppositePathAlgebraElement"); # should be operation with args [ IsElementOfPathAlgebra ]
DeclareOperation( "OppositeRelations", [ IsDenseList ] );
