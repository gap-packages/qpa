# GAP Implementation
# This file was generated from 
# $Id: functors.gi,v 1.3 2010/11/06 14:48:40 sunnyquiver Exp $
InstallMethod ( DualOfPathAlgebraMatModule,
    "for a representation of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )
#
#   M     = a representation of the quiver Q over K with rels
#
        local A_op, KQ, KQ_op, dim_vect_M, num_vert, vertices, i, j, 
  	      count, mat_M, mat_M_op, X, mats, arrows, a, 
	      vert_in_quiver, origin, target, N;

    A_op  := OppositePathAlgebra(RightActingAlgebra(M));
    KQ_op := OriginalPathAlgebra(A_op);
    dim_vect_M := DimensionVector(M); 
    num_vert := Length(DimensionVector(M));
    mat_M    := MatricesOfPathAlgebraMatModule(M); 
    arrows   := ArrowsOfQuiver(QuiverOfPathAlgebra(KQ_op));
    mat_M_op := List(mat_M, X -> TransposedMat(X));
    vert_in_quiver := VerticesOfQuiver(QuiverOfPathAlgebra(KQ_op));
    mats := [];
    for a in arrows do
    	i := Position(arrows,a); 
	origin := Position(vert_in_quiver,SourceOfPath(a));
	target := Position(vert_in_quiver,TargetOfPath(a)); 
 	if ( dim_vect_M[origin] <> 0 ) and ( dim_vect_M[target] <> 0 ) then 
	   Add(mats,[a,mat_M_op[i]]);
	else
	   Add(mats,[a,[dim_vect_M[origin],dim_vect_M[target]]]);
        fi;
    od;

    N:= RightModuleOverPathAlgebra(KQ_op,mats);
    SetDualOfPathAlgebraMatModule(N,M);
    return N;
#
#  Want to have:
#    return RightModuleOverPathAlgebra(A_op,mats);
#
end
);

