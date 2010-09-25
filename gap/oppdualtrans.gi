InstallMethod ( OppositeOfQuiver,
# Written by Oeyvind Solberg, 08.06.2010
# Dependent on: VerticesOfQuiver, ArrowsOfQuiver, SourceOfPath, 
#               TargetOfPath, Quiver
#
    "for a quiver",
    [IsQuiver], 0,
    function( Q )
        local arrows, vertices, names_of_opposite_vertices, 
        opposite_vertex, names_of_opposite_arrows, opposite_arrow, 
        list_of_opposite_arrows, source, target, new_entry, i;
        
        if not (  IsQuiver(Q) )  then
            Error("Need to enter a quiver as an argument!");
        fi;
#
# Giving names to the vertices in the opposite quiver. 
# Name in opposite quiver:  "name in original quiver"+"_op"
#
        names_of_opposite_vertices := [];
        vertices := VerticesOfQuiver(Q);
        for i in [1..Length(vertices)] do
            opposite_vertex := Concatenation(String(vertices[i]),"_op");
            Add(names_of_opposite_vertices, opposite_vertex);
        od;
#
# Giving names to the arrows in the opposite quiver. 
# Name in opposite quiver:  "name in original quiver"+"_op"
#
        names_of_opposite_arrows := [];
        list_of_opposite_arrows := [];
        arrows:= ArrowsOfQuiver(Q);
        for i in [1..Length(arrows)] do
            opposite_arrow := Concatenation(String(arrows[i]),"_op");
            source := 
                Concatenation(String(TargetOfPath(arrows[i])),"_op");
            target := 
                Concatenation(String(SourceOfPath(arrows[i])),"_op");
            new_entry := [];
            Add(new_entry,String(source));
            Add(new_entry,String(target));
            Add(new_entry,String(opposite_arrow));
            Add(list_of_opposite_arrows,new_entry);
            Add(names_of_opposite_arrows, opposite_arrow);
        od;

        return Quiver(names_of_opposite_vertices,list_of_opposite_arrows);
    end
);

InstallMethod ( OppositeOfRelations,
# Written by Oeyvind Solberg, 09.06.2010
# Dependent on: ArrowsOfQuiver, QuiverOfPathAlgebra, Tip, TipMonomial, 
#               TipCoefficient, WalkOfPath. 
#
    "for a path algebra and a set of relations",
    [IsPathAlgebra, IsPathAlgebra, IsList], 0,
    function( kQ_op, kQ, rels )
#
#   kQ_op  = path algebra of the opposite quiver of Q
#   kQ   = path algebra of Q 
#   rels = relations in kQ
#
        local arrows, opposite_arrows, i, j, r, rTip, rMon, rCoeff, 
                rWalk, op_rel, opposite_relations;
        
        opposite_arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(kQ_op)); 
        arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(kQ));
        opposite_relations := [];
        for i in [1..Length(rels)] do
            r := rels[i];
            op_rel := Zero(kQ_op); 
            while r <> Zero(kQ) do
                rTip := Tip(r);
                rMon := TipMonomial(rTip);
                rCoeff := TipCoefficient(rTip); 
                rWalk := WalkOfPath(rMon); 
                rWalk := Reversed(rWalk); 
                for j in [1..Length(rWalk)] do
                    rWalk[j] := opposite_arrows[Position(arrows,rWalk[j])]; 
                od;
                op_rel := op_rel + 
                          One(kQ_op)*rCoeff*Product(rWalk{[1..Length(rWalk)]});
                r := r - rTip;
            od;
            Add(opposite_relations,op_rel);
        od;

        return opposite_relations;
    end
);

InstallMethod ( OppositeOfQuotientOfPathAlgebra,
# Written by Oeyvind Solberg, 09.06.2010
# Dependent on: OppositeOffQuiver, QuiverOfPathAlgebra, LeftActingDomain, 
#                PathAlgebra, OppositeOfRelations, Ideal. 
#
    "for a path algebra and a set of relations",
    [ IsPathAlgebra, IsList ], 0,
    function( KQ, rels )
#
#   KQ   = path algebra of Q 
#   rels = relations in kQ
#
        local Q_op, K, KQ_op, opposite_relations, I_op;
        
        Q_op := OppositeOfQuiver(QuiverOfPathAlgebra(KQ));
        K    := LeftActingDomain(KQ);
        KQ_op:= PathAlgebra(K,Q_op);
        opposite_relations := OppositeOfRelations(KQ_op,KQ,rels);
        I_op := Ideal(KQ_op,opposite_relations);

        return KQ_op/I_op;

    end
);

InstallMethod ( MatricesOfPathAlgebraMatModule,
# Written by Oeyvind Solberg, 09.06.2010
# Dependent on: OppositeOffQuiver, QuiverOfPathAlgebra, LeftActingDomain, 
#                PathAlgebra, OppositeOfRelations, Ideal. 
#
    "for a representation of a quiver",
    [ IsPathAlgebraMatModule ], 0,
    function( M )
#
#   M     = a representation of the quiver Q over K
#
        local basis_M, A, Q, num_vert, vertices, arrows_as_paths, 
                dim_vect_M, i, j, count, dim_size, s, t, v, dom_a, 
                im_a, mat_op, mat, pd, pi, m, a;
        
        basis_M  := Basis(M);
        A        := RightActingAlgebra(M);
        Q        := QuiverOfPathRing(A);
        num_vert := Length(VerticesOfQuiver(Q));
        vertices:=[];
        for i in [1..num_vert]  do
            Add(vertices,GeneratorsOfAlgebra(A)[i]);
        od;
        arrows_as_paths:=[];
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            Add(arrows_as_paths,GeneratorsOfAlgebra(A)[i+num_vert]);
        od; 
#
# Finding the dimension vector for M
#
        dim_vect_M := []; 
        for i in [1.. num_vert] do
            count := 0;
            for j in [1..Length(basis_M)] do
                if basis_M[j]^vertices[i] <> Zero(M) then
                    count := count + 1;   
                fi;
            od;
            Add(dim_vect_M,count);    
        od;
#
# Finding intervals for the basis vectors for the different vertices
#
        dim_size := [];
        s := 1;
        t := 0;
        for i in [1..Length(vertices)] do
            t := dim_vect_M[i] + t;
            Add(dim_size,[s..t]);
            s := t + 1;
        od;  
#
# Finding the matrices of the linear transformations for the different
# arrows in the quiver Q. 
#
        mat_op:=[];
        for a in arrows_as_paths do
            mat := [];
            for v in vertices do
                if v*a <> Zero(A) then
                    dom_a:= v;
                fi;
            od;
            for v in vertices do
                if a*v <> Zero(A) then
                    im_a:= v;
                fi;
            od; 
        
            pd := Position(vertices,dom_a);
            pi := Position(vertices,im_a);
        
            for m in basis_M{dim_size[pd]} do
                Add(mat,Coefficients(basis_M,m^a){dim_size[pi]}); 
            od;
            Add(mat_op,mat);
        od;

        return mat_op;
    end
);

InstallMethod ( DualOfPathAlgebraMatModule,
# Written by Oeyvind Solberg, 09.06.2010
# Dependent on: MatricesOfPathAlgebraMatModule, RightModuleOverPathAlgebra. 
#
    "for a representation of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebra ], 0,
    function( M, KQ_op )
#
#   M     = a representation of the quiver Q over K with rels
#   KQ_op = opposite algebra of KQ
#
        local num_vert, vertices, basis, dim_vect_M, i, j, count, 
           mat_M, mat_M_op, X, new_op_mat, arrows, a, vert_in_quiver, 
	   origin, target;

#
# The dimension vector for M
#
    num_vert := Length(VerticesOfQuiver(QuiverOfPathAlgebra(KQ_op)));
    vertices := GeneratorsOfAlgebra(RightActingAlgebra(M)){[1..num_vert]};
    basis    := Basis(M); 
    dim_vect_M:=[]; 
    for i in [1.. num_vert] do
    	count:=0;
    	for j in [1..Length(basis)] do
       	    if basis[j]^vertices[i] <> Zero(M) then
               count := count + 1;   
       	    fi;
    	od;
    	Add(dim_vect_M,count);    
    od;

    mat_M    := MatricesOfPathAlgebraMatModule(M); 
    arrows   := ArrowsOfQuiver(QuiverOfPathAlgebra(KQ_op));
    mat_M_op := List(mat_M, X -> TransposedMat(X));
    vert_in_quiver := VerticesOfQuiver(QuiverOfPathAlgebra(KQ_op));
    new_op_mat := [];
    for a in arrows do
    	i := Position(arrows,a); 
	if not ( IsMatrix(mat_M_op[i]) ) then 
	   origin := Position(vert_in_quiver,SourceOfPath(a));
	   target := Position(vert_in_quiver,TargetOfPath(a)); 
	   Add(new_op_mat,[a,[dim_vect_M[origin],dim_vect_M[target]]]); 
	else 
	   Add(new_op_mat,[a,mat_M_op[i]]);
        fi;
    od;

    return RightModuleOverPathAlgebra(KQ_op,new_op_mat);
end
);
