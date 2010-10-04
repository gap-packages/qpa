# GAP Implementation
# $Id: specialreps.gi,v 1.2 2010/10/04 07:07:35 sunnyquiver Exp $

# specialreps.gi: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 

InstallMethod ( IndecomposableProjectiveRepresentations, 
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra, IsList ], 0,
    function( KQ, rels )
        local I, groebner, groebner_basis, basis_A, Q, num_vert, num_arrows, i,
              vertices, arrows_as_paths, indec_proj, j, indec_proj_list,
                indec_proj_rep, l, A, arrow, 
                P, gens, length_B, B, intervals_of_basis, mat, a, source,
                target, vertices_Q, vector, source_index, target_index, mat_index,
                p, mat_list, zero_vertices, rows, cols, K, partial_mat; 
#
#    A = KQ/<rels>
#
#
#    Testing input
#        
        if not (  IsPathAlgebra(KQ) )  then
            Error("Need to enter a path algebra as an argument!");
        fi;
        I := Ideal(KQ,rels);
        groebner := GBNPGroebnerBasis(rels,KQ);
        groebner_basis := GroebnerBasis(I,groebner);
        if not ( AdmitsFinitelyManyNontips(groebner_basis) ) then
            Error("Need to have a finite dimensional quotient of a path algebra.");
        fi;
#
#    Finding vertices and arrows as elements of the algebra  A  for later use. 
#
        Q := QuiverOfPathAlgebra(KQ); 
        K := LeftActingDomain(KQ);
        A := KQ/I;
        num_vert   := Length(VerticesOfQuiver(Q));
        num_arrows := Length(ArrowsOfQuiver(Q)); 
        vertices   := GeneratorsOfAlgebra(A){[2..num_vert+1]};
        arrows_as_paths := 
            GeneratorsOfAlgebra(A){[num_vert+2..num_arrows+num_vert+1]}; 
#
#
#
        mat_list := [];
        for p in [1..num_vert] do
            P := RightIdeal(A,[vertices[p]]);
            B := CanonicalBasis(P);
            length_B := Length(B);
            intervals_of_basis := [];
            for i in [1..num_vert] do
                Add(intervals_of_basis,[]);
            od;
            for i in [1..Length(B)] do
                for j in [1..num_vert] do
                    if ( B[i]*vertices[j] <> Zero(A) ) then
                        Add(intervals_of_basis[j],i);
                    fi;
                od;
            od;
            zero_vertices := [];
            for i in [1..num_vert] do
                if ( Length(intervals_of_basis[i]) = 0 ) then
                    Add(zero_vertices,i);
                fi;
            od;
            vertices_Q := VerticesOfQuiver(Q); 
            mat := [];
            for a in arrows_as_paths do
                partial_mat := [];
                mat_index := Position(arrows_as_paths,a);
                source := SourceOfPath(TipMonomial(a));
                target := TargetOfPath(TipMonomial(a));
                source_index := Position(vertices_Q,source);
                target_index := Position(vertices_Q,target);
                rows := Length(intervals_of_basis[source_index]);
                cols := Length(intervals_of_basis[target_index]);
                arrow := TipMonomial(a);
                if ( rows = 0 or cols = 0 ) then
                    partial_mat := [arrow,[rows,cols]];
                    Add(mat,partial_mat);
                else 
                    for i in intervals_of_basis[source_index] do
                        vector := Coefficients(Basis(P), Basis(P)[i]*a){intervals_of_basis[target_index]};
                        Add(partial_mat,vector);
                    od;
                    Add(mat,[arrow,partial_mat]);
                fi;
            od;
            Add(mat_list,mat);
        od;
        indec_proj_list := [];
        for i in [1..num_vert] do
            Add(indec_proj_list,RightModuleOverPathAlgebra(KQ,mat_list[i]));
        od;
        return indec_proj_list;
    end
);

InstallMethod ( IndecomposableInjectiveRepresentations, 
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra, IsList ], 0,
    function( KQ, rels )
        local Q, K, Q_op, KQ_op, rels_op, P_op, num_vert, indec_inj_list, i; 
#
        Q := QuiverOfPathAlgebra(KQ); 
        KQ_op := OppositePathAlgebra(KQ);
        rels_op := OppositeRelations(rels);
        P_op := IndecomposableProjectiveRepresentations(KQ_op,rels_op);
        num_vert   := Length(VerticesOfQuiver(Q));        
        indec_inj_list := [];
        for i in [1..num_vert] do
            Add(indec_inj_list,DualOfPathAlgebraMatModule(P_op[i]));
        od;

        return indec_inj_list;
    end
);

InstallMethod ( VertexSimpleRepresentations, 
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra, IsList ], 0,
    function( KQ, rels )
        local arrows, vertices, v, mats, a, simple_rep;
#
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(KQ));
    vertices := VerticesOfQuiver(QuiverOfPathAlgebra(KQ));
    simple_rep := [];
    for v in vertices do
        mats := [];
        for a in arrows do
            if SourceOfPath(TipMonomial(One(KQ)*a)) <> v then
                if TargetOfPath(TipMonomial(One(KQ)*a)) <> v then
                    Add(mats,[a,[0,0]]);
                else
                    Add(mats,[a,[0,1]]);
                fi;
            else
                if TargetOfPath(TipMonomial(One(KQ)*a)) <> v then
                    Add(mats,[a,[1,0]]);
                else
                    Add(mats,[a,[[0]]]);
                fi;
            fi;
        od;
        Add(simple_rep,RightModuleOverPathAlgebra(KQ,mats));
    od;
        
    return simple_rep;
end
);

InstallMethod ( ZeroRepresentation, 
# Written by Oeyvind Solberg, 25.09.2010
# Dependent on: ????
#
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra, IsList ], 0,
    function( KQ, rels )
        local arrows, vertices, v, mats, a, simple_rep;
#
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(KQ));
    mats := [];
    for a in arrows do
       Add(mats,[a,[0,0]]);
    od;
        
    return RightModuleOverPathAlgebra(KQ,mats);;
end
);

