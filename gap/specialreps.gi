# GAP Implementation
# $Id: specialreps.gi,v 1.11 2012/05/15 06:58:20 sunnyquiver Exp $

# specialreps.gi: Provides special representations of a quiver, 
# 		  indecomposble projective, indecomposable injective, 
# 		  and vertex simple representations. 
#
InstallMethod ( IndecProjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A ) 
    local which_proj, fam, KQ, I, Q, num_vert, num_arrows, i, vertices, 
          arrows_as_paths, indec_proj, j, indec_proj_list, 
          indec_proj_rep, l, arrow, P, gens, length_B, B, 
          mat, a, source, target, vertices_Q, vector, intervals_of_basis, 
          source_index, target_index, mat_index, p, mat_list, 
          zero_vertices, rows, cols, K, partial_mat, list_of_min_gen; 
#
#    A = KQ/<rels>
#
#    Testing input
#       
  KQ := OriginalPathAlgebra(A);
  Q := QuiverOfPathAlgebra(KQ);
  num_vert := Length(VerticesOfQuiver(Q));
  which_proj := [1..num_vert]; 
  if not ForAll(which_proj, x -> x in [1..num_vert]) then 
     Print("The range of projectives entered is wrong.");
     return fail;
  else
     fam := ElementsFamily(FamilyObj(A));
     if HasGroebnerBasisOfIdeal(fam!.ideal) and 
          AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then 

#
#    Finding vertices and arrows as elements of the algebra  A  for later use. 
#

        K := LeftActingDomain(KQ);

        num_arrows := Length(ArrowsOfQuiver(Q)); 
        vertices   := GeneratorsOfAlgebra(A){[2..num_vert+1]};
        arrows_as_paths := 
            GeneratorsOfAlgebra(A){[num_vert+2..num_arrows+num_vert+1]}; 
#
#
#
        mat_list := [];
        list_of_min_gen := [];
        for p in which_proj do
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
            for i in [1..Length(intervals_of_basis)] do
               if intervals_of_basis[i] = [] then
                  intervals_of_basis[i] := [Zero(K)];
               else 
                  for j in [1..Length(intervals_of_basis[i])] do
                     if intervals_of_basis[i][j] > 1 then 
                        intervals_of_basis[i][j] := Zero(K);
                     else
                        intervals_of_basis[i][j] := One(K);
                     fi;
                  od;
               fi;
            od;
            Add(list_of_min_gen,intervals_of_basis);
        od;
        indec_proj_list := [];
        for i in [1..Length(which_proj)] do
            Add(indec_proj_list,RightModuleOverPathAlgebra(A,mat_list[i]));
            list_of_min_gen[i] := PathModuleElem(FamilyObj(Zero(indec_proj_list[i])![1]),list_of_min_gen[i]); 
            
            list_of_min_gen[i] := Objectify( TypeObj( Zero(indec_proj_list[i]) ), [ list_of_min_gen[i] ] );
            SetMinimalGeneratingSetOfModule(indec_proj_list[i],[list_of_min_gen[i]]);

        od;
        return indec_proj_list;
     else
        Print("Compute a Groebner basis of the ideal you are factoring out with before you form the quotient algebra, or you have entered an algebra which is not finite dimensional.\n");
        return fail;
     fi;
  fi;
end
);

InstallOtherMethod ( IndecProjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra ], 0,
    function( A )
        local which_proj, KQ, rels, I, groebner, groebner_basis, Q, num_vert, 
              num_arrows, i, vertices, arrows_as_paths, j, 
              indec_proj_list, l, arrow, P, gens, 
              length_B, B, intervals_of_basis, mat, a, source, target, 
              vertices_Q, vector, source_index, target_index, mat_index,
              p, mat_list, zero_vertices, rows, cols, K, partial_mat,
              list_of_min_gen; 
#
#    A = KQ
#
#    Testing input
#       
     Q := QuiverOfPathAlgebra(A); 
     num_vert := Length(VerticesOfQuiver(Q));
     which_proj := [1..num_vert];
     if not ForAll(which_proj, x -> x in [1..num_vert]) then 
        Print("The range of projectives entered is wrong.\n");
        return fail;
     elif not IsAcyclicQuiver(Q)  then
        Print("Need to have a finite dimensional path algebra as argument.\n");
        return fail;
     else
#
#    Finding vertices and arrows as elements of the algebra  A  for later use. 
#
        K := LeftActingDomain(A);
        num_vert   := Length(VerticesOfQuiver(Q));
        num_arrows := Length(ArrowsOfQuiver(Q)); 
        vertices   := GeneratorsOfAlgebra(A){[1..num_vert]};
        arrows_as_paths := 
            GeneratorsOfAlgebra(A){[num_vert+1..num_arrows+num_vert]}; 
#
#
#
        mat_list := [];
        list_of_min_gen := [];
        for p in which_proj do
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
            for i in [1..Length(intervals_of_basis)] do
               if intervals_of_basis[i] = [] then
                  intervals_of_basis[i] := [Zero(K)];
               else 
                  for j in [1..Length(intervals_of_basis[i])] do
                     if intervals_of_basis[i][j] > 1 then 
                        intervals_of_basis[i][j] := Zero(K);
                     else
                        intervals_of_basis[i][j] := One(K);
                     fi;
                  od;
               fi;
            od;
            Add(list_of_min_gen,intervals_of_basis);
        od;
        indec_proj_list := [];
        for i in [1..Length(which_proj)] do
            Add(indec_proj_list,RightModuleOverPathAlgebra(A,mat_list[i]));
            list_of_min_gen[i] := PathModuleElem(FamilyObj(Zero(indec_proj_list[i])![1]),list_of_min_gen[i]); 
            
            list_of_min_gen[i] := Objectify( TypeObj( Zero(indec_proj_list[i]) ), [ list_of_min_gen[i] ] );
            SetMinimalGeneratingSetOfModule(indec_proj_list[i],[list_of_min_gen[i]]);
        od;
        return indec_proj_list;
     fi;
end
);

InstallMethod ( IndecInjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsAlgebra ], 0,
    function( A )
        local A_op, P_op, Q, num_vert, indec_inj_list, i; 


        A_op := OppositeAlgebra(A);
        Q    := QuiverOfPathAlgebra(A); 
        num_vert := Length(VerticesOfQuiver(Q));        
        P_op := IndecProjectiveModules(A_op);
        indec_inj_list := List([1..num_vert], x -> DualOfModule(P_op[i]));

        return indec_inj_list;
    end
);

InstallMethod ( SimpleModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsAlgebra ], 0,
    function( A )

    local KQ, num_vert, simple_rep, zero, v, temp;
#
    if ( not IsPathAlgebra(A) ) and ( not IsQuotientOfPathAlgebra(A) ) then 
       Error("argument entered is a not (a quotient of) a path algebra,\n");
    fi;
 
    KQ := OriginalPathAlgebra(A); 
    num_vert := NumberOfVertices(QuiverOfPathAlgebra(A));
    simple_rep := [];
    zero := List([1..num_vert], x -> 0);
    for v in [1..num_vert] do
    	temp := ShallowCopy(zero);
        temp[v] := 1; 
        Add(simple_rep,RightModuleOverPathAlgebra(A,temp,[]));
    od;
        
    return simple_rep;
end
);

InstallMethod ( ZeroModule, 
    "for a finite dimensional quotient of a path algebra",
    [ IsAlgebra ], 0,
    function( A )

    local KQ, num_vert, simple_rep, zero, v, temp;
#
    if ( not IsPathAlgebra(A) ) and ( not IsQuotientOfPathAlgebra(A) ) then 
       Error("argument entered is a not (a quotient of) a path algebra,\n");
    fi;
 
    KQ := OriginalPathAlgebra(A); 
    num_vert := NumberOfVertices(QuiverOfPathAlgebra(A));
    zero := List([1..num_vert], x -> 0);

    return RightModuleOverPathAlgebra(A,zero,[]);
end
);

