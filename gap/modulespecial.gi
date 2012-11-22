# GAP Implementation
# $Id: specialreps.gi,v 1.12 2012/05/22 14:54:46 sunnyquiver Exp $

#######################################################################
##
#A  IndecProjectiveModules( <A> )
##
##  This function constructs all the indecomposable projective modules
##  over a finite dimensional path algebra or quotient of path algebra. 
##  It returns all the indecomposable projective modules as a list of 
##  modules corresponding to the numbering of the vertices. 
##
InstallMethod ( IndecProjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A ) 
    local fam, KQ, I, Q, num_vert, num_arrows, i, vertices, 
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
    fam := ElementsFamily(FamilyObj(A));
    if HasGroebnerBasisOfIdeal(fam!.ideal) and 
       AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then 
#
#    Finding vertices and arrows as elements of the algebra  A  for later use. 
#
        K := LeftActingDomain(KQ);
        num_arrows := Length(ArrowsOfQuiver(Q)); 
        vertices := GeneratorsOfAlgebra(A){[2..num_vert+1]};
        arrows_as_paths :=  GeneratorsOfAlgebra(A){[num_vert+2..num_arrows+num_vert+1]}; 
#
#
#
        mat_list := [];
        list_of_min_gen := [];
        for p in [1..num_vert] do
            #
            #   Finding a K-basis for indecomposable projective associated with vertex  p.
            #
            P := RightIdeal(A,[vertices[p]]);
            B := CanonicalBasis(P);
            length_B := Length(B);
            intervals_of_basis := List([1..num_vert], x -> []);
            for i in [1..Length(B)] do
                for j in [1..num_vert] do
                    if ( B[i]*vertices[j] <> Zero(A) ) then
                        Add(intervals_of_basis[j],i);
                    fi;
                od;
            od;
            #  
            #   Finding where the indecomposable projective has no support.
            #
            zero_vertices := [];
            for i in [1..num_vert] do
                if ( Length(intervals_of_basis[i]) = 0 ) then
                    Add(zero_vertices,i);
                fi;
            od;
            vertices_Q := VerticesOfQuiver(Q); 
            # 
            #   Finding the matrices defining the representation
            # 
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
            #
            #   Constructing a generator for each indecomposable projective
            #
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
        #
        #   Construct the indecomposable projectives and set the generating set
        # 
        indec_proj_list := [];
        for i in [1..num_vert] do

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
end
);

InstallOtherMethod ( IndecProjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsPathAlgebra ], 0,
    function( A )
    local KQ, rels, I, groebner, groebner_basis, Q, num_vert, num_arrows, 
          i, vertices, arrows_as_paths, j, indec_proj_list, l, arrow, P, 
          gens, length_B, B, intervals_of_basis, mat, a, source, target, 
          vertices_Q, vector, source_index, target_index, mat_index, p, 
          mat_list, zero_vertices, rows, cols, K, partial_mat,
          list_of_min_gen; 
#
#    A = KQ
#
#    Testing input
#       
    Q := QuiverOfPathAlgebra(A); 
    num_vert := Length(VerticesOfQuiver(Q));
    if not IsAcyclicQuiver(Q)  then
        Print("Need to have a finite dimensional path algebra as argument.\n");
        return fail;
    else
#
#    Finding vertices and arrows as elements of the algebra  A  for later use. 
#
        K := LeftActingDomain(A);
        num_vert := Length(VerticesOfQuiver(Q));
        num_arrows := Length(ArrowsOfQuiver(Q)); 
        vertices := GeneratorsOfAlgebra(A){[1..num_vert]};
        arrows_as_paths := GeneratorsOfAlgebra(A){[num_vert+1..num_arrows+num_vert]}; 
        mat_list := [];
        list_of_min_gen := [];
        for p in [1..num_vert] do
            #
            #   Finding a K-basis for indecomposable projective associated with vertex  p.
            #
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
            #  
            #   Finding where the indecomposable projective has no support.
            #
            zero_vertices := [];
            for i in [1..num_vert] do
                if ( Length(intervals_of_basis[i]) = 0 ) then
                    Add(zero_vertices,i);
                fi;
            od;
            vertices_Q := VerticesOfQuiver(Q); 
            # 
            #   Finding the matrices defining the representation
            #             
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
            #
            #   Constructing a generator for each indecomposable projective
            #        
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
        #
        #   Construct the indecomposable projectives and set the generating set
        #    
        indec_proj_list := [];
        for i in [1..Length([1..num_vert])] do
            Add(indec_proj_list,RightModuleOverPathAlgebra(A,mat_list[i]));
            list_of_min_gen[i] := PathModuleElem(FamilyObj(Zero(indec_proj_list[i])![1]),list_of_min_gen[i]);            
            list_of_min_gen[i] := Objectify( TypeObj( Zero(indec_proj_list[i]) ), [ list_of_min_gen[i] ] );
            SetMinimalGeneratingSetOfModule(indec_proj_list[i],[list_of_min_gen[i]]);
        od;
        return indec_proj_list;
    fi;
end
);

#######################################################################
##
#A  IndecInjectiveModules( <A> )
##
##  This function constructs all the indecomposable injective modules
##  over a finite dimensional path algebra or quotient of path algebra. 
##  It returns all the indecomposable injective modules as a list of 
##  modules corresponding to the numbering of the vertices. 
##
InstallMethod ( IndecInjectiveModules, 
    "for a finite dimensional quotient of a path algebra",
    [ IsAlgebra ], 0,
    function( A )
        local A_op, P_op; 

        A_op := OppositeAlgebra(A);
        P_op := IndecProjectiveModules(A_op);
        
        return List(P_op, x -> DualOfModule(x));
    end
);
    
#######################################################################
##
#A  SimpleModules( <A> )
##
##  This function constructs all the simple modules over a finite 
##  dimensional path algebra or quotient of path algebra. It returns 
##  all the simple modules as a list of modules corresponding to the 
##  numbering of the vertices. 
##
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

#######################################################################
##
#A  ZeroModule( <A> )
##
##  This function constructs all the zero module over a finite 
##  dimensional path algebra or quotient of path algebra. 
##
InstallMethod ( ZeroModule, 
    "for a finite dimensional quotient of a path algebra",
    [ IsAlgebra ], 0,
    function( A )

    local KQ, num_vert, zero;
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

#######################################################################
##
#O  BasisOfProjectives(<A>)
##
##  Given a finite dimensional qoutient of a path algebra  <A>, this 
##  function finds in the basis of each indecomposable projective 
##  in terms of paths (nontips of the ideal  I defining  A).
##
InstallMethod ( BasisOfProjectives, 
    "for a finite dimensional quotient of a path algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )

    local Q, num_vert, fam, vertices, basis_of_projs, 
          v, basis_of_proj, P, B, i, j;
#
#    Testing input if finite dimensional.
#       
   fam := ElementsFamily(FamilyObj(A));
   if HasGroebnerBasisOfIdeal(fam!.ideal) and 
          AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then 
#
      Q := QuiverOfPathAlgebra(OriginalPathAlgebra(A));
      num_vert := NumberOfVertices(Q); 
      vertices := List(VerticesOfQuiver(Q), x -> x*One(A));
      basis_of_projs := [];
      for v in vertices do
         basis_of_proj := List([1..num_vert], x -> []);
         P := RightIdeal(A,[v]);
         B := CanonicalBasis(P);
         for i in [1..Length(B)] do
            for j in [1..num_vert] do
               if ( B[i]*vertices[j] <> Zero(A) ) then
                   Add(basis_of_proj[j],B[i]);
               fi;
            od;
         od;
         Add(basis_of_projs,basis_of_proj);
      od;
      return basis_of_projs;
   else
      Print("Need to have a finite dimensional quotient of a path algebra as argument.\n");
      return fail;      
   fi;
end
);

#######################################################################
##
#O  BasisOfProjectives(<A>)
##
##  Given a finite dimensional path algebra  <A>, this function finds 
##  in the basis of each indecomposable projective in terms of paths.
##
InstallOtherMethod ( BasisOfProjectives, 
   "for a finite dimensional quotient of a path algebra",
   [ IsPathAlgebra ], 0,
   function( A )
   
   local Q, num_vert, fam, vertices, basis_of_projs, 
          v, basis_of_proj, P, B, i, j;
#
#  Testing input if finite dimensional.
#       
   Q := QuiverOfPathAlgebra(A); 
   num_vert := NumberOfVertices(Q);
   if not IsAcyclicQuiver(Q)  then
      Print("Need to have a finite dimensional path algebra as argument.\n");
      return fail;
   else
      vertices := List(VerticesOfQuiver(Q), x -> x*One(A));
      basis_of_projs := [];
      for v in vertices do
         basis_of_proj := List([1..num_vert], x -> []);
         P := RightIdeal(A,[v]);
         B := CanonicalBasis(P);
         for i in [1..Length(B)] do
            for j in [1..num_vert] do
               if ( B[i]*vertices[j] <> Zero(A) ) then
                  Add(basis_of_proj[j],B[i]);
               fi;
            od;
         od;
         Add(basis_of_projs,basis_of_proj);
      od;
   fi;
  
   return basis_of_projs;
end
);
