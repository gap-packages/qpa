# GAP Implementation
# This file was generated from 
# $Id: functors.gi,v 1.6 2012/02/27 12:26:34 sunnyquiver Exp $
InstallMethod ( DualOfModule,
    "for a representation of a quiver",
    [ IsPathAlgebraModule ], 0,
    function( M )
#
#   M     = a representation of the quiver Q over K with rels
#
        local A_op, KQ, KQ_op, dim_vect_M, num_vert, vertices, i, j, 
  	      count, mat_M, mat_M_op, X, mats, arrows, a, 
	      vert_in_quiver, origin, target, N;

A_op  := OppositeAlgebra(RightActingAlgebra(M));
dim_vect_M := DimensionVector(M); 
if Sum(dim_vect_M) = 0 then
   return ZeroModule(A_op);
else
    KQ_op := OriginalPathAlgebra(A_op);
    num_vert := Length(DimensionVector(M));
    mat_M    := MatricesOfPathAlgebraModule(M); 
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

    if IsPathAlgebra(A_op) then 
       N:= RightModuleOverPathAlgebra(A_op,mats);
    else 
       N:= RightModuleOverPathAlgebra(A_op,mats);
    fi;
    SetDualOfModule(N,M);
    return N;
fi;
end
);

InstallMethod( TransposeOfModule,
   "for a path algebra",
   [ IsPathAlgebraModule ], 0,
   function( M ) 

   local A, B, Q, fam, KQ, gens, K, 
         num_vert, vertices,  verticesinalg, 
         arrows, BB, i, j, B_M, G, MM, projective_cover, r, 
         test, vertex_list, P, PI_list, zero_in_P, g, temp, 
         matrix, solutions, Solutions, projective_vector, s, 
         a, syzygyoverKQ, U, RG, rgb, new_rtgbofI, run_time, rad, 
         vector, part, newpart, zero, radicalspace, V, W, VW, f, 
         topofsyzygy, x, Topofsyzygy, PP, P_list, P1_list, P_0, P_1, 
         min_gen_list, inc_list, supp, Aop, dim, PPop, P_0op_list, P_0op,
         matrix_op;

   A := RightActingAlgebra(M);
   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
   fam := ElementsFamily(FamilyObj( UnderlyingLeftModule( B ) ));
   KQ := OriginalPathAlgebra(A);
   K := LeftActingDomain(A);

   num_vert := NumberOfVertices(Q);
   vertices := VerticesOfQuiver(Q);
   if IsPathAlgebra(A) then 
      verticesinalg := GeneratorsOfAlgebra(A){[1..num_vert]};
      arrows   := GeneratorsOfAlgebra(A){[1+num_vert..num_vert+NumberOfArrows(Q)]};  
   else 
      verticesinalg := GeneratorsOfAlgebra(A){[2..1+num_vert]};
      arrows   := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+NumberOfArrows(Q)]};
   fi;
#
# Finding a basis for all the indecomposable projective A-modules as
# right ideals.
#
   BB := [];
   for i in [1..num_vert] do 
      BB[i] := [];
   od;
   for i in [1..num_vert] do 
      for j in [1..Length(B)] do 
         if verticesinalg[i]*B[j] <> Zero(A) then
            Add(BB[i],B[j]);
         fi;
      od;
   od;

   B_M := CanonicalBasis(M);
   G   := MinimalGeneratingSetOfModule(M);
#
# Computing the product of all the elements in the 
# minimal generating set G of M, with the basis of the
# indecomposable projective covering that generator.
#
   MM := [];
   projective_cover := [];
   for i in [1..Length(G)] do
      r := 0;
      repeat 
         r := r + 1;
         test := G[i]^verticesinalg[r] <> Zero(M);
      until test;
      projective_cover[i] := r;
      for j in [1..Length(BB[r])] do 
         Add(MM,G[i]^BB[r][j]); 
      od;
   od;

   matrix := NullMat(Length(MM),Dimension(M),K);
   for i in [1..Length(MM)] do
      matrix[i] := Flat(ExtRepOfObj(MM[i])![1]);
   od;
#
#  Finding the kernel of the projective cover as a vectorspace 
#  
   solutions := NullspaceMat(matrix);
#
# Finding the kernel of the projective cover as elements of the
# direct sum of the indecomposable projective modules in the 
# projective cover. 
#
   Solutions := [];
   for i in [1..Length(solutions)] do
      projective_vector := [];
      s := 0;
      for j in [1..Length(projective_cover)] do
         a := Zero(A);
         for r in [1..Length(BB[projective_cover[j]])] do
            a := a + BB[projective_cover[j]][r]*solutions[i][r+s];
         od;
         Add(projective_vector,a);
         s := s + Length(BB[projective_cover[j]]);
      od;
      Add(Solutions,projective_vector);
   od;
#
# Finding the top of the syzygy, first finding generators of the 
# radical of the syzygy, rad.
#
   zero := ListWithIdenticalEntries(Length(projective_cover),Zero(A));
   rad := [];
   for s in Solutions do
      for a in arrows do
         if s*a <> zero then  
            Add(rad,s*a);
         fi; 
      od;
   od;
#
# Computing the radical as a vectorspace by flattening the vectors.
#
   radicalspace := [];
   for s in rad do 
      vector := [];
      for i in [1..Length(projective_cover)] do
         part := Coefficients(B,s[i]);
         newpart := [];
         for j in [1..Length(BB[projective_cover[i]])] do 
            newpart[j] := part[Position(B,BB[projective_cover[i]][j])];
         od;
         Add(vector,newpart);
      od;
      vector := Flat(vector);
      Add(radicalspace,vector);
   od;
#
# Computing a basis for syzygy modulo radical of syzygy as a 
# vectorspace and pulling it back to the syzygy.
#
if Length(solutions) = 0 then 
   Print("The entered module is projective.\n");
   return ZeroModule(OppositeAlgebra(A));
else
   V := VectorSpace(K,solutions);
   W := Subspace(V,radicalspace);
   f := NaturalHomomorphismBySubspace(V,W);
   VW := CanonicalBasis(Range(f));
   topofsyzygy := [];
   for x in VW do
      Add(topofsyzygy,PreImagesRepresentative(f,x));
   od;
#
# Finding the minimal generators of the syzygy as elements of the
# direct sum of the indecomposable projective modules in the 
# projective cover. 
# 
   Topofsyzygy := [];
   for i in [1..Length(topofsyzygy)] do
      projective_vector := [];
      s := 0;
      for j in [1..Length(projective_cover)] do
         a := Zero(A);
         for r in [1..Length(BB[projective_cover[j]])] do
            a := a + BB[projective_cover[j]][r]*topofsyzygy[i][r+s];
         od;
         Add(projective_vector,a);
         s := s + Length(BB[projective_cover[j]]);
      od;
      Add(Topofsyzygy,projective_vector);
   od;
#
# Finding the indecomposable projective modules in the projective 
# cover of M.
#
   PP:= IndecProjectiveModules(A); 
   P_list := [];
   for i in [1..Length(projective_cover)] do
      Add(P_list,PP[projective_cover[i]]);
   od;
   P_0:= DirectSumOfModules(P_list); 
#
# Finding the indecomposable projective modules in the projective
# cover of the syzygy of M.
#
   P1_list := [];
   for i in [1..Length(Topofsyzygy)] do
      supp := [];
      for j in [1..Length(vertices)] do
          if Topofsyzygy[i]*(vertices[j]*One(A)) <> Topofsyzygy[i]*Zero(A) then
             Add(supp,j);
         fi;
      od;      
      if Length(supp) = 1 then 
         Add(P1_list,supp[1]);
      else
         Print("Error: Minimal generators of the syzygy are not uniform.\n");
         return fail;
      fi;
   od;
#
# Transposing the presentation matrix of the module M, 
# and taking the corresponding opposite element in the 
# opposite algebra of all entries in the presentation matrix.
#
   Aop := OppositeAlgebra(A);
   matrix := ShallowCopy(TransposedMat(Topofsyzygy));
   dim := DimensionsMat(matrix);
   for i in [1..dim[1]] do
      matrix[i] := List(matrix[i],x->OppositePathAlgebraElement(x));
   od;
# 
# Finding all indecomposable projective modules over the opposite algebra.
#
   PPop := IndecProjectiveModules(Aop);
#
# Constructing the projective cover of the transpose of M.
#
   P_0op_list := List(P1_list,x->PPop[x]);
   P_0op := DirectSumOfModules(P_0op_list);
#
# Constructing the transpose as the coker of the inclusion of the image
# of the presentation matrix of the transpose of M.
#
   gens := [];
   inc_list := DirectSumInclusions(P_0op);
   min_gen_list := [];
   for i in [1..Length(P_0op_list)] do
      Add(min_gen_list,ImageElm(inc_list[i],MinimalGeneratingSetOfModule(P_0op_list[i])[1]));
   od;
   for i in [1..dim[1]] do
      temp := Zero(P_0op);
      for j in [1..Length(min_gen_list)] do
         temp := temp +  min_gen_list[j]^matrix[i][j];
      od;
      Add(gens,temp);      
   od;
   f := SubRepresentationInclusion(P_0op,gens);
   U := CoKernel(f);

   return U;
fi;
end
);

InstallMethod( DTr,
   "for a path algebra module",
   [ IsPathAlgebraModule ], 0,
   function( M );

   return DualOfModule(TransposeOfModule(M));
end
);

InstallMethod( TrD,
   "for a path algebra module",
   [ IsPathAlgebraModule ], 0,
   function( M );

   return TransposeOfModule(DualOfModule(M));
end
);

InstallOtherMethod( DTr,
   "for a path algebra module",
   [ IsPathAlgebraModule, IS_INT ], 0,
   function( M, n )

   local U, i;

   U := M;
   if n < 0 then 
      for i in [1..-n] do
         Print("Computing step ",i,"...\n");
         U := TrD(U);
      od;
   else 
      if n >= 1 then 
         for i in [1..n] do
            Print("Computing step ",i,"...\n");
            U := DTr(U);
         od;
      fi;
   fi;

   return U;
end
);

InstallOtherMethod( TrD,
   "for a path algebra module",
   [ IsPathAlgebraModule, IS_INT ], 0,
   function( M, n )

   local U, i;

   U := M;
   if n < 0 then 
      for i in [1..-n] do
         Print("Computing step ",i,"...\n");
         U := DTr(U);
      od;
   else 
      if n >= 1 then 
         for i in [1..n] do
         Print("Computing step ",i,"...\n");
            U := TrD(U);
         od;
      fi;
   fi;

   return U;
end
);
