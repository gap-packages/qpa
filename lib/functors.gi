# GAP Implementation
# This file was generated from 
# $Id: functors.gi,v 1.9 2012/05/15 07:12:00 sunnyquiver Exp $


#######################################################################
##
#O  DualOfModule(<M>)
##
##  This function computes the dual DM of the module <M>, that is, 
##  it computes Hom_k(M, k). If <M> is a module over A, then DM is a
##  module over the opposite algebra of A.
##
InstallMethod ( DualOfModule,
    "for a representation of a quiver",
    [ IsPathAlgebraMatModule ], 0,
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
    if HasIsIndecomposableModule(M) and IsIndecomposableModule(M) then
        SetIsIndecomposableModule(N, true); 
    fi;
    if HasIsProjectiveModule(M) and IsProjectiveModule(M) then
        SetIsInjectiveModule(N, true); 
    fi;
    if HasIsInjectiveModule(M) and IsInjectiveModule(M) then
        SetIsProjectiveModule(N, true); 
    fi;
    return N;
fi;
end
);

#######################################################################
##
#O  TransposeOfModule(<M>)
##
##  This function computes the transpose Tr(M) of the module <M>, that is, 
##  if  P1 ---> P0 ---> M ---> 0 is a projective presentation of <M>,
##  and <M> is a module over an algebra A, then Tr(M) is the cokernel
##  of the map
##
##    Hom_A(P0,A) ---> Hom_A(P1,A)
##
InstallMethod( TransposeOfModule,
   "for a path algebra",
   [ IsPathAlgebraMatModule ], 0,
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
   
   if IsProjectiveModule(M) then
       return ZeroModule(OppositeAlgebra(A));
   fi;
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
      Add(P_list,ShallowCopy(PP[projective_cover[i]]));
   od;
   P_0:= DirectSumOfQPAModules(P_list); 
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
   P_0op_list := List(P1_list,x -> ShallowCopy(PPop[x]));
   P_0op := DirectSumOfQPAModules(P_0op_list);
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
   if Dimension(U) <> 0 and HasIsIndecomposableModule(M) and IsIndecomposableModule(M) then
       SetIsIndecomposableModule(U, true);
   fi;
   
   return U;
fi;
end
);

#######################################################################
##
#O  DTr(<M>)
##
##  This function returns the dual of the transpose of a module <M>,
##  sometimes denoted by "tau". It uses the previous methods DualOfModule
##  and TransposeOfModule.
##  
InstallMethod( DTr,
    "for a path algebra module",
    [ IsPathAlgebraMatModule ], 0,
    function( M );

    return DualOfModule(TransposeOfModule(M));
end
);

#######################################################################
##
#O  TrD(<M>)
##
##  This function returns the transpose of the dual of a module <M>,
##  sometimes denoted by "tau inverse". It uses the previous methods DualOfModule
##  and TransposeOfModule.
##  
InstallMethod( TrD,
    "for a path algebra module",
    [ IsPathAlgebraMatModule ], 0,
    function( M );

    return TransposeOfModule(DualOfModule(M));
end
);

#######################################################################
##
#O  DTr(<M>, <n>)
##
##  This function returns returns DTr^n(<M>) of a module <M>.
##  It will also print "Computing step i ..." when computing the ith 
##  DTr. <n> must be an integer.
##  
InstallOtherMethod( DTr,
   "for a path algebra module",
   [ IsPathAlgebraMatModule, IS_INT ], 0,
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

#######################################################################
##
#O  TrD(<M>, <n>)
##  
##  This function returns returns TrD^n(<M>) of a module <M>.
##  It will also print "Computing step i ..." when computing the ith 
##  TrD. <n> must be an integer.
##  
##
InstallOtherMethod( TrD,
   "for a path algebra module",
   [ IsPathAlgebraMatModule, IS_INT ], 0,
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

#######################################################################
##
#A  StarOfModule( <M> )
##
##  This function takes as an argument a module over an algebra  A  and
##  computes the module  Hom_A(M,A)  over the opposite of  A. 
##  
InstallMethod( StarOfModule, 
    "for a matrix",
    [ IsPathAlgebraMatModule ],
        
    function( M )
    local A, Aop, K, P, V, BV, BasisVflat, b, dimvector, Bproj, arrows, 
          arrowsop, leftmultwitharrows, a, source, target, matrices, temp;
    # 
    # Setting up necessary infrastructure.
    # 
    A := RightActingAlgebra(M);
    Aop := OppositeAlgebra(A);
    K := LeftActingDomain(M); 
    P := IndecProjectiveModules(A);
    #
    # Finding the vector spaces in each vertex in the representation of 
    # Hom_A(M,A)  and making a vector space of the flatten matrices in 
    # each basis vector in  Hom_A(M,e_iA).
    #
    V := List(P, p -> HomOverAlgebra(M,p));
    BV := List(V, v -> List(v, x -> 
                  Flat(MatricesOfPathAlgebraMatModuleHomomorphism(x))));
    BasisVflat := [];
    for b in BV do
        if Length(b) <> 0 then 
            Add(BasisVflat, 
                Basis(Subspace(FullRowSpace(K, Length(b[1])), b, "basis"))); 
        else 
            Add(BasisVflat, []);
        fi;
    od;
    # 
    # Finding the dimension vector of  Hom_A(M,A)  and computing the maps 
    # induced by left multiplication by arrows  alpha : i ---> j  from 
    # e_jA ---> e_iA.
    #
    dimvector := List(V, Length);
    Bproj := BasisOfProjectives(A);
    Bproj := List([1..Length(dimvector)], 
                  i -> Subspace(A, Flat(Bproj[i]), "basis"));
    arrows := ArrowsOfQuiver(QuiverOfPathAlgebra(A));
    arrowsop := ArrowsOfQuiver(QuiverOfPathAlgebra(Aop));    
    leftmultwitharrows := [];
    for a in arrows do
        source := VertexIndex(SourceOfPath(a));
        target := VertexIndex(TargetOfPath(a));
        Add(leftmultwitharrows, 
            HomFromProjective(MinimalGeneratingSetOfModule(P[source])[1]^(a*One(A)), P[source]));
    od;
    #
    # Finally, finding the linear maps in the representation given by  Hom_A(M,A)  and 
    # constructing the module/representation over  A^op.
    # 
    matrices := [];
    for a in arrows do
        source := VertexIndex(SourceOfPath(a));
        target := VertexIndex(TargetOfPath(a));
        if dimvector[source] <> 0 and dimvector[target] <> 0 then 
            temp := V[target]*leftmultwitharrows[Position(arrows, a)];
            temp := List(temp, t -> Coefficients(BasisVflat[source], 
                            Flat(MatricesOfPathAlgebraMatModuleHomomorphism(t))));
            Add(matrices, [String(arrowsop[Position(arrows, a)]), temp]);
        fi;
    od;
    
    return RightModuleOverPathAlgebra(Aop, dimvector, matrices);    
end
);
    
#######################################################################
##
#A  StarOfModuleHomomorphism( <f> )
##
##  This function takes as an argument a homomorphism  f  between two 
##  modules  M  and  N  over an algebra  A  and computes the induced 
##  homomorphism from the module  Hom_A(N, A)  to the module  
##  Hom_A(M, A)  over the opposite of  A. 
##
InstallMethod( StarOfModuleHomomorphism, 
    "for a matrix",
    [ IsPathAlgebraMatModuleHomomorphism ],
        
    function( f )
    local M, N, starM, starN, A, Aop, K, P, VN, dimvectorstarN, BfVN, 
          VM, dimvectorstarM, BVM, BasisVMflat, b, maps, i, map, dimrow, 
          dimcol;
    
    M := Source(f);
    N := Range(f);
    # 
    # Finding the source and range of  StarOfModuleHomomorphism(f).
    #
    starM := StarOfModule(M);
    starN := StarOfModule(N); 
    #
    # Setting up necessary infrastructure.
    #
    A := RightActingAlgebra(M);
    Aop := OppositeAlgebra(A);
    K := LeftActingDomain(M); 
    P := IndecProjectiveModules(A);
    VN := List(P, p -> HomOverAlgebra(N,p));
    dimvectorstarN := List(VN, Length); 
    #
    # Finding the image of the map  StarOfModuleHomomorphism(f).
    # 
    VN := List(VN, v -> f*v);
    BfVN := List(VN, v -> List(v, x -> 
                    Flat(MatricesOfPathAlgebraMatModuleHomomorphism(x))));
    #
    # Computing the basis of StarOfModule(Source(f))  used when the module 
    # was constructed. 
    # 
    VM := List(P, p -> HomOverAlgebra(M,p));
    dimvectorstarM := List(VM, Length); 
    BVM := List(VM, v -> List(v, x -> 
                   Flat(MatricesOfPathAlgebraMatModuleHomomorphism(x))));
    BasisVMflat := [];
    for b in BVM do
        if Length(b) <> 0 then 
            Add(BasisVMflat, 
                Basis(Subspace(FullRowSpace(K, Length(b[1])), b, "basis"))); 
        else 
            Add(BasisVMflat, []);
        fi;
    od;
    #
    # Constructing the map induced by  f  from StarOfModule(Range(f))  to 
    # StarOfModule(Source(f)). 
    # 
    maps := [];
    for i in [1..Length(BfVN)] do
        if Length(BfVN[i]) <> 0 then 
            map := List(BfVN[i], t -> Coefficients(BasisVMflat[i], t));
        else
            if dimvectorstarN[i] = 0 then 
                dimrow := 1;
            else 
                dimrow := dimvectorstarN[i];
            fi;
            if dimvectorstarM[i] = 0 then
                dimcol := 1;
            else 
                dimcol := dimvectorstarM[i];
            fi;
            map := NullMat(dimrow,dimcol,K);
        fi;
        Add(maps, map);
    od;
    
    return RightModuleHomOverAlgebra(starN, starM, maps);
end
);

#######################################################################
##
#A  NakayamaFunctorOfModule( <M> )
##
##  This function takes as an argument a module over an algebra  A  and
##  computes the image of the Nakayama functor, that is, the module  
##  Hom_K(Hom_A(M, A), K)  over  A. 
##  
InstallMethod( NakayamaFunctorOfModule, 
    "for a matrix",
    [ IsPathAlgebraMatModule ],
        
    function( M );
    
    return DualOfModule(StarOfModule(M));
end
  );

#######################################################################
##
#A  NakayamaFunctorOfModuleHomomorphism( <f> )
##
##  This function takes as an argument a homomorphism  f  between two 
##  modules  M  and  N  over an algebra  A  and computes the induced 
##  homomorphism from the module  Hom_K(Hom_A(N, A), K)  to the module  
##  Hom_K(Hom_A(M, A), K)  over  A. 
##
InstallMethod( NakayamaFunctorOfModuleHomomorphism, 
    "for a matrix",
    [ IsPathAlgebraMatModuleHomomorphism ],
        
    function( f ); 
    
    return DualOfModuleHomomorphism(StarOfModuleHomomorphism(f));
end
);