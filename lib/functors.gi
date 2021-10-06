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
	   Add( mats, [ String( a ), mat_M_op[ i ] ] );
	fi;
    od;
    N := RightModuleOverPathAlgebra( A_op, dim_vect_M, mats );

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
                Basis(Subspace(FullRowSpace(K, Length(b[1])), b, "basis"), b)); 
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
                Basis(Subspace(FullRowSpace(K, Length(b[1])), b, "basis"), b)); 
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
        if Length(BfVN[i]) <> 0 and Length( BasisVMflat[ i ] ) <> 0 then 
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
    
    return DualOfModule( StarOfModule( M ) );
end
  );

#######################################################################
##
#A  OppositeNakayamaFunctorOfModule( <M> )
##
##  This function takes as an argument a module over an algebra  A  and
##  computes module Hom_A(Hom_K( M, K ), A )  over  A. 
##  
InstallMethod( OppositeNakayamaFunctorOfModule, 
    "for a matrix",
    [ IsPathAlgebraMatModule ],
        
    function( M );
    
    return StarOfModule( DualOfModule( M ) );
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
    
    return DualOfModuleHomomorphism( StarOfModuleHomomorphism( f ) );
end
);

#######################################################################
##
#A  OppositeNakayamaFunctorOfModuleHomomorphism( <f> )
##
##  This function takes as an argument a homomorphism  f  between two 
##  modules  M  and  N  over an algebra  A  and computes the induced 
##  homomorphism from the module  Hom_K(Hom_A(N, A), K)  to the module  
##  Hom_K(Hom_A(M, A), K)  over  A. 
##
InstallMethod( OppositeNakayamaFunctorOfModuleHomomorphism, 
    "for a matrix",
    [ IsPathAlgebraMatModuleHomomorphism ],
        
    function( f ); 
    
    return StarOfModuleHomomorphism( DualOfModuleHomomorphism( f ) );
end
);

#######################################################################
##
#O  TensorProductOfModules( <M>, <N> )
##
##  Given two representations  <Arg>M</Arg> and <Arg>N</Arg>, where
##  <Arg>M</Arg> is a right module over  <M>A</M> and  <Arg>N</Arg> is
##  a right module over the opposite of <M>A</M>, then this function
##  computes the tensor product <M>M\otimes_A N</M> as a vector space
##  and a function <M>M\times N \to M\otimes_A N</M>.
##
InstallMethod ( TensorProductOfModules, 
    "for two representations of a fin. dim. quotient of a path algebra",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 
    0,
    function( M, N )

  local   A,  F,  topofM,  dimvectorN,  V,  elementarytensor,  
          projres,  pi,  d1,  inclusionsP1,  projectionsP0,  
          imagesofgensP1,  t,  projectivesinP0,  verticesofquiver,  
          projectivesinP1,  matrix,  topofP1,  topofP0,  dimvectorM,  
          size_matrix,  bigmatrix,  i,  tempmat,  j,  m,  s,  startN,  
          temp,  d1tensoronemap,  rowsinblock,  r,  dim,  W,  h;
        
   A := RightActingAlgebra( M );
   if A <> OppositeAlgebra( RightActingAlgebra( N ) ) then
       Error("The entered modules are not over the appropriate algebras,\n");
   fi;
   F := LeftActingDomain( M );
#
# An easy check if the tensor product is zero, if so, return zero.
#
   topofM := DimensionVector( TopOfModule( M ) );
   dimvectorN := DimensionVector( N );
   if topofM*dimvectorN = 0 then
     V := F^1;
     elementarytensor := function( m, n ); 
        return Zero(V); 
     end;

     return [ TrivialSubspace( V ), elementarytensor ];  
   fi;
#
# Computing M \otimes N via taking a projective resolution of M and then tensoring
# with N. 
#
   projres := ProjectiveResolution( M );
   ObjectOfComplex( projres, 1 );
   pi := DifferentialOfComplex( projres, 0 ); 
   d1 := DifferentialOfComplex( projres, 1 ); 
   
   projectionsP0 := DirectSumProjections( Range( d1 ) );
   t := Length( projectionsP0 );
   projectivesinP0 := List( projectionsP0, p -> SupportModuleElement( MinimalGeneratingSetOfModule( Range( p ) )[ 1 ] ) );
   if IsQuotientOfPathAlgebra( A ) then 
       projectivesinP0 := List( projectivesinP0, t -> CoefficientsAndMagmaElements( t![ 1 ]![ 1 ] )[ 1 ] );
   elif IsPathAlgebra( A ) then
       projectivesinP0 := List( projectivesinP0, t -> CoefficientsAndMagmaElements( t![ 1 ] )[ 1 ] );
   else
       Error("Something went wrong, inform the maintainers of the QPA-package,\n");
   fi;  
   verticesofquiver := VerticesOfQuiver( QuiverOfPathAlgebra( A ) );
   projectivesinP0 := List( projectivesinP0, v -> Position( verticesofquiver, v ) );
#
# If P_1 - d1 -> P_0 --> M --> 0 is a minimal projective presentation of M and 
# P_1\otimes N = (0), then  M\otimes N \simeq P_0\otimes N.
#
   if DimensionVector( Source( d1 ) ) * dimvectorN = 0 then
       dim := topofM * dimvectorN;
       V := F^dim;
       h := NaturalHomomorphismBySubspace( V, TrivialSubspace( V ) );
       elementarytensor := function( m, n ) 
           local  p0, psum, nsum, v;
           
           p0 := PreImagesRepresentative( pi, m );
           psum := List( [ 1..t ], n -> ElementInIndecProjective( A, ImageElm( projectionsP0[ n ], p0 ), projectivesinP0[ n ] ) );
           if IsQuotientOfPathAlgebra( A ) then
               psum := List( psum, x -> ElementOfQuotientOfPathAlgebra( ElementsFamily( FamilyObj( OppositeAlgebra( A ) ) ), OppositePathAlgebraElement( x ), true ) );
           elif IsPathAlgebra( A ) then
               psum := List( psum, x -> OppositePathAlgebraElement( x ) );
           else
               Error("Something went wrong, inform the maintainers of the QPA-package,\n");
           fi;
           nsum := List( psum, x -> n^x );
           v := Flat( List( [ 1..size_matrix[ 2 ] ], l -> nsum[ l ]![ 1 ]![ 1 ][ projectivesinP0[ l ] ] ) );
                
           return ImageElm( h, v );
       end;
       
       return [ Range( h ), elementarytensor ];
   fi;
   
   inclusionsP1 := DirectSumInclusions( Source( d1 ) );   
   imagesofgensP1 := List( inclusionsP1, f -> ImageElm( f, MinimalGeneratingSetOfModule( Source( f ) )[ 1 ] ) );
   imagesofgensP1 := List( imagesofgensP1, m -> ImageElm( d1, m ) );
   projectivesinP1 := List( inclusionsP1, i -> SupportModuleElement( MinimalGeneratingSetOfModule( Source( i ) )[ 1 ] ) );
   projectivesinP1 := List( projectivesinP1, t -> CoefficientsAndMagmaElements( t![ 1 ]![ 1 ] )[ 1 ] );
   projectivesinP1 := List( projectivesinP1, v -> Position( verticesofquiver, v ) );
   
   matrix := List( imagesofgensP1, m -> List( [ 1..t ], n ->  
                                                          ElementInIndecProjective( A, ImageElm( projectionsP0[ n ], m ), projectivesinP0[ n ] ) ) ); 
   if IsQuotientOfPathAlgebra( A ) then 
       matrix := List( matrix, m -> List( m, x -> ElementOfQuotientOfPathAlgebra( ElementsFamily( FamilyObj( OppositeAlgebra( A ) ) ), OppositePathAlgebraElement( x ), true ) ) );
   elif IsPathAlgebra( A ) then
       matrix := List( matrix, m -> List( m, x -> OppositePathAlgebraElement( x ) ) );
   else
       Error("Something went wrong, inform the maintainers of the QPA-package,\n");
   fi;
   
   topofP1 := DimensionVector( TopOfModule( Source( d1 ) ) );
   topofP0 := DimensionVector( TopOfModule( Range( d1 ) ) );
   dimvectorM := DimensionVector( M );
   size_matrix := DimensionsMat( matrix );
   
   bigmatrix := NullMat( size_matrix[ 1 ], size_matrix[ 2 ], F );   
   for i in [ 1..size_matrix[ 1 ] ] do
       for j in [ 1..size_matrix[ 2 ] ] do
           m := matrix[ i ][ j ];
           s := projectivesinP0[ j ];
           t := projectivesinP1[ i ];
           startN := Sum( dimvectorN{ [ 1..t - 1 ] } );
           if Length( Basis( N ){ [ startN + 1..startN + dimvectorN[ t ] ] } ) > 0 then 
               temp := List( Basis( N ){ [ startN + 1..startN + dimvectorN[ t ] ] }, b -> b^m );
               tempmat := List( temp, t -> t![ 1 ]![ 1 ][ s ] );
           else
               tempmat := NullMat( 1, Length( Basis( N )[ 1 ]![ 1 ]![ 1 ][ s ] ), F );
           fi;
           bigmatrix[ i ][ j ] := tempmat;
       od;
   od;
   
   d1tensoronemap := [];
   for i in [ 1..size_matrix[ 1 ] ] do
       rowsinblock := DimensionsMat( bigmatrix[ i ][ 1 ] )[ 1 ];
       for r in [ 1..rowsinblock ] do      
           temp := [];
           for j in [ 1..size_matrix[ 2 ] ] do
               Add( temp, bigmatrix[ i ][ j ][ r ] );
           od;
           temp := Flat(temp);
           Add( d1tensoronemap, temp );
       od;
   od;

   dim := Length( d1tensoronemap[ 1 ] );
   V := F^dim;
   W := Subspace( V, d1tensoronemap );
   h := NaturalHomomorphismBySubspace( V, W );
   
   elementarytensor := function( m, n ) 
       local  p0, psum, nsum, v;
       
       p0 := PreImagesRepresentative( pi, m );
       psum := List( [ 1..t ], n -> ElementInIndecProjective( A, ImageElm( projectionsP0[ n ], p0 ), projectivesinP0[ n ] ) );
       if IsQuotientOfPathAlgebra( A ) then
           psum := List( psum, x -> ElementOfQuotientOfPathAlgebra( ElementsFamily( FamilyObj( OppositeAlgebra( A ) ) ), OppositePathAlgebraElement( x ), true ) );
       elif IsPathAlgebra( A ) then
           psum := List( psum, x -> OppositePathAlgebraElement( x ) );
       else
           Error("Something went wrong, inform the maintainers of the QPA-package,\n");
       fi;
       nsum := List( psum, x -> n^x );
       v := Flat( List( [ 1..size_matrix[ 2 ] ], l -> nsum[ l ]![ 1 ]![ 1 ][ projectivesinP0[ l ] ] ) );
                
       return ImageElm( h, v );
   end;
   
   return [ Range( h ), elementarytensor ];
end);

InstallMethod( TransposeOfModuleHomomorphism,
    "for a homomorphism a path algebra module",
    [ IsPathAlgebraMatModuleHomomorphism ], 0,
    function( f ) 

    local   B,  C,  d0,  d0prime,  kerd0,  kerd0prime,  d1temp,  
            d1primetemp,  d1,  d1prime,  f0,  f1prime,  f1,  homstar,  
            Dhomstar,  h;

    B := Source( f );
    C := Range( f );
    d0 := ProjectiveCover( B );
    d0prime := ProjectiveCover( C );
    kerd0 := KernelInclusion( d0 );
    kerd0prime := KernelInclusion( d0prime );
    d1temp := ProjectiveCover( Source( kerd0 ) ); 
    d1primetemp := ProjectiveCover( Source( kerd0prime ) );
    d1 := d1temp * kerd0;
    d1prime := d1primetemp * kerd0prime;
    f0 := LiftingMorphismFromProjective( d0prime, d0 * f );
    f1prime := MorphismOnKernel( d0, d0prime, f0, f );
    f1 := LiftingMorphismFromProjective( d1primetemp, d1temp * f1prime );
    
    homstar := List( [ d1, d1prime, f1, f0 ], StarOfModuleHomomorphism );
    h := MorphismOnCoKernel( homstar[ 2 ], homstar[ 1 ], homstar[ 4 ], homstar[ 3 ] );
    
    return h;
end
  );

InstallMethod( TrD,
    "for a homomorphism a path algebra module",
    [ IsPathAlgebraMatModuleHomomorphism ], 0,
    function( f )
    
    local   h,  B,  C;
    
    h := TransposeOfModuleHomomorphism( DualOfModuleHomomorphism( f ) );
    B := TrD( Source( f ) );
    C := TrD( Range( f ) );
    h := RightModuleHomOverAlgebra( B, C, h!.maps );
    
    return h;
end
  );

InstallMethod( DTr,
    "for a homomorphism a path algebra module",
    [ IsPathAlgebraMatModuleHomomorphism ], 0,
    function( f )
    
    local   h,  B,  C;
    
    h := DualOfModuleHomomorphism( TransposeOfModuleHomomorphism( f ) );
    B := DTr( Source( f ) );
    C := DTr( Range( f ) );
    h := RightModuleHomOverAlgebra( B, C, h!.maps );
    
    return h;

    return DualOfModuleHomomorphism( TransposeOfModuleHomomorphism( f ) );
end
  );

#######################################################################
##
#O  DTr( <f>, <n> )
##
##  This function returns returns DTr^n( <f> ) of a module homomorphism <f>.
##  It will also print "Computing step i ..." when computing the ith 
##  DTr. <n> must be an integer.
##  
InstallOtherMethod( DTr,
   "for a path algebra module homomorphism",
   [ IsPathAlgebraMatModuleHomomorphism, IS_INT ], 0,
   function( f, n )

   local g, i;

   g := f;
   if n < 0 then 
      for i in [ 1..-n ] do
         Print( "Computing step ",i,"...\n" );
         g := TrD( g );
      od;
   else 
      if n >= 1 then 
         for i in [ 1..n ] do
            Print( "Computing step ",i,"...\n" );
            g := DTr( g );
         od;
      fi;
   fi;

   return g;
end
  );

#######################################################################
##
#O  TrD( <f>, <n> )
##  
##  This function returns returns TrD^n( <f> ) of a module homomorphism <f>.
##  It will also print "Computing step i ..." when computing the ith 
##  TrD. <n> must be an integer.
##  
##
InstallOtherMethod( TrD,
   "for a path algebra module homomorphism",
   [ IsPathAlgebraMatModuleHomomorphism, IS_INT ], 0,
   function( f, n )

   local g, i;

   g := f;
   if n < 0 then 
      for i in [ 1..-n ] do
         Print( "Computing step ",i,"...\n" );
         g := DTr( g );
      od;
   else 
      if n >= 1 then 
         for i in [ 1..n ] do
         Print( "Computing step ",i,"...\n" );
            g := TrD( g );
         od;
      fi;
   fi;

   return g;
end
  ); 