#######################################################################
##
#O  PushOut(<f>, <g>)
##
##  This function finds the pushout of the homomorphisms 
##                    f
##              A ---------> B 
##              |            |
##              | g          | g'
##              |            |
##              V    f'      V
##              C ---------> E
##  in that it returns the homomorphisms  [f', g']. 
##
InstallMethod( PushOut,
   "for two homomorphisms starting in a common module over a path algebra",
   [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 0,
   function( f, g ) 

   local B, C, CplusB, inclusions, projections, h;

   if Source(f) = Source(g) then 
      B := Range(f);
      C := Range(g);
      CplusB := DirectSumOfQPAModules([C,B]);
      inclusions := DirectSumInclusions(CplusB);
      h := f*inclusions[2]-g*inclusions[1];
      h := CoKernelProjection(h);

      return [inclusions[1]*h,inclusions[2]*h];
   else
      Print("Error: The two maps entered don't start in the same module.\n");
      return fail;
   fi;
end
);

#######################################################################
##
#O  PullBack(<f>, <g>)
##
##  This function finds the pullback of the homomorphisms 
##                    f'
##              E ---------> C 
##              |            |
##              | g'         | g
##              |            |
##              V     f      V
##              A ---------> B
##  in that it returns the homomorphisms  [f', g']. 
##
InstallMethod( PullBack,
   "for two homomorphisms ending in a common module over a path algebra",
   [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ], 0,
   function( f, g ) 

   local A, C, AplusC, projections, h;

   if Range(f) = Range(g) then 
      A := Source(f);
      C := Source(g);
      AplusC := DirectSumOfQPAModules([A,C]);
      projections := DirectSumProjections(AplusC);
      h := projections[1]*f-projections[2]*g;
      h := KernelInclusion(h);

      return [h*projections[2],h*projections[1]];
   else
      Print("Error: The two maps entered don't end in the same module.\n");
      return fail;
   fi;
end
);

#######################################################################
##
#O  IsOmegaPeriodic( <M>, <n> )
##
##  This function tests if the module  <M>  is \Omega-periodic, that is,
##  if  M \simeq \Omega^i(M)  when  i  ranges over the set {1,2,...,n}.
##  Otherwise it returns false.
##
InstallMethod( IsOmegaPeriodic, 
   "for a path algebra matmodule and an integer",
   [ IsPathAlgebraMatModule, IS_INT  ], 0,
   function( M, n ) 

   local N0, N1, i;
 
   N0 := M;
   for i in [1..n] do
      Print("Computing syzygy number: ",i,"\n");
      N1 := 1stSyzygy(N0);
      if IsomorphicModules(M,N1) then
         return i;
      else
         N0 := N1;
      fi;
   od;
   return false;
end
);

#######################################################################
##
#O  IsTauPeriodic( <M>, <n> )
##
##  This function tests if the module  <M>  is \tau-periodic, that is,
##  if  M \simeq \tau^i(M)  when  i  ranges over the set {1,2,...,n}.
##  Otherwise it returns false.
##
InstallMethod( IsTauPeriodic, 
   "for a path algebra matmodule and an integer",
   [ IsPathAlgebraMatModule, IS_INT  ], 0,
   function( M, n ) 

   local N0, N1, i;
 
   N0 := M;
   for i in [1..n] do
      N1 := DTr(N0);
      if IsomorphicModules(M,N1) then
         return i;
      else
         N0 := N1;
      fi;
   od;
   return false;
end
);

#######################################################################
##
#O  1stSyzygy( <M> )
##
##  This function computes the first syzygy of the module  <M>  by first
##  finding a minimal set of generators for  <M>, then finding the 
##  projective cover and computing the first syzygy as a submodule of 
##  this module.
##
InstallMethod( 1stSyzygy,
   "for a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M );

   if Dimension( M ) = 0 then 
     return ZeroModule( RightActingAlgebra( M ) );
   else
     return Kernel( ProjectiveCover( M ) );
   fi;
end
  );

#######################################################################
##
#O  NthSyzygy( <M>, <n> )
##
##  This function computes the <n>-th syzygy of the module <M>. 
##
InstallMethod( NthSyzygy,
   "for a path algebra module and a positive integer",
   [ IsPathAlgebraMatModule, IsInt ], 0,
   function( M, n ) 

  local projres, diff;

  if n < 0 then
     Error( "Entered integer is negative.\n" );
  fi;
  if n = 0 then
     return M;
  fi;
  projres := ProjectiveResolution( M ); 
  diff := DifferentialOfComplex( projres, n - 1 ); 
  
  return Kernel( diff );
end
  );

#######################################################################
##
#O  RightApproximationByAddM( <M>, <C> )
##
##  This function computes a right add<M>-approximation of the module  
##  <C>, and the approximation is not necessarily minimal.
##
InstallMethod ( RightApproximationByAddM, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( M, C )

   local K, HomMC, EndM, radEndM, radHomMC, i, j, FlatHomMC, 
         FlatradHomMC, V, BB, W, f, VoverW, B, gens, approx, approxmap;

   if RightActingAlgebra(M) <> RightActingAlgebra(C) then
       Error(" the two modules entered into MinimalRightApproximation are not modules over the same algebra.");
       return fail;
   fi;
   if Dimension(C) = 0 then
       return ZeroMapping(ZeroModule(RightActingAlgebra(M)),C);
   fi;
   K := LeftActingDomain(M);
   HomMC := HomOverAlgebra(M,C);
   if Length(HomMC) = 0 then 
       return ZeroMapping(ZeroModule(RightActingAlgebra(M)),C);
   else  
       EndM  := EndOverAlgebra(M);
       radEndM := RadicalOfAlgebra(EndM);
       radEndM := BasisVectors(Basis(radEndM));
       radEndM := List(radEndM, x -> FromEndMToHomMM(M,x));
       radHomMC := [];
       for i in [1..Length(HomMC)] do
           for j in [1..Length(radEndM)] do
               Add(radHomMC,radEndM[j]*HomMC[i]);
           od;
       od;
       FlatHomMC := List(HomMC, x -> Flat(x!.maps));
       FlatradHomMC := List(radHomMC, x -> Flat(x!.maps));
       V := VectorSpace(K,FlatHomMC,"basis");
       BB := Basis(V,FlatHomMC);
       W := Subspace(V,FlatradHomMC);
       f := NaturalHomomorphismBySubspace( V, W );
       VoverW := Range(f);
       B := BasisVectors(Basis(VoverW));
       gens := List(B, x -> PreImagesRepresentative(f,x)); 
       gens := List(gens, x -> Coefficients(BB,x));
       gens := List(gens, x -> LinearCombination(HomMC,x));         
       approx := List(gens, x -> Source(x));
       approx := DirectSumOfQPAModules(approx);
       approxmap := ShallowCopy(DirectSumProjections(approx));
       approxmap := List([1..Length(approxmap)], x -> approxmap[x]*gens[x]);         
       
       approxmap := Sum(approxmap);
       return approxmap;
   fi;
end
);

#######################################################################
##
#O  MinimalRightApproximation( <M>, <C> )
##
##  This function computes the minimal right add<M>-approximation of the
##  module  <C>.  TODO/CHECK: If one can modify the algorithm as indicated
##  below with ####.
##
InstallMethod ( MinimalRightApproximation, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( M, C )

    local   f;
   
   f := RightApproximationByAddM( M, C );
   
   return RightMinimalVersion( f )[ 1 ];
end
);

#######################################################################
##
#O   LeftApproximationByAddM( <C>, <M> )
##
##  This function computes a left add<M>-approximation of the module
##  <C>.
##
InstallMethod ( LeftApproximationByAddM, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( C, M )

    local   K,  HomCM,  EndM,  radEndM,  radHomCM,  i,  j,  FlatHomCM,  
            FlatradHomCM,  V,  BB,  W,  f,  VoverW,  B,  gens,  
            approx,  approxmap;

   if RightActingAlgebra(M) <> RightActingAlgebra(C) then
       Error(" the two modules entered into LeftApproximationByAddM are not modules over the same algebra.");
       return fail;
   fi;
   if Dimension(C) = 0 then
       return ZeroMapping(C, ZeroModule(RightActingAlgebra(M)));
   fi;
   K := LeftActingDomain(M);
   HomCM := HomOverAlgebra(C,M);
   if Length(HomCM) = 0 then 
       return ZeroMapping(C,ZeroModule(RightActingAlgebra(M)));
   else  
       EndM  := EndOverAlgebra(M);
       radEndM := RadicalOfAlgebra(EndM);
       radEndM := BasisVectors(Basis(radEndM));
       radEndM := List(radEndM, x -> FromEndMToHomMM(M,x));
       radHomCM := [];
       for i in [1..Length(HomCM)] do
           for j in [1..Length(radEndM)] do
               Add(radHomCM,HomCM[i]*radEndM[j]);
           od;
       od;
       FlatHomCM := List(HomCM, x -> Flat(x!.maps));
       FlatradHomCM := List(radHomCM, x -> Flat(x!.maps));
       V := VectorSpace(K,FlatHomCM,"basis");
       BB := Basis(V,FlatHomCM);
       W := Subspace(V,FlatradHomCM);
       f := NaturalHomomorphismBySubspace( V, W );
       VoverW := Range(f);
       B := BasisVectors(Basis(VoverW));
       gens := List(B, x -> PreImagesRepresentative(f,x)); 
       gens := List(gens, x -> Coefficients(BB,x));
       gens := List(gens, x -> LinearCombination(HomCM,x));
       approx := List(gens, x -> Range(x));
       approx := DirectSumOfQPAModules(approx);
       approxmap := ShallowCopy(DirectSumInclusions(approx));
       for i in [1..Length(approxmap)] do
           approxmap[i] := gens[i]*approxmap[i];
       od;
       approxmap := Sum(approxmap);
       
       return approxmap;
   fi;
end
  );



#######################################################################
##
#O   MinimalLeftApproximation( <C>, <M> )
##
##  This function computes the minimal left add<M>-approximation of the
##  module  <C>.  
##
InstallMethod ( MinimalLeftApproximation, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( C, M )

    local   f;
   
   f := LeftApproximationByAddM( C, M );
   
   return LeftMinimalVersion( f )[ 1 ];
end
  );


#######################################################################
##
#O  ProjectiveCover(<M>)
##
##  This function finds the projective cover P(M) of the module  <M>  
##  in that it returns the map from P(M) ---> M. 
##
InstallMethod ( ProjectiveCover, 
   "for a PathAlgebraMatModule",
   true,
   [ IsPathAlgebraMatModule ],
   0,
   function( M )

   local mingen, maps, PN, projections;

   if Dimension(M) = 0 then 
      return ZeroMapping(ZeroModule(RightActingAlgebra(M)), M);
   else 
      mingen := MinimalGeneratingSetOfModule(M);
      maps := List(mingen, x -> HomFromProjective(x,M));
      PN := List(maps, x -> Source(x));
      PN := DirectSumOfQPAModules(PN);
      projections := DirectSumProjections(PN);

      return projections*maps;;
   fi;
end
);

#######################################################################
##
#A  InjectiveEnvelope(< M >)
##
##  This function finds the injective envelope I(M) of the module  <M>  
##  in that it returns the map from M ---> I(M). 
##
InstallMethod ( InjectiveEnvelope, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule ],
    0,
    function( M );

    return DualOfModuleHomomorphism( ProjectiveCover( DualOfModule( M ) ) );
end
);

#######################################################################
##
#O  ExtOverAlgebra(<M>,<N>)
##
##  This function returns a list of three elements: (1) the kernel of 
##  the projective cover  Omega(<M>) --> P(M), (2) a basis of 
##  Ext^1(<M>,<N>)  inside  Hom(Omega(<M>),<N>)  and (3) a function 
##  that takes as an argument a homomorphism in  Hom(Omega(<M>),<N>)
##  and returns the coefficients of this element when written in 
##  terms of the basis of  Ext^1(<M>,<N>), if the group is non-zero. 
##  Otherwise it returns an empty list. 
##
InstallMethod( ExtOverAlgebra, 
   "for two PathAlgebraMatModule's",
   true, 
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
   function( M, N )

   local index, fromflatHomToHom, K, f, g, PM, syzygy, G, H, Img1, zero, 
         genssyzygyN, VsyzygyN, Img, gensImg, VImg, pi, ext, preimages, 
         homvecs, dimsyz, dimN, vec, t, l, i, H2, coefficients;
#
# Defining functions for later use.
#
    index := function( n )
   	if n = 0 then 
            return 1;
        else
            return n;
        fi;
    end;

    fromflatHomToHom := function( flat, dim, dim2 )
        local totalmat, matrix, i, start, j, d, d2;

        d := List(dim, index);
        d2 := List(dim2, index);
    	start := 0;
    	totalmat := [];
    	for i in [1..Length(dim)] do
            matrix := [];
            for j in [1..d[i]] do
            	Add(matrix, flat{[start+1..start+d2[i]]});
            	start := start+d2[i];
            od;
            Add(totalmat,matrix);
    	od;
        return totalmat;
    end;
#
# Test of input.
#
   if RightActingAlgebra(M) <> RightActingAlgebra(N) then 
      Error(" the two modules entered are not modules over the same algebra.\n"); 
   else  
      K := LeftActingDomain(M);
#
# creating a short exact sequence 0 -> Syz(M) -> P(M) -> M -> 0
# f: P(M) -> M, g: Syz(M) -> P(M)  
#
      f := ProjectiveCover(M);
      g := KernelInclusion(f);
      PM := Source(f);
      syzygy := Source(g);
#
# using Hom(-,N) on the s.e.s. above
#
      G := HomOverAlgebra(PM,N);
      H := HomOverAlgebra(syzygy,N);
#
# Making a vector space of Hom(Syz(M),N)
# by first rewriting the maps as vectors
#
      genssyzygyN := List(H, x -> Flat(x!.maps));
      if Length(genssyzygyN) = 0 then
         return [ g, [ ], [ ] ];
      else
         VsyzygyN := VectorSpace(K, genssyzygyN);
#
# finding a basis for im(g*)
# first, find a generating set of im(g*)
# 
         Img1 := g*G;
#
# removing 0 maps by comparing to zero = Zeromap(syzygy,N)
#
         zero := ZeroMapping(syzygy,N);
         Img  := Filtered(Img1, x -> x <> zero);
#
# Rewriting the maps as vectors
#
         gensImg := List(Img, x -> Flat(x!.maps));
#
# Making a vector space of <Im g*>
         VImg := Subspace(VsyzygyN, gensImg);  
#
# Making the vector space Ext1(M,N)
#
         pi := NaturalHomomorphismBySubspace(VsyzygyN, VImg);
         ext := Range(pi);
         if Dimension(Range(pi)) = 0 then 
            return [ g, [ ], [ ] ];
         else 
#
# Sending elements of ext back to Hom(Syz(M),N)
#
            preimages := List(BasisVectors(Basis(ext)), x -> PreImagesRepresentative(pi,x));
#
# need to put the parentheses back in place
#
            homvecs := []; # to store all lists of matrices, one list for each element (homomorphism)
            dimsyz := DimensionVector(syzygy);
            dimN := DimensionVector(N);
            homvecs := List(preimages, x -> fromflatHomToHom(x,dimsyz,dimN));
#
# Making homomorphisms of the elements
#
            H2 := List(homvecs, x -> RightModuleHomOverAlgebra(syzygy,N,x));

            coefficients := function( map ) 
               local vector, B;
       
               vector := ImageElm(pi,Flat(map!.maps)); 
               B := Basis(Range(pi));
               
               return Coefficients(B,vector);
            end;

            return [g,H2,coefficients];
         fi;
      fi;
   fi;
end
);

#######################################################################
##
#O  ExtAlgebraGenerators(<M>,<n>)
##
##  This function computes a set of generators of the Ext-algebra Ext^*(M,M)
##  up to degree  <n>. It returns a list of three elements, where the  
##  first element is the dimensions of Ext^[0..n](M,M), the second element
##  is the number of a set generators in the degrees [0..n], and the 
##  third element is the generators in these degrees. TODO: Create a 
##  minimal set of generators. Needs to take the radical of the degree
##  zero part into account.
##
InstallMethod( ExtAlgebraGenerators, 
    "for a module over a quotient of a path algebra",
    [ IsPathAlgebraMatModule, IS_INT ], 0,
    function( M, n ) 

    local N, projcovers, f, i, EndM, J, gens, extgroups, dim_ext_groups, 
          generators, productelements, j, induced, k, liftings, products, 
          l, m, productsinbasis, W, V, K, extalggenerators, p, tempgens, 
          I, g, templist, idealsquare; 

    K := LeftActingDomain(M);
    N := M;
    projcovers := [];
    for i in [1..n] do
        f := ProjectiveCover(N);
        Add(projcovers,f);
        N := Kernel(f);
    od;
    extgroups := [];
    #
    #   Computing Ext^i(M,M) for i = 0, 1, 2,...., n. 
    #
    for i in [0..n] do 
        if i = 0 then 
            EndM := EndOverAlgebra(M); 
            J := RadicalOfAlgebra(EndM);
            gens := GeneratorsOfAlgebra(J);
            Add(extgroups, [[],List(gens, x -> FromEndMToHomMM(M,x))]);
        elif i = 1 then 
            Add(extgroups, ExtOverAlgebra(M,M)); 
        else
            Add(extgroups, ExtOverAlgebra(Kernel(projcovers[i-1]),M)); 
        fi;
    od;
    dim_ext_groups := List(extgroups, x -> Length(x[2]));
#    Print(Dimension(EndM) - Dimension(J)," generators in degree 0.\n");    
    #
    #   Computing Ext^j(M,M) x Ext^(i-j)(M,M) for j = 1..i-1 and i = 2..n.
    #
    generators := List([1..n + 1], x -> 0);
    extalggenerators := List([1..n + 1], x -> []);
    for i in [1..n] do
        productelements := [];
        for j in [0..i] do
            if ( dim_ext_groups[i+1] <> 0 ) and ( dim_ext_groups[j+1] <> 0 ) and ( dim_ext_groups[i-j+1] <> 0 ) then 
                induced := ShallowCopy(extgroups[i-j+1][2]);
                for k in [1..j] do
                    liftings := List([1..dim_ext_groups[i-j+1]], x -> projcovers[i-j+k]*induced[x]);
                    liftings := List([1..dim_ext_groups[i-j+1]], x -> LiftingMorphismFromProjective(projcovers[k],liftings[x]));
                    induced  := List([1..dim_ext_groups[i-j+1]], x -> MorphismOnKernel(projcovers[i-j+k],projcovers[k],liftings[x],induced[x]));
                od;
                products := [];
                for l in [1..dim_ext_groups[j+1]] do
                    for m in [1..dim_ext_groups[i-j+1]] do 
                        Add(products,induced[m]*extgroups[j+1][2][l]);
                    od;
                od;
                productsinbasis := List(products, x -> extgroups[i+1][3](x));
                productsinbasis := Filtered(productsinbasis, x -> x <> Zero(x));
                if Length(productsinbasis) <> 0 then
                    Append(productelements,productsinbasis);
                fi;
            fi;
        od;
        if Length(productelements) <> 0 then 
            W := FullRowSpace(K,Length(extgroups[i+1][2]));            
            V := Subspace(W,productelements);
            if dim_ext_groups[i+1] > Dimension(V) then
                generators[i+1] := Length(extgroups[i+1][2]) - Dimension(V); 
                p := NaturalHomomorphismBySubspace(W,V);
                tempgens := List(BasisVectors(Basis(Range(p))), x -> PreImagesRepresentative(p,x));
                extalggenerators[i+1] := List(tempgens, x -> LinearCombination(extgroups[i+1][2],x));                 
#                Print(Length(extgroups[i+1][2]) - Dimension(V)," new generator(s) in degree ",i,".\n");
            fi;
        elif dim_ext_groups[i+1] <> 0 then 
            generators[i+1] := Length(extgroups[i+1][2]); 
            extalggenerators[i+1] := extgroups[i+1][2];                 
#            Print(Length(extgroups[i+1][2])," new generator(s) in degree ",i,".\n");
        fi;
    od; 
    dim_ext_groups[1] := Dimension(EndM);
    templist := [];
    for i in [1..Length(gens)] do
        for j in [1..Dimension(EndM)] do
            Add(templist,BasisVectors(Basis(EndM))[j]*gens[i]);
        od;
    od;
    templist := Filtered(templist, x -> x <> Zero(x));
    idealsquare := [];
    for i in [1..Length(templist)] do
        for j in [1..Length(gens)] do
            Add(idealsquare,gens[j]*templist[i]);
        od;
    od;
    I := Ideal(EndM,idealsquare);
    g := NaturalHomomorphismByIdeal(EndM,I);
    extalggenerators[1] := List(BasisVectors(Basis(Range(g))), x -> FromEndMToHomMM(M,PreImagesRepresentative(g,x)));
    generators[1] := Length(extalggenerators[1]);
    return [dim_ext_groups,generators,extalggenerators];
end
);

#######################################################################
##
#O  PartialIyamaGenerator( <M> )
##
##  Given a module  <M>  this function returns the submodule of  <M>
##  given by the radical of the endomorphism ring of  <M>  times <M>.
##  If  <M>  is zero, then  <M>  is returned.
##
InstallMethod( PartialIyamaGenerator,
    "for a path algebra",
    [ IsPathAlgebraMatModule ], 0,
    function( M ) 

    local B, EndM, radEndM, Brad, subgens, b, r; 

    if Dimension(M) = 0 then 
        return M;
    fi;
    B := CanonicalBasis(M); 
    EndM := EndOverAlgebra(M);
    radEndM := RadicalOfAlgebra(EndM);  
    Brad := BasisVectors(Basis(radEndM));
    Brad := List(Brad, x -> FromEndMToHomMM(M,x));
    subgens := [];
    for b in B do
        for r in Brad do
            Add(subgens, ImageElm(r,b));
        od;
    od;

    return SubRepresentation(M,subgens);
end
);

#######################################################################
##
#O  IyamaGenerator( <M> )
##
##  Given a module  <M> this function returns a module  N  such that
##  <M>  is a direct summand of  N  and such that the global dimension
##  of the endomorphism ring of  N  is finite. 
##
InstallMethod( IyamaGenerator,
    "for a path algebra",
    [ IsPathAlgebraMatModule ], 0,
    function( M ) 

    local iyamagen, N, IG, L; 

    iyamagen := [];
    N := M;
    repeat 
        Add(iyamagen,N);
        N := PartialIyamaGenerator(N);
    until
        Dimension(N) = 0;

    IG := DirectSumOfQPAModules(iyamagen);
    if IsFinite(LeftActingDomain(M)) then 
        L := DecomposeModuleWithMultiplicities(IG);
        IG := DirectSumOfQPAModules(L[1]);
    fi;

    return IG;
end
);

#######################################################################
##
#O  GlobalDimensionOfAlgebra( <A>, <n> )
##
##  Returns the global dimension of the algebra  <A>  if it is less or
##  equal to  <n>, otherwise it returns false.
##  
InstallMethod ( GlobalDimensionOfAlgebra, 
    "for a finite dimensional quotient of a path algebra",
    true,
    [ IsQuiverAlgebra, IS_INT ], 
    0,
    function( A, n )

    local simples, S, projres, dimension, j;
    
    if not IsFiniteDimensional(A) then
        TryNextMethod();
    fi;
    if HasGlobalDimension(A) then
        return GlobalDimension(A);
    fi;
    if Trace(AdjacencyMatrixOfQuiver(QuiverOfPathAlgebra(A))) > 0 then
        SetGlobalDimension(A,infinity);
        return infinity;
    fi;
    simples := SimpleModules(A);
    dimension := 0;
    for S in simples do
        projres := ProjectiveResolution(S);
        j := 0;
        while ( j < n + 1 ) and ( Dimension(ObjectOfComplex(projres,j)) <> 0 ) do 
            j := j + 1;
        od; 
        if ( j < n + 1 ) then 
            dimension := Maximum(dimension, j - 1 );
        else
            if ( Dimension(ObjectOfComplex(projres,n + 1)) <> 0 ) then
                return false;
            else
                dimension := Maximum(dimension, n );
            fi;
        fi;
    od;
    
    SetGlobalDimension(A,dimension);
    return dimension;
end 
);

#######################################################################
##
#O  DominantDimensionOfModule( <M>, <n> )
##
##  Returns the dominant dimension of the module  <M>  if it is less or
##  equal to  <n>. If the module  <M>  is injective and projective, then 
##  it returns infinity. Otherwise it returns false.
##  
InstallMethod ( DominantDimensionOfModule, 
    "for a finite dimensional quotient of a path algebra",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )

    local A, Mop, resop, dimension, j;
    
    A := RightActingAlgebra(M);
    #
    # If the algebra  A  is not a finite dimensional, try another method.
    #
    if not IsFiniteDimensional(A) then
        TryNextMethod();
    fi;
    # 
    # If the module  <M>  is injective and projective, then the dominant dimension of 
    # <M>  is infinite.
    if IsInjectiveModule( M ) and IsProjectiveModule( M ) then 
        return infinity;
    fi;
    
    Mop := DualOfModule( M );
    #
    # Check how far out the modules in the minimal projective resolution of Mop
    # is injective.
    #
    resop := ProjectiveResolution( Mop );
    dimension := 0;
    j := 0;
    while ( j < n ) and ( IsInjectiveModule( ObjectOfComplex( resop, j ) ) ) do 
        j := j + 1;
    od; 
    if ( j < n ) then # if this happens, then j-th projective is not injective and domdim equal to j for <M>.
        dimension := Maximum( dimension, j );
    else                    # if this happens, then j = n. 
        if not IsInjectiveModule( ObjectOfComplex( resop, n  ) ) then # then domdim is  n  of  <M>.
            dimension := Maximum( dimension, n );
        else                                                      # then domdim is > n  for  <M>.
            return false;
        fi;
    fi;

#    SetDominantDimension(M,dimension);
    return dimension;
end 
);


#######################################################################
##
#O  DominantDimensionOfAlgebra( <A>, <n> )
##
##  Returns the dominant dimension of the algebra  <A>  if it is less or
##  equal to  <n>. If the algebra  <A>  is selfinjectiv, then it returns
##  infinity. Otherwise it returns false.
##  
InstallMethod ( DominantDimensionOfAlgebra, 
    "for a finite dimensional quotient of a path algebra",
    true,
    [ IsQuiverAlgebra, IS_INT ], 
    0,
    function( A, n )

    local   P,  domdimlist,  pos_infinite,  pos_false,  pos_finite;
    #
    # If the algebra  <A>  is not a finite dimensional, try another method.
    #
    if not IsFiniteDimensional( A ) then
        TryNextMethod( );
    fi;
    # 
    # If the algebra  <A>  is selfinjective, then the dominant dimension of 
    # <A>  is infinite.
    if IsSelfinjectiveAlgebra( A ) then 
        return infinity;
    fi;
    P := IndecProjectiveModules( A );
    domdimlist := List( P, p -> DominantDimensionOfModule( p, n ) );   # Checking the dominant dimension of each indec. projective module. 
    pos_infinite := Positions( domdimlist, infinity ); # positions of the indec. projective modules with infinite domdim.
    pos_false := Positions( domdimlist, false ); # positions of the indec. projective modules with bigger dom.dim. than n.
    pos_finite := [ 1..Length( P ) ];     
    SubtractSet( pos_finite, Union( pos_infinite, pos_false ) ); # positions of the indec. projective modules with domdim less or equal to n.
    
    if Length( pos_finite ) = 0 then 
        return false;
    else
        return Minimum( domdimlist{ pos_finite } );
    fi;
end
  );

#######################################################################
##
#O  ProjDimensionOfModule( <M>, <n> )
##
##  Returns the projective dimension of the module  <M>  if it is less
##  or equal to  <n>, otherwise it returns false.
##  
InstallMethod ( ProjDimensionOfModule, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )

    local projres, i;
    
    if HasProjDimension(M) then
        if ProjDimension(M) > n then
            return false;
        else
            return ProjDimension(M);
        fi;
    fi;
    projres := ProjectiveResolution(M);
    i := 1;
    while i < n + 2 do
        if Dimension(ObjectOfComplex(projres,i)) = 0 then # the projective dimension is  i - 1.
            SetProjDimension(M, i - 1);
            return i-1;
        fi;
        i := i + 1;
    od; 
    
    return false;
end 
  );

#######################################################################
##
#O  GorensteinDimensionOfAlgebra( <A>, <n> )
##
##  Returns the Gorenstein dimension of the algebra  <A>  if it is less
##  or equal to  <n>, otherwise it returns false.
##  
InstallMethod ( GorensteinDimensionOfAlgebra,
    "for a quiver algebra",
    true,
    [ IsQuiverAlgebra, IS_INT ], 
    0,
    function( A, n )

    local Ilist, dimension_right, I, projdim, Ilistop, dimension_left;
    
    #
    #  If  <A>  has finite global dimension, then the Gorenstein 
    #  dimension of  <A>  is equal to the global dimension of  <A>. 
    #
    if HasGlobalDimension(A) then
        if GlobalDimension(A) <> infinity then 
            SetGorensteinDimension(A, GlobalDimension(A));
            return GlobalDimension(A);
        fi;
    fi;
    if HasGorensteinDimension(A) then
        return GorensteinDimension(A);
    fi;
    #
    #  First checking the injective dimension of the indecomposable
    #  projective left  <A>-modules.
    #
    Ilist := IndecInjectiveModules(A);
    dimension_right := 0;
    for I in Ilist do
        projdim := ProjDimensionOfModule(I,n); 
        if projdim = false then   # projective dimension of  I  is bigger than  n.
            return false;
        else                      # projective dimension of  I  is less or equal to  n.
            dimension_right := Maximum(dimension_right, projdim );
        fi;
    od;
    #
    #  Secondly checking the injective dimension of the indecomposable
    #  projective right  <A>-modules.
    #  
    Ilistop := IndecInjectiveModules(OppositeAlgebra(A));    
    dimension_left := 0;
    for I in Ilistop do
        projdim := ProjDimensionOfModule(I,n); 
        if projdim = false then   # projective dimension of  I  is bigger than  n.
            return false;
        else                      # projective dimension of  I  is less or equal to  n.
            dimension_left := Maximum(dimension_left, projdim );
        fi;
    od;    
    if dimension_left <> dimension_right then
        Print("You have a counterexample to the Gorenstein conjecture!\n\n");
    fi;
    
    SetGorensteinDimension(A, Maximum(dimension_left, dimension_right));
    return Maximum(dimension_left, dimension_right);    
end 
  );

#######################################################################
##
#O  N_RigidModule( <M>, <n> )
##
##  Returns true if the module  <M>  is n-rigid, otherwise false.
##  
InstallMethod ( N_RigidModule, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )

    local i, N;
    
    i := 1;
    N := M;
    while i < n + 1 do
        if Length(ExtOverAlgebra(N,M)[2]) > 0 then # Ext^i(M,M) <> 0.
            return false;
        else                                       # Ext^i(M,M) = 0.
            i := i + 1;                            
            N := NthSyzygy(N,1);
        fi;
    od;
    
    return true;
end 
  );

#######################################################################
##
#O  LeftFacMApproximation( <C>, <M> )
##
##  This function computes a left, not necessarily left minimal, 
##  Fac<M>-approximation of the module  <C>  by doing the following:
##                  p
##           P(C) -------> C  (projective cover)
##             |           |
##           f |           | g - the returned homomorphism
##             V           V
##          M^{P(C)} ----> E
##
InstallMethod( LeftFacMApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( C, M )

    local p, f; 
    #
    # Checking if the modules  <C>  and  <M>  are modules over the same algebra.
    #
    if RightActingAlgebra(C) <> RightActingAlgebra(M) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;    
    p := ProjectiveCover(C);
    f := MinimalLeftAddMApproximation(Source(p),M);
    
    return PushOut(p,f)[2];
end
  );


#######################################################################
##
#O  MinimalLeftFacMApproximation( <C>, <M> )
##
##  This function computes a minimal left Fac<M>-approximation of the 
##  module  <C>. 
##
InstallMethod( MinimalLeftFacMApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( C, M );

    return LeftMinimalVersion(LeftFacMApproximation(C,M))[1];
end
  );


#######################################################################
##
#O  RightSubMApproximation( <M>, <C> )
##
##  This function computes a right, not necessarily a right minimal,
##  Sub<M>-approximation of the module  <C>  by doing the following:
##
##             E -----> M^{I(C)}  (minimal right Add<M>-approximation)
##             |           |
##           g |           | f     g - the returned homomorphism
##             V    i      V
##             C -------> I(C)  (injective envelope)
##
InstallMethod( RightSubMApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( M, C )

    local iop, i, f; 
    #
    # Checking if the modules  <M>  and  <C>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(C) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;      
    iop := ProjectiveCover(DualOfModule(C));
    i := DualOfModuleHomomorphism(iop);
    f := MinimalRightAddMApproximation(M, Range(i));
    
    return PullBack(i,f)[2];
end
  );


#######################################################################
##
#O  MinimalRightSubMApproximation( <M>, <C> )
##
##  This function computes a minimal right Sub<M>-approximation of the module 
##  <C>. 
##
InstallMethod( MinimalRightSubMApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( M, C );

    return RightMinimalVersion(RightSubMApproximation(M,C))[1];
end
  );


#######################################################################
##
#O  LeftSubMApproximation( <C>, <M> )
##
##  This function computes a minimal left Sub<M>-approximation of the 
##  module  <C>  by doing the following: 
##            f
##       C -------> M^{C} = the minimal left Add<M>-approxiamtion
##
##   Returns the natural projection  C --------> Im(f). 
##
InstallMethod( LeftSubMApproximation,
    "for a representations of a quiver",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
    function( C, M )

    local f; 
    #
    # Checking if the modules  <C>  and  <M>  are modules over the same algebra.
    #
    if RightActingAlgebra(C) <> RightActingAlgebra(M) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;      
    f := MinimalLeftAddMApproximation(C,M);
    
    return ImageProjection(f);
end
  );

#######################################################################
##
#O  InjDimensionOfModule( <M>, <n> )
##
##  Returns the injective dimension of the module  <M>  if it is less
##  or equal to  <n>, otherwise it returns false.
##  
InstallMethod ( InjDimensionOfModule, 
    "for a PathAlgebraMatModule and a positive integer",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )
    
    local DM, projresop, i;
    
    if HasInjDimension(M) then 
        if InjDimension(M) > n then
            return false;
        else 
            return InjDimension(M);
        fi;
    fi;
    DM := DualOfModule(M);     
    if HasProjDimension(DM) then
        return ProjDimension(DM);
    fi;
    projresop := ProjectiveResolution(DM);
    i := 1;
    while i < n + 2 do
        if Dimension(ObjectOfComplex(projresop,i)) = 0 then # the projective dimension of DM is  i - 1.
            SetInjDimension(M, i - 1);
            return i-1;
        fi;
        i := i + 1;
    od; 
    
    return false;
end 
  );

#######################################################################
##
#O  HaveFiniteCoresolutionInAddM( <N>, <M>, <n> )
##
##  This function checks if the module  <N>  has a finite coresolution
##  in  add<M>  of length at most  <n>.  If it does, then this 
##  coresoultion is returned, otherwise false is returned. 
##
InstallMethod( HaveFiniteCoresolutionInAddM,
    "for two PathAlgebraMatModules and a positive integer",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule, IS_INT ],
    
    function( N, M, n )
    local U, differentials, g, f, i, cat, coresolution;
    #
    # Checking if  <N>  and  <M>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;
    #
    # If n = 0, then this is the same as N being in add M.
    # 
    # Computing successive minimal left add<M>-approximation to produce
    # a coresolution of  <N>  in  add<M>.
    #
    cat := CatOfRightAlgebraModules(RightActingAlgebra(M));
    U := N;
    differentials := [];
    g := IdentityMapping(U);
    f := MinimalLeftAddMApproximation(U,M);
    if n = 0 then 
      if IsIsomorphism( f ) then 
        return FiniteComplex(cat, 0, [ f ] );
      else
        return false;
      fi;
    fi;
    for i in [0..n] do
        if not IsInjective(f) then 
            return false;
        fi;
        Add(differentials, g*f); 
        g := CoKernelProjection(f);
        if Dimension(Range(g)) = 0 then
            break;
        fi;
        f := MinimalLeftAddMApproximation(Range(g),M); 
    od;
    differentials := Reversed(differentials);
    coresolution := FiniteComplex(cat, -Length(differentials) + 1, differentials);
    return coresolution;
end
  );


#######################################################################
##
#O  TiltingModule( <M>, <n> )
##
##  This function checks if the module  <M>  is a tilting module of 
##  projective dimension at most  <n>.
##
InstallMethod( TiltingModule,
    "for a PathAlgebraMatModule and a positive integer",
    [ IsPathAlgebraMatModule, IS_INT ],
    
    function( M, n )
    local m, N, i, P, t, coresolution, temp;
    #
    # Checking if the module has projective dimension at most  <n>.
    #
    if HasProjDimension(M) then 
        if ProjDimension(M) > n then
            return false;
        fi;
    else
        if not IS_INT(ProjDimensionOfModule(M,n)) then
            return false;
        fi;
    fi;
    #
    # Now we know that the projective dimension of  <M>  is at most  <n>.
    # Next we check if the module  <M>  is selforthogonal. 
    # 
    m := ProjDimension(M); 
    N := M;
    i := 0;
    repeat 
        if Length(ExtOverAlgebra(N,M)[2]) <> 0 then
            return false;
        fi;
        i := i + 1;
        if i < m then 
            N := NthSyzygy(N,1); 
        fi;
    until i = m + 1;
    #
    # Now we know that the projective dimension of  <M>  is at most  <n>, 
    # and that the module  <M>  is selforthogonal.  Next we check if all the 
    # indecomposable projectives can be coresolved in  add<M>.
    #
    P := IndecProjectiveModules(RightActingAlgebra(M));
    t := Length(P);
    coresolution := [];
    for i in [1..t] do
        temp := HaveFiniteCoresolutionInAddM(P[i], M, m); 
        if temp = false then
            return false;
        fi;
        Add(coresolution, temp);
    od;
    #
    # Now we know that the module  <M>  is a tilting module of projective 
    # dimension  m. 
    # 
    SetIsTiltingModule(M, true);
    return [m, coresolution];
end
  );

#######################################################################
##
#O  HaveFiniteResolutionInAddM( <N>, <M>, <n> )
##
##  This function checks if the module  <N>  has a finite resolution
##  in  add<M>  of length at most  <n>.  If it does, then this 
##  resolution is returned, otherwise false is returned. 
##
InstallMethod( HaveFiniteResolutionInAddM,
    "for two PathAlgebraMatModules and a positive integer",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule, IS_INT ],
    
    function( N, M, n )
    local U, differentials, g, f, i, cat, resolution;
    #
    # Checking if  <N>  and  <M>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(N) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;
    #
    # If n = 0, then this is the same as N being in add M.
    # 
    # Computing successive minimal right add<M>-approximation to produce
    # a resolution of  <N>  in  add<M>.
    #
    U := N;
    differentials := [];
    g := IdentityMapping(U);
    f := MinimalRightAddMApproximation(M,U);
    if n = 0 then 
      if IsIsomorphism( f ) then 
        return true;
      else
        return false;
      fi;
    fi;
    for i in [0..n] do
        if not IsSurjective(f) then
            return false;
        fi;
        Add(differentials, f*g);
        g := KernelInclusion(f);
        if Dimension(Source(g)) = 0 then
            break;
        fi;
        f := MinimalRightAddMApproximation(M, Source(g));
    od;
    cat := CatOfRightAlgebraModules(RightActingAlgebra(M));
    resolution := FiniteComplex(cat, 1, differentials);
    return resolution;
end
  );


#######################################################################
##
#O  CotiltingModule( <M>, <n> )
##
##  This function checks if the module  <M>  is a cotilting module of 
##  projective dimension at most  <n>.
##
InstallMethod( CotiltingModule,
    "for a PathAlgebraMatModule and a positive integer",
    [ IsPathAlgebraMatModule, IS_INT ],
    
    function( M, n )
    local m, N, i, I, t, resolution, temp;
    #
    # Checking if the module has injective dimension at most  <n>.
    #
    if HasInjDimension(M) then 
        if InjDimension(M) > n then
            return false;
        fi;
    else
        if not IS_INT(InjDimensionOfModule(M,n)) then
            return false;
        fi;
    fi;
    #
    # Now we know that the injective dimension of  <M>  is at most  <n>.
    # Next we check if the module  <M>  is selforthogonal. 
    # 
    m := InjDimension(M); 
    N := M;
    i := 0;
    repeat 
        if Length(ExtOverAlgebra(N,M)[2]) <> 0 then
	    SetIsCotiltingModule(M, false);
            return false;
        fi;
        i := i + 1;
        if i < m then 
            N := NthSyzygy(N,1); 
        fi;
    until i = m + 1;
    #
    # Now we know that the injective dimension of  <M>  is at most  <n>, 
    # and that the module  <M>  is selforthogonal.  Next we check if all the 
    # indecomposable injectives can be resolved in  add<M>.
    #
    I := IndecInjectiveModules(RightActingAlgebra(M));
    t := Length(I);
    resolution := [];
    for i in [1..t] do
        temp := HaveFiniteResolutionInAddM(I[i], M, m); 
        if temp = false then
            return false;
        fi;
        Add(resolution, temp);
    od;
    #
    # Now we know that the module  <M>  is a cotilting module of injective 
    # dimension  <M>. 
    # 
    SetIsCotiltingModule(M, true);    
    return [m, resolution];
end
  );


#######################################################################
##
#O  AllComplementsOfAlmostCompleteTiltingModule( <M>, <X> )
##
##  This function constructs all complements of an almost complete 
##  tilting module  <M>  given a complement  <X>  of  <M>.  The 
##  complements are returned as a long exact sequence (whenever possible)
##  of minimal left and minimal right  add<M>-approximations.
##
InstallMethod( AllComplementsOfAlmostCompleteTiltingModule, 
    "for two PathAlgebraMatModules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    
    function( M, X)
    local U, leftdifferentials, g, f, cat, resolution, 
          rightdifferentials, coresolution;
    #
    # Checking if  <M>  and  <X>  are modules over the same algebra.
    #
    if RightActingAlgebra(M) <> RightActingAlgebra(X) then
        Error("the entered modules are not modules over the same algebra,\n");
    fi;
    # 
    # Computing successive minimal right add<M>-approximation to produce
    # a resolution of  <X>  in  add<M>.
    #
    U := X;
    leftdifferentials := [];
    g := IdentityMapping(U);
    f := MinimalRightAddMApproximation(M,U);
    while IsSurjective(f) do
        Add(leftdifferentials, f*g);
        g := KernelInclusion(f);
        if Dimension(Source(g)) = 0 then
            break;
        fi;
        f := MinimalRightAddMApproximation(M, Source(g));
    od;
    cat := CatOfRightAlgebraModules(RightActingAlgebra(M));
    if Length(leftdifferentials) = 0 then
        resolution := [];
    else
        Add(leftdifferentials, KernelInclusion(leftdifferentials[Length(leftdifferentials)]));         
        resolution := FiniteComplex(cat, 1, leftdifferentials);
    fi;
    # 
    # Computing successive minimal left add<M>-approximation to produce
    # a coresolution of  <X>  in  add<M>.
    #
    U := X;
    rightdifferentials := [];
    g := IdentityMapping(U);
    f := MinimalLeftAddMApproximation(U,M);    
    while IsInjective(f) do
        Add(rightdifferentials, g*f); 
        g := CoKernelProjection(f);
        f := MinimalLeftAddMApproximation(Range(g),M); 
    od;

    if Length(rightdifferentials) = 0 then
        coresolution := [];
    else
        Add(rightdifferentials, CoKernelProjection(rightdifferentials[Length(rightdifferentials)])); 
        rightdifferentials := Reversed(rightdifferentials);
        coresolution := FiniteComplex(cat, -Length(rightdifferentials) + 1, rightdifferentials);    
    fi;
    
    return [resolution,coresolution];
end
  );

#######################################################################
##
#A  FaithfulDimension( <M> )
##
##  This function computes the faithful dimension of a module  <M>.
##
InstallMethod( FaithfulDimension,
    "for a PathAlgebraMatModules",
    [ IsPathAlgebraMatModule ],
    
    function( M )
    local P, lengths, p, tempdim, U, g, f;
    
    P := IndecProjectiveModules(RightActingAlgebra(M));
    # 
    # For all indecomposable projective  A-modules computing successive 
    # minimal left add<M>-approximation to produce find the faithful
    # dimension  <M>  has with respect to that indecomposable 
    # projective. 
    #
    lengths := [];
    for p in P do
        tempdim := 0; 
        U := p;
        f := MinimalLeftAddMApproximation( U, M );
	if not IsInjective( f ) then
	   return 0;
	fi;
        while IsInjective( f ) do 
            tempdim := tempdim + 1;
            g := CoKernelProjection( f );
            if Dimension( Range( g ) ) = 0 then
                Add( lengths, infinity );
                break;
            else
                f := MinimalLeftAddMApproximation( Range( g ), M ); 
            fi;
        od;
        if Dimension( Range( g ) ) <> 0 then 
            Add( lengths, tempdim );
        fi;
    od; 

    return Minimum( lengths );
end
  );

#######################################################################
##
#O  NumberOfComplementsOfAlmostCompleteTiltingModule( <M> )
##
##  This function computes the number complements of an almost 
##  complete tilting/cotilting module  <M>, assuming that  <M>
##  is an almost complete tilting module.
##
InstallMethod( NumberOfComplementsOfAlmostCompleteTiltingModule,
    "for a PathAlgebraMatModules",
    [ IsPathAlgebraMatModule ],
    
    function( M );
    
    return FaithfulDimension(M) + 1; 
end
  );

#######################################################################
##
#O  LeftMutationOfTiltingModuleComplement( <M>, <N> )
##
##  This function computes the left mutation of a complement  <N>  of
##  an almost complete tilting module  <M>, assuming that  <M>  is an 
##  almost complete tilting module.  If it doesn't exist, then the
##  function returns false.
##
InstallMethod( LeftMutationOfTiltingModuleComplement,
    "for a PathAlgebraMatModules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
        
    function( M, N )
    
    local f;
    
    f := MinimalLeftAddMApproximation(N, M);
    if IsInjective(f) then
        return CoKernel(f);
    else
        return false;
    fi;
end
  );

#######################################################################
##
#O  RightMutationOfTiltingModuleComplement( <M>, <N> )
##
##  This function computes the right mutation of a complement  <N>  of
##  an almost complete tilting module  <M>, assuming that  <M>  is an 
##  almost complete tilting module.  If it doesn't exist, then the
##  function returns false.
##
InstallMethod( RightMutationOfTiltingModuleComplement,
    "for a PathAlgebraMatModules",
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
        
    function( M, N )
    
    local f;
    
    f := MinimalRightAddMApproximation(M, N);
    if IsSurjective(f) then
        return Kernel(f);
    else
        return false;
    fi;
end
  );

##########################################################################
##
#P DeclareProperty( "IsHereditaryAlgebra", [ A ] )
##
## The function is defined for an admissible quotient  <A>  of a path
## algebra, and it returns true if  <A>  is a hereditary algebra and 
## false otherwise.
## 
InstallMethod( IsHereditaryAlgebra,
    "for an algebra",
    [ IsAdmissibleQuotientOfPathAlgebra ],
        
    function( A )
    local test;
    
    if HasGlobalDimension(A) then
        return GlobalDimension(A) < 2;
    fi;
    
    test := GlobalDimensionOfAlgebra(A, 1);
    #
    # test = false => gldim(A) > 1
    # test = infinity => gldim(A) = infinity
    # test = integer > 1 => gldim(A) > 1
    # if all these are not true, then gldim(A) < 1, and  A  is hereditary.
    #
    if IsBool(test) or ( test = infinity ) or ( IS_INT(test) and test > 1 ) then
        return false;
    else
        return true;
    fi;
end 
  );

#######################################################################
##
#O  RightApproximationByPerpT( <T>, <M> )
##
##  Returns the minimal right $\widehat{\add T}$-approximation of the 
##  module  <M>.  It checks if  <T>  is a cotilting module, and if not
##  it returns an error message. 
## 
InstallMethod ( RightApproximationByPerpT, 
    "for a cotilting PathAlgebraMatModule and a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    0,
    function( T, M )

    local   n,  projres,  exactsequences,  f,  i,  beta,  g,  h,  
            alpha,  fprime;
    
    if not IsCotiltingModule( T ) then
        Error("the first argument is not a cotilting module,\n");
    fi;
    if IsZero( M ) then
      return ZeroMapping( ZeroModule( RightActingAlgebra( M ) ), M );
    fi;
    n := InjDimension( T );
    projres := ProjectiveResolution( M );
    ObjectOfComplex( projres, n);
    exactsequences := List( [ 0..n - 1 ], i -> [ KernelInclusion(DifferentialOfComplex(projres, i)), 
                              ImageProjection(DifferentialOfComplex(projres, i)) ]); 
    f := MinimalLeftApproximation( Source( exactsequences[ n ][ 1 ] ), T );    
    for i in [ 1..n ] do
      beta := IsomorphismOfModules( Source( exactsequences[ n + 1 - i ][ 1 ] ), Source( f ) );
      f := beta * f;
      g := PushOut( exactsequences[ n + 1 - i ][ 1 ], f )[ 1 ];
      h := CoKernelProjection( g );
      alpha := IsomorphismOfModules( Range( h ), Range( exactsequences[ n + 1 - i ][ 2 ] ) );
      h := h * alpha;
      if i < n then
        fprime := MinimalLeftApproximation( Source( h ), T );
        f := PushOut( fprime, h )[ 1 ]; 
      fi;
    od;
    
    return RightMinimalVersion( h )[ 1 ];
end
);


#######################################################################
##
#O  LeftApproximationByAddTHat( <T>, <M> )
##
##  Returns the minimal left $\widehat{\add T}$-approximation of the 
##  module  <M>.  It checks if  <T>  is a cotilting module, and if not
##  it returns an error message. 
## 
InstallMethod ( LeftApproximationByAddTHat, 
    "for a cotilting PathAlgebraMatModule and a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
    0,
    function( T, M )

    local   n,  projres,  exactsequences,  currentsequences,  f,  i,  g,  
            h,  alpha,  t, fprime;
    
    if not IsCotiltingModule( T ) then
        Error("the first argument is not a cotilting module,\n");
    fi;
    if IsZero( M ) then
      return ZeroMapping( M, ZeroModule( RightActingAlgebra( M ) ) );
    fi;
    n := InjDimension( T );
    projres := ProjectiveResolution( M );
    ObjectOfComplex( projres, n);
    exactsequences := List( [ 0..n - 1 ], i -> [ KernelInclusion(DifferentialOfComplex(projres, i)), 
                              ImageProjection(DifferentialOfComplex(projres, i)) ]); 
    f := MinimalLeftApproximation( Source( exactsequences[ n ][ 1 ] ), T );    
    for i in [ 1..n ] do
      g := PushOut( exactsequences[ n + 1 - i ][ 1 ], f )[ 1 ];
      h := CoKernelProjection( g );
      alpha := IsomorphismOfModules( Range( h ), Range( exactsequences[ n + 1 - i ][ 2 ] ) );
      h := h * alpha;
      fprime := MinimalLeftApproximation( Source( h ), T );
      f := PushOut( fprime, h )[ 1 ]; 
    od;
        
    return LeftMinimalVersion( f )[ 1 ];
end
  );

InstallMethod( IsNthSyzygy, 
    "for a PathAlgebraMatModule",
    true,
    [ IsPathAlgebraMatModule, IS_INT ], 
    0,
    function( M, n )

    local N;
    
    N := NthSyzygy( DualOfModule( NthSyzygy( DualOfModule( M ), n ) ), n );
    N := DirectSumOfQPAModules( [ N, Source( ProjectiveCover( M ) ) ] );
    
    return IsDirectSummand( M, N ); 
end
  );

#######################################################################
##
#O  RightApproximationByAddM( <L>, <C> )
##
##  Given a list of module <L> over a finite dimensional quotient A of
##  a path algebra and a module <C> over  A, this function computes
##  a right approximation of <C> in the additive closure of the modules
##  in the list <L>.
##
InstallMethod ( RightApproximationByAddM,
    "for a list of modules and one module",
    true,
    [ IsList, IsPathAlgebraMatModule ],
    0,
    function( L, C )

    local   A,  K,  approximation,  i,  homL_iC,  homL_ilessthanL_i,  
            homL_iapproxC,  endL_i,  radendL_i,  RadendL_i,  radmaps,  
            generators,  radgenerators,  generatorshomL_iC,  g,  V,  
            t,  M,  projections,  f;

    A := RightActingAlgebra( C );
    if Length( L ) = 0 then
      return ZeroMapping( ZeroModule( A ), C );
    fi;
    if not ForAll( L, l -> RightActingAlgebra( l ) = A ) then
      Error( "Not all modules in the list of modules entered are modules over the same algebra.\n" );
    fi;
    
    K := LeftActingDomain( A );
    approximation := [ ];
    for i in [ 1..Length( L ) ] do
      homL_iC := HomOverAlgebra( L[ i ], C );
      if Length( homL_iC ) > 0 then
        homL_ilessthanL_i := List( [ 1..Length(approximation) ], r -> HomOverAlgebra( L[ i ], Source( approximation[ r ] ) ) );
        homL_iapproxC := Flat( List( [ 1..Length(approximation) ], r -> homL_ilessthanL_i[ r ] * approximation[ r ] ) );
        homL_iapproxC := Filtered( homL_iapproxC, h -> not IsZero( h ) );
        endL_i := EndOverAlgebra( L[ i ] );
        radendL_i := RadicalOfAlgebra( endL_i );
        RadendL_i := List( BasisVectors( Basis( radendL_i ) ), x -> FromEndMToHomMM( L[i], x) );
        radmaps := Flat( List( homL_iC, h -> RadendL_i * h ) );
        generators := List( homL_iapproxC, h -> Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) ) );
        radgenerators := List( radmaps, h -> Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) ) );
        Append( generators, radgenerators );
        generatorshomL_iC := List( homL_iC, h -> Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) ) );
        if Length( generators ) = 0 then
          generators := [ Flat( MatricesOfPathAlgebraMatModuleHomomorphism( ZeroMapping( L[ i ],C ) ) ) ];
        fi;
        for g in generatorshomL_iC do
          V := VectorSpace( K, generators );
          if not g in V then
            Add( generators, g );
            t := Position( generatorshomL_iC, g );
            Add( approximation, homL_iC[ t ] );
          fi;
        od;
      fi;
    od;
    
    M := DirectSumOfQPAModules( List( approximation, Source ) );
    projections := DirectSumProjections( M );
    f := Sum( List( [ 1..Length( projections ) ], i -> projections[ i ] * approximation[ i ] ) );

    return f;
end
  );

InstallMethod ( RadicalRightApproximationByAddM, 
    "for a list of modules and a module",
    true,
    [ IsList, IsInt ],
    0,
    function( modulelist, t )
        
    local   length,  moduleset,  N,  f,  K,  M,  endoN,  radendoN,  
            BradendoN,  gen_NtoM,  radhomNN,  Vgen_NtoM,  VradhomNN,  
            U,  pi,  addhoms,  i,  newM,  projections,  homs;
    
    length := Length( modulelist );
    if not t in [ 1..length ] then
        Error( "The entered integer is not in the correct interval.\n" );
    fi;
    
    moduleset := [ 1..length ];
    RemoveSet( moduleset, t );
    N := modulelist[ t ];
    f := RightApproximationByAddM( modulelist{ moduleset }, N );
    K := LeftActingDomain( N );
    M := Source( f );
    
    endoN := EndOverAlgebra( N );
    radendoN := RadicalOfAlgebra( endoN );
    BradendoN := Basis( radendoN );
    gen_NtoM := List( HomOverAlgebra( N, M ), h -> h * f );
    radhomNN := List( BradendoN, b -> FromEndMToHomMM( N, b ) );
    if Length( radhomNN ) > 0 then 
        Vgen_NtoM := List( gen_NtoM, h -> Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) ) );    
        VradhomNN := List( radhomNN, h -> Flat( MatricesOfPathAlgebraMatModuleHomomorphism( h ) ) );
        U := FullRowSpace( K, Length( VradhomNN[ 1 ] ) );
        pi := NaturalHomomorphismBySubspace( U, Subspace( U, Vgen_NtoM ) );
        addhoms := [ ];
        for i in [ 1..Length( VradhomNN ) ] do
            if not IsZero( ImageElm( pi, VradhomNN[ i ] ) ) then
                Add( addhoms, i );
            fi;
        od;
        newM := [ M ];
        for i in addhoms do
            Add( newM, N );
        od;
        M := DirectSumOfQPAModules( newM );
        projections := DirectSumProjections( M );
        homs := [ f ];
        for i in addhoms do
            Add( homs, radhomNN[ i ] );
        od;
        
        f := Sum( List( [ 1..Length( projections ) ], i -> projections[ i ] * homs[ i ] ) );
    fi;
        
    return f;
end
  );

InstallMethod ( ProjectiveResolutionOfSimpleModuleOverEndo, 
    "for a list of modules, an index of the list, and an integer",
    true,
    [ IsList, IsInt, IsInt ],
    0,
    function( modulelist, t, length )
    
    local   f,  syzygies,  n,  U,  m,  test;
    
    if length < 0 then
        Error( "The entered length is less than zero.\n" );
    fi;
    f := RadicalRightApproximationByAddM( modulelist, t );    
    syzygies := [ ];
    for n in [ 1..length - 2 ] do
        if Dimension( Source( f ) ) = 0 then
            return [ n - 1, syzygies ];
        fi;
        if IsInjective( f ) then
            return [ n, syzygies ];
        fi;    
        U := Kernel( f );
        for m in modulelist do
            repeat
                test := CommonDirectSummand( U, m );
                if test <> false then 
                    U := test[ 2 ];
                fi;
            until
              test = false;
        od;
        if Dimension( U ) = 0 then
            return [ n + 1, syzygies ];
        fi;
        Add( syzygies, U );
        f := RightApproximationByAddM( modulelist, U );
    od;
    
    return [ Concatenation( "projdim > ", String( t ) ), syzygies ];
end
  );