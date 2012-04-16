InstallMethod( IsSelfinjectiveAlgebra, 
   "for a finite dimension quotient of a path algebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local fam, KQ, rels, I, B, gb, gbb, Inj, T, num_vert, total, i;

   fam := ElementsFamily(FamilyObj(A));
   if HasGroebnerBasisOfIdeal(fam!.ideal) and 
          AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then 
      Inj := IndecInjectiveModules(A);
      T := List(Inj, M -> DimensionVector(TopOfModule(M)));
      num_vert := Length(T);
      total := List([1..num_vert], x -> 0);
      for i in [1..num_vert] do
         total := total + T[i];
      od;
   
      if ( num_vert = Sum(total) ) and ( ForAll(total, x -> x > 0) )  then
         return true;
      else
         return false;
      fi;
   else
      return fail;
   fi;   
end
);

InstallOtherMethod( IsSelfinjectiveAlgebra,
   "for a path algebra",
   [ IsPathAlgebra ], 0,
   function( A )

   local Q;

   Q := QuiverOfPathAlgebra(A);
   if IsAcyclicQuiver(Q) then 
      if NumberOfArrows(Q) > 0 then 
         return false;
      else
         return true;
      fi;
   else
      return fail;
   fi;
end
);

InstallMethod( LoewyLength, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, i;

   N := M;
   i := 0;
   if Dimension(M) = 0 then 
      return 0;
   else
      repeat 
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      return i;
   fi;
end
);

InstallOtherMethod( LoewyLength, 
   "for a SubalgebraFpPathAlgebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local fam, gb, N;

   fam := ElementsFamily(FamilyObj(A));
   if HasGroebnerBasisOfIdeal(fam!.ideal) and
          AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(fam!.ideal)) then
      gb := GroebnerBasisOfIdeal(fam!.ideal);
      N := List(Nontips(gb), x -> LengthOfPath(x));
      return Maximum(N)+1;
   else
      return fail;
   fi;
end
);

InstallOtherMethod( LoewyLength, 
   "for a PathAlgebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

   local N, i;

   if IsAcyclicQuiver(QuiverOfPathAlgebra(A)) then 
      N := IndecProjectiveModules(A);
      N := List(N, x -> LoewyLength(x));
      return Maximum(N);
   else
      return fail;
   fi;
end
);

InstallOtherMethod( CartanMatrix, 
   "for a PathAlgebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   P := IndecProjectiveModules(A);
   C := [];
   for i in [1..Length(P)] do
      Add(C,DimensionVector(P[i]));
   od;

   return C;
end
);

InstallOtherMethod( CartanMatrix, 
   "for a SubalgebraFpPathAlgebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   P := IndecProjectiveModules(A);
   C := [];
   for i in [1..Length(P)] do
      Add(C,DimensionVector(P[i]));
   od;

   return C;
end
);

InstallMethod( CoxeterMatrix, 
   "for a PathAlgebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   C := CartanMatrix(A);
   if DeterminantMat(C) <> 0 then 
      return (-1)*C^(-1)*TransposedMat(C);
   else
      return fail;
   fi;
end
);

InstallOtherMethod( CoxeterMatrix, 
   "for a PathAlgebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local C;

   C := CartanMatrix(A);
   if DeterminantMat(C) <> 0 then 
      return (-1)*C^(-1)*TransposedMat(C);
   else
      return fail;
   fi;
end
);

InstallMethod( CoxeterPolynomial, 
   "for a PathAlgebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   C := CartanMatrix(A);
   if C <> fail then 
      return CharacteristicPolynomial(C);
   else
      return fail;
   fi;
end
);

InstallOtherMethod( CoxeterPolynomial, 
   "for a PathAlgebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   C := CartanMatrix(A);
   if C <> fail then 
      return CharacteristicPolynomial(C);
   else
      return fail;
   fi;
end
);


InstallMethod( RadicalSeries, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, radlayers, i;

   if Dimension(M) = 0 then 
      return DimensionVector(M);
   else
      radlayers := [];
      i := 0;
      N := M;
      repeat     
         Add(radlayers,DimensionVector(TopOfModule(N)));
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      return radlayers;
   fi;
end
);

InstallMethod( SocleSeries, 
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local N, series, socleseries, i, n;

   if Dimension(M) = 0 then 
      return DimensionVector(M);
   else
      series := [];
      i := 0;
      N := DualOfModule(M);
      repeat     
         Add(series,DimensionVector(TopOfModule(N)));
         N := RadicalOfModule(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      n := Length(DimensionVector(M));
      socleseries := [];
      for n in [1..i] do
         socleseries[n] := series[i-n+1];
      od;
      return socleseries;
   fi;
end
);

InstallOtherMethod( Dimension,
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M );

   return Sum(DimensionVector(M));
end
);

InstallMethod( Centre,
   "for a path algebra",
   [ IsQuotientOfPathAlgebra ], 0,
   function( A ) 

   local Q, K, num_vert, vertices, arrows, B, cycle_list, i, j, b, 
         pos, commutators, center, zero_count, a, c, x, cycles, matrix,
         data, coeffs, fam, elms, solutions, gbb;

   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
   K := LeftActingDomain(A);
   num_vert := NumberOfVertices(Q);
   vertices := VerticesOfQuiver(Q);
   arrows   := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+NumberOfArrows(Q)]};
  
   cycle_list := [];
   for b in B do
      if SourceOfPath(TipMonomial(b)) = TargetOfPath(TipMonomial(b)) then 
         Add(cycle_list,b);
      fi;
   od;
   commutators := [];
   center := [];
   cycles := Length(cycle_list);
   for c in cycle_list do
      for a in arrows do
         Add(commutators,c*a-a*c);
      od;
   od;

   matrix := NullMat(cycles,Length(B)*NumberOfArrows(Q),K);

   for j in [1..cycles] do
      for i in [1..NumberOfArrows(Q)] do
         matrix[j]{[Length(B)*(i-1)+1..Length(B)*i]} := Coefficients(B,commutators[i+(j-1)*NumberOfArrows(Q)]); 
      od;
   od;

   solutions := NullspaceMat(matrix);

   center := [];
   for i in [1..Length(solutions)] do
      x := Zero(A);
      for j in [1..cycles] do 
         if solutions[i][j] <> Zero(K) then 
            x := x + solutions[i][j]*cycle_list[j];
         fi;
      od;
      center[i] := x;
   od;
   return center;
end
);

InstallMethod( IsProjectiveModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local top, dimension, i, P; 

   top := TopOfModule(M); 
   dimension := 0; 
   for i in [1..Length(DimensionVector(M))] do 
      if DimensionVector(top)[i] <> 0 then 
         P := IndecProjectiveModules(RightActingAlgebra(M),[i]); 
         dimension := dimension + Dimension(P[1])*DimensionVector(top)[i];
      fi;
   od; 

   if dimension = Dimension(M) then 
      return true;
   else 
      return false;
   fi;
end
);

InstallMethod( IsInjectiveModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   return IsProjectiveModule(DualOfModule(M));
end
);

InstallMethod( IsSimpleModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   if Dimension(M) = 1 then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod( IsSemisimpleModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   if Dimension(RadicalOfModule(M)) = 0 then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod( TipMonomialandCoefficientOfVector, 
   "for a path algebra",
   [ IsAlgebra, IsCollection ], 0,
   function( A, x ) 

   local pos, tipmonomials, sortedtipmonomials, tippath, i, n;

   pos := 1;
   tipmonomials := List(x,TipMonomial);
   sortedtipmonomials := ShallowCopy(tipmonomials);
   Sort(sortedtipmonomials,\<);

   if Length(tipmonomials) > 0 then
      tippath := sortedtipmonomials[Length(sortedtipmonomials)];
      pos := Minimum(Positions(tipmonomials,tippath));
   fi;

   return [pos,TipMonomial(x[pos]),TipCoefficient(x[pos])];
end
);


InstallMethod( TipReduceVectors, 
   "for a path algebra",
   [ IsAlgebra, IsCollection ], 0,
   function( A, x ) 

   local i, j, k, n, y, s, t, pos_m, z, stop;

   n := Length(x);
   if n > 0 then 
      for k in [1..n] do
      	 for j in [1..n] do
            if j <> k then  
               s := TipMonomialandCoefficientOfVector(A,x[k]);
               t := TipMonomialandCoefficientOfVector(A,x[j]); 
               if ( s <> fail ) and ( t <> fail ) then
                  if ( s[1] = t[1] ) and ( s[2] = t[2] ) and 
                     ( s[3] <> Zero(LeftActingDomain(A)) ) then 
               	     x[j] := x[j] - (t[3]/s[3])*x[k];
                  fi;
               fi;
            fi;
      	 od;
      od;
      return x;
   fi;
end
);

#
# This has a bug !!!!!!!!!!!!!!!!!!!
#
InstallMethod( CoefficientsOfVectors, 
   "for a path algebra",
   [ IsAlgebra, IsCollection, IsList ], 0,
   function( A, x, y ) 

   local i, j, n, s, t, tiplist, vector, K, zero_vector, z;

   tiplist := [];
   for i in [1..Length(x)] do 
      Add(tiplist,TipMonomialandCoefficientOfVector(A,x[i]){[1..2]});
   od;

   vector := []; 
   K := LeftActingDomain(A);
   for i in [1..Length(x)] do
      Add(vector,Zero(K));
   od;
   zero_vector := [];
   for i in [1..Length(y)] do
      Add(zero_vector,Zero(A));
   od;

   z := y;
   if z = zero_vector then
      return vector;
   else 
      repeat
         s := TipMonomialandCoefficientOfVector(A,z);
         j := Position(tiplist,s{[1..2]});
         if j = fail then 
            return [x,y];
         fi;
         t := TipMonomialandCoefficientOfVector(A,x[j]); 
         z := z - (s[3]/t[3])*x[j];
         vector[j] := vector[j] + (s[3]/t[3]);
      until 
         z = zero_vector;
      return vector;
   fi;
end
);

InstallMethod( 1stSyzygy,
   "for a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) 

   local A, Q, K, num_vert, vertices, verticesinalg, arrows, B, BB, B_M, 
         G, MM, r, test, projective_cover, s, projective_vector, Solutions, 
         syzygy, zero_in_P, big_mat, mat, a, arrow, arrows_of_quiver, dom_a,
	 im_a, pd, pi, dim_vect, first_syzygy, V, BV, m, v, 
         cycle_list, i, j, b,  
         pos, commutators, center, zero_count, c, x, cycles, matrix,
         data, coeffs, fam, elms, solutions, gbb, run_time, BU, FlatBasisSyzygy,
         partial,temp,V_list,B_list;

   A:= RightActingAlgebra(M);   
   if Dimension(M) = 0 then 
       return ZeroModule(A);
   fi;
   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
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
#  Finding a basis of each indecomposable right projective A-module
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
#
# Finding a basis for the module M and a set of minimal generators
#
   B_M := CanonicalBasis(M);
   G   := MinimalGeneratingSetOfModule(M);
#
#  Assuming that the generators G of M is uniform.
#  Finding generators multiplied with all basis elements of A.
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
#
# Finding the matrix with the rows being the generators of M times the 
# the basis of the corresponding projective.
#
   matrix := NullMat(Length(MM),Dimension(M),K);
   for i in [1..Length(MM)] do
      matrix[i] := Flat(ExtRepOfObj(MM[i])![1]);
   od;
#
#  Finding the kernel of the projective cover as a vectorspace 
#  
   solutions := NullspaceMat(matrix);
#
# Finding the kernel as a submodule of the projective cover.
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
# Finding a uniform basis for the syzygy.
# 
   syzygy := [];
   for i in [1..num_vert] do 
      syzygy[i] := [];
   od;
   zero_in_P := [];
   for i in [1..Length(G)] do
      Add(zero_in_P,Zero(A));
   od;
   for i in [1..num_vert] do 
      for j in [1..Length(Solutions)] do 
         if Solutions[j]*verticesinalg[i] <> zero_in_P then
            Add(syzygy[i],Solutions[j]*verticesinalg[i]);
         fi;
      od;
   od;

   for i in [1..num_vert] do
      if Length(syzygy[i]) > 0 then
         syzygy[i] := TipReduceVectors(A,syzygy[i]);
      fi;
   od;

   arrows_of_quiver := GeneratorsOfQuiver(Q){[1+num_vert..num_vert+Length(ArrowsOfQuiver(Q))]};
#
# Finding the dimension vector of the syzygy.
#
   dim_vect := List(syzygy,x->Length(x));
#
# Finding a basis for the vectorspace in each vertex for the syzygy.
#
   BU := Basis(UnderlyingLeftModule(B));   
   FlatBasisSyzygy := [];
   for i in [1..num_vert] do
      if dim_vect[i] > 0 then
         partial := [];
         for j in [1..dim_vect[i]] do
            Add(partial,Flat(List(syzygy[i][j],x->Coefficients(BU,x))));
         od;
         Add(FlatBasisSyzygy,partial);
      else
         Add(FlatBasisSyzygy,[]);
      fi;
   od;
#
# Constructing the vectorspace in each vertex in the syzygy and 
# saying that the basis is the generators we have found.
#
  V_list := []; 
  for i in [1..num_vert] do
     if FlatBasisSyzygy[i] <> [] then 
        Add(V_list,VectorSpace(K,FlatBasisSyzygy[i],"basis"));
     else
        Add(V_list,Subspace(K,[]));
     fi;
   od;
   B_list := [];
   for i in [1..num_vert] do 
      Add(B_list,Basis(V_list[i],FlatBasisSyzygy[i]));
   od;
   Print("Dimension vector for syzygy: ",dim_vect,"\n");
#
#  Finding the 1st syzygy as a representation of the quiver.
#
   big_mat:=[];
   for a in arrows do
      mat := [];
      for v in verticesinalg do
         if v*a <> Zero(A) then
              dom_a := v;
         fi;
      od;
      for v in verticesinalg do
         if a*v <> Zero(A) then
            im_a := v;
         fi;
      od; 
        
      pd := Position(verticesinalg,dom_a);
      pi := Position(verticesinalg,im_a);
        
      arrow := arrows_of_quiver[Position(arrows,a)];
      if ( dim_vect[pd] = 0 ) or ( dim_vect[pi] = 0 ) then 
         mat := [dim_vect[pd],dim_vect[pi]];
      else 
         for m in syzygy[pd] do
            temp := Flat(List(m*a,x->Coefficients(BU,x)));
            temp := Coefficients(B_list[pi],temp);
            Add(mat,temp);
         od;
      fi;
      Add(big_mat,[arrow,mat]);
   od;
#
# Creating the syzygy as a representation of the quiver.
#
   if IsPathAlgebra(A) then 
      first_syzygy := RightModuleOverPathAlgebra(A,big_mat);
   else
      first_syzygy := RightModuleOverPathAlgebra(A,big_mat); 
   fi;      

   return first_syzygy;
end
);

InstallMethod( NthSyzygy,
   "for a path algebra module and a positive integer",
   [ IsPathAlgebraMatModule, IS_INT ], 0,
   function( M, n ) 

   local i, result;
 
   result := ShallowCopy(M);
   if IsProjectiveModule(M) then 
      Print("The module entered is projective.\n");
   else 
      for i in [1..n] do
         Print("Computing syzygy number: ",i," ...\n");
         result := 1stSyzygy(result);
         Print("Top of the ",i,"th syzygy: ",DimensionVector(TopOfModule(result)),"\n");
         if IsProjectiveModule(result) then 
            Print("The module has projective dimension ",i,".\n");
            break;
         fi;
      od;
   fi;

   return result;
end
);

InstallMethod( NthSyzygyNC,
   "for a path algebra module and a positive integer",
   [ IsPathAlgebraMatModule, IS_INT ], 0,
   function( M, n ) 

   local i, result;
 
   result := ShallowCopy(M);
   for i in [1..n] do
      Print("Computing syzygy number: ",i," ....\n");
      result := 1stSyzygy(result);
      if Dimension(result) = 0 then 
         Print("The module has projective dimension ",i-1,".\n");
         break;
      fi;
   od;

   return result;
end
);

InstallMethod( DirectSumOfModules,
   "for a list of modules over a path algebra",
   [ IsList ], 0,
   function( L ) 

   local n, A, K, Q, arrows, vertices, dim_list, dim_vect, i, 
         temp, j, big_mat, a, mat, origin, target, r, row_pos,
         col_pos, direct_sum, dimension, list_of_projs, 
         list_of_incls, map, maps; 

   n := Length(L);
   if n > 0 then 
      A := RightActingAlgebra(L[1]);
      K := LeftActingDomain(A);
      if ForAll(L,IsPathAlgebraMatModule) and ForAll(L,x -> RightActingAlgebra(x) = A) then 
         Q := QuiverOfPathAlgebra(OriginalPathAlgebra(A));
         arrows := ArrowsOfQuiver(Q);
         vertices := VerticesOfQuiver(Q);
         dim_list := List(L,DimensionVector);
         dim_vect := [];

         for i in [1..Length(dim_list[1])] do
            temp := 0; 
            for j in [1..Length(dim_list)] do
               temp := temp + dim_list[j][i];
            od;
            Add(dim_vect,temp);
         od;
         big_mat := [];
         i := 0;
         for a in arrows do
            i := i + 1;
            origin := Position(vertices,SourceOfPath(a));
            target := Position(vertices,TargetOfPath(a));
            mat := [];
            if dim_vect[origin] <> 0 and dim_vect[target] <> 0 then 
               mat := NullMat(dim_vect[origin],dim_vect[target], K);
               row_pos := 1;
               col_pos := 1;
               for r in [1..n] do
                  if dim_list[r][origin] <> 0 and dim_list[r][target] <> 0 then
                     mat{[row_pos..(row_pos+dim_list[r][origin]-1)]}{[col_pos..(col_pos+dim_list[r][target]-1)]} := 
                     MatricesOfPathAlgebraModule(L[r])[i];
                  fi;
                  row_pos := row_pos + dim_list[r][origin];
                  col_pos := col_pos + dim_list[r][target]; 
               od;
               mat := [a,mat];
            else 
               mat := [a,[dim_vect[origin],dim_vect[target]]];
            fi;
            Add(big_mat,mat);
         od;
      fi;
      if IsPathAlgebra(A) then 
         direct_sum := RightModuleOverPathAlgebra(A,big_mat);
      else
         direct_sum := RightModuleOverPathAlgebra(A,big_mat); 
      fi; 

      list_of_projs := [];
      for i in [1..n] do 
         map := [];
         maps := [];
         for j in [1..Length(dim_vect)] do
            if dim_vect[j] = 0 then 
               if dim_list[i][j] = 0 then 
                  map := NullMat(1,1,K);
               else
                  map := NullMat(1,dim_list[i][j],K);
               fi;
            else
               if dim_list[i][j] = 0 then 
                  map := NullMat(dim_vect[j],1,K);
               else
                  dimension := 0;
                  if i > 1 then 
                     for r in [1..i-1] do 
                        dimension := dimension + dim_list[r][j];
                     od;
                  fi;
                  map := NullMat(dim_vect[j],dim_list[i][j],K);
                  map{[dimension + 1..dimension + dim_list[i][j]]}{[1..dim_list[i][j]]} := IdentityMat(dim_list[i][j],K);
               fi;
            fi;
            Add(maps,map);
         od;
         Add(list_of_projs,RightModuleHomOverAlgebra(direct_sum,L[i],maps));
      od;

      list_of_incls := [];
      for i in [1..n] do 
         map := [];
         maps := [];
         for j in [1..Length(dim_vect)] do
            if dim_vect[j] = 0 then 
               map := NullMat(1,1,K);
            else
               if dim_list[i][j] = 0 then 
                  map := NullMat(1,dim_vect[j],K);
               else
                  dimension := 0;
                  if i > 1 then 
                     for r in [1..i-1] do 
                        dimension := dimension + dim_list[r][j];
                     od;
                  fi;
                  map := NullMat(dim_list[i][j],dim_vect[j],K);
                  map{[1..dim_list[i][j]]}{[dimension + 1..dimension + dim_list[i][j]]} := IdentityMat(dim_list[i][j],K);
               fi;
            fi;
            Add(maps,map);
         od;
         Add(list_of_incls,RightModuleHomOverAlgebra(L[i],direct_sum,maps));
      od;

      if Sum(dim_vect) <> 0 then 
         SetIsDirectSumOfModules(direct_sum,true);
         SetDirectSumProjections(direct_sum,list_of_projs);
         SetDirectSumInclusions(direct_sum,list_of_incls);
      else 
         SetIsDirectSumOfModules(direct_sum,false);         
      fi;
      return direct_sum;
   fi;

   return fail;
end
);
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
      CplusB := DirectSumOfModules([C,B]);
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
      AplusC := DirectSumOfModules([A,C]);
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

InstallMethod ( FromHomMMToEndM, 
   "for a subset of EndOverAlgebra to HomOverAlgebra",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local K, dim_vect, end_f, r, j;

   K := LeftActingDomain(Source(f));
   dim_vect := DimensionVector(Source(f));
   end_f := NullMat(Dimension(Source(f)),Dimension(Source(f)),K);
   r := 1; 
   for j in [1..Length(dim_vect)] do 
      if dim_vect[j] <> 0 then 
         end_f{[r..r+dim_vect[j]-1]}{[r..r+dim_vect[j]-1]} := f!.maps[j];
         r := r + dim_vect[j];
      fi;
   od; 

   return end_f;
end
);

InstallMethod ( FromEndMToHomMM, 
   "for a subset of EndOverAlgebra to HomOverAlgebra",
   true,
   [ IsPathAlgebraMatModule, IsMatrix ],
   0,
   function( M, mat )

   local K, dim_vect, maps, i, r;

   K := LeftActingDomain(M); 
   dim_vect := DimensionVector(M);

   maps := [];
   r := 1;
   for i in [1..Length(dim_vect)] do
      if dim_vect[i] = 0 then 
         Add(maps,NullMat(1,1,K));
      else
         Add(maps,mat{[r..r+dim_vect[i]-1]}{[r..r+dim_vect[i]-1]});
         r := r + dim_vect[i];
      fi;
   od;

   return RightModuleHomOverAlgebra(M,M,maps);
end
);

InstallMethod ( IsRightMinimal, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local B, C, BB, mat, Ann_f, radEndB;

   B   := Source(f);
   C   := Range(f);
   BB  := HomOverAlgebra(B,B);
   mat := List(BB, x -> x*f);
   mat   := List(mat,x -> Flat(x!.maps));
   Ann_f := NullspaceMat(mat);
   Ann_f := List(Ann_f,x -> LinearCombination(BB,x));
   radEndB := RadicalOfAlgebra(EndOverAlgebra(B)); 
   Ann_f := List(Ann_f, x -> FromHomMMToEndM(x));

   if ForAll(Ann_f, x -> x in radEndB) then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod ( IsLeftMinimal, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local A, B, BB, mat, Ann_f, radEndB;

   A   := Source(f);
   B   := Range(f);
   BB  := HomOverAlgebra(B,B);
   mat := List(BB, x -> f*x );
   mat   := List(mat,x -> Flat(x!.maps));
   Ann_f := NullspaceMat(mat);
   Ann_f := List(Ann_f,x -> LinearCombination(BB,x));
   radEndB := RadicalOfAlgebra(EndOverAlgebra(B)); 
   Ann_f := List(Ann_f, x -> FromHomMMToEndM(x));

   if ForAll(Ann_f, x -> x in radEndB) then 
      return true;
   else
      return false;
   fi;
end
);

InstallMethod ( IsSplitMonomorphism, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local B, C, CB, id_B, mat, flat_id_B, split_f;

   B   := Source(f);
   C   := Range(f);
   if IsInjective(f) then 
      CB  := HomOverAlgebra(C,B);
      if Length(CB) = 0 then 
         return false;
      else 
         mat := [];
         mat := List(CB, x -> f*x);
         id_B := IdentityMapping(B); 
         mat   := List(mat,x -> Flat(x!.maps));
         flat_id_B := Flat(id_B!.maps); 
         split_f := SolutionMat(mat,flat_id_B);

         if split_f <> fail then 
            split_f := LinearCombination(CB,split_f);
            SetIsSplitEpimorphism(split_f,f);
            SetIsSplitMonomorphism(f,split_f);
            return split_f;
         else
            SetIsSplitMonomorphism(f,false);
            return false;
         fi;
      fi;
   else
      return false;
   fi;
end
);

InstallMethod ( IsSplitEpimorphism, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local B, C, CB, id_C, mat, flat_id_C, split_f;

   B   := Source(f);
   C   := Range(f);
   if IsSurjective(f) then 
      CB  := HomOverAlgebra(C,B);
      if Length(CB) = 0 then 
         return false;
      else 
         mat := List(CB, x -> x*f );
         id_C := IdentityMapping(C); 
         mat   := List(mat,x -> Flat(x!.maps));
         flat_id_C := Flat(id_C!.maps); 
         split_f := SolutionMat(mat,flat_id_C);

         if split_f <> fail then 
            split_f := LinearCombination(CB,split_f);
            SetIsSplitMonomorphism(split_f,f);
            SetIsSplitEpimorphism(f,split_f);
            return split_f;
         else
            SetIsSplitEpimorphism(f,false);
            return false;
         fi;
      fi;
   else
      return false;
   fi;
end
);

InstallMethod ( MoreRightMinimalVersion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local B, g, A, HomBA, HomAA, n, i, j, gg, hh;

   B := Source(f);
   g := KernelInclusion(f);
   A := Source(g);
   HomBA  := HomOverAlgebra(B,A);
   if Length(HomBA) = 0 then 
      return f;
   else 
      HomAA  := HomOverAlgebra(A,A);
      n := Maximum(Concatenation(DimensionVector(A),DimensionVector(B)));
      for i in [1..Length(HomBA)] do 
         for j in [1..Length(HomAA)] do
            gg := (g*HomBA[i]*HomAA[j])^n;
            if gg <> ZeroMapping(A,A) then
               hh :=  (HomBA[i]*HomAA[j]*g)^n;
               return [KernelInclusion(hh)*f,ImageInclusion(hh)*f];
            fi;
         od;
      od;
      return f;
   fi;
end
);

InstallMethod ( MoreLeftMinimalVersion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local B, g, C, HomCB, HomCC, n, i, j, gg, hh, t;

   B := Range(f);
   g := CoKernelProjection(f);
   C := Range(g);
   HomCB  := HomOverAlgebra(C,B);
   if Length(HomCB) = 0 then 
      return f;
   else 
      HomCC  := HomOverAlgebra(C,C);
      n := Maximum(Concatenation(DimensionVector(B),DimensionVector(C)));
      for i in [1..Length(HomCC)] do 
         for j in [1..Length(HomCB)] do
            gg := (HomCC[i]*HomCB[j]*g)^n;
            if gg <> ZeroMapping(C,C) then
               hh :=  (g*HomCC[i]*HomCB[j])^n;
               t := IsSplitMonomorphism(KernelInclusion(hh));
               return [f*t,f*ImageProjection(hh)];
            fi;
         od;
      od;
      return f;
   fi;
end
);

InstallMethod ( RightMinimalVersion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local Bprime, g, L;

   Bprime := [];
   g := f;
   repeat
      L := MoreRightMinimalVersion(g);
      if L <> g then 
         g := L[1];
         Add(Bprime,Source(L[2]));
      fi;
   until 
      L = g;

   SetIsRightMinimal(g,true);

   return [g,Bprime];
end
);

InstallMethod ( LeftMinimalVersion, 
   "for a PathAlgebraMatModuleMap",
   true,
   [ IsPathAlgebraMatModuleHomomorphism ],
   0,
   function( f )

   local Bprime, g, L;

   Bprime := [];
   g := f;
   repeat
      L := MoreLeftMinimalVersion(g);
      if L <> g then 
         g := L[1];
         Add(Bprime,Range(L[2]));
      fi;
   until 
      L = g;

   SetIsLeftMinimal(g,true);

   return [g,Bprime];
end
);


InstallMethod ( MinimalRightApproximation, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( M, C )

   local K, HomMC, EndM, radEndM, radHomMC, i, j, FlatHomMC, 
         FlatradHomMC, V, BB, W, f, VoverW, B, gens, approx, approxmap;

   if RightActingAlgebra(M) = RightActingAlgebra(C) then 
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
         approx := DirectSumOfModules(approx);
         approxmap := ShallowCopy(DirectSumProjections(approx));
         for i in [1..Length(approxmap)] do
            approxmap[i] := approxmap[i]*gens[i];
         od;
         approxmap := Sum(approxmap);
         return RightMinimalVersion(approxmap)[1];
      fi;
   else
      Error(" the two modules entered into MinimalRightApproximation are not modules over the same algebra.");
      return fail;
   fi;
end
);

InstallMethod ( MinimalLeftApproximation, 
   "for two PathAlgebraMatModules",
   true,
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ],
   0,
   function( C, M )

   local K, HomCM, EndM, radEndM, radHomCM, i, j, FlatHomCM, 
         FlatradHomCM, V, BB, W, f, VoverW, B, gens, approx, approxmap;

   if RightActingAlgebra(M) = RightActingAlgebra(C) then 
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
         approx := DirectSumOfModules(approx);
         approxmap := ShallowCopy(DirectSumInclusions(approx));
         for i in [1..Length(approxmap)] do
            approxmap[i] := gens[i]*approxmap[i];
         od;
         approxmap := Sum(approxmap);
         return LeftMinimalVersion(approxmap)[1];
      fi;
   else
      Error(" the two modules entered into MinimalLeftApproximation are not modules over the same algebra.");
      return fail;
   fi;
end
);

#######################################################################
##
#O  SupportModuleElement(<m>)
##
##  Given an element  <m>  in a PathAlgebraMatModule, this function
##  finds in which vertices this element is support, that is, has 
##  non-zero coordinates. 
##
InstallMethod ( SupportModuleElement, 
   "for an element in a PathAlgebraMatModule",
   true,
   [ IsRightAlgebraModuleElement ],
   0,
   function( m )

   local fam, A, Q, idempotents, support;
#
#  Finding the algebra over which the module  M, in which  m is an 
#  element, is a module over. 
# 
   if IsPathModuleElem(ExtRepOfObj(m)) then 
      fam := CollectionsFamily(FamilyObj(m))!.rightAlgebraElementsFam;   
      if "wholeAlgebra" in NamesOfComponents(fam) then 
         A := fam!.wholeAlgebra;
      elif "pathRing" in NamesOfComponents(fam) then 
         A := fam!.pathRing;
      fi;
   else 
      return fail;
   fi;
#
#  Finding the support of m.
#
   Q := QuiverOfPathAlgebra(OriginalPathAlgebra(A));
   idempotents := List(VerticesOfQuiver(Q), x -> x*One(A));
   support := Filtered(idempotents, x -> m^x <> Zero(m));

   return support;
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

#######################################################################
##
#O  VertexPosition(<elm>)
##
##  This function assumes that the input is a residue class of a trivial
##  path in finite dimensional quotient of a path algebra, and it finds  
##  the position of this trivial path/vertex in the list of vertices for
##  the quiver used to define the algebra.
##
InstallMethod ( VertexPosition, 
   "for an element in a quotient of a path algebra",
   true,
   [ IsElementOfQuotientOfPathAlgebra ],
   0,
   function( elm )

   return elm![1]![2][1]!.gen_pos;
end
);

#######################################################################
##
#O  VertexPosition(<elm>)
##
##  This function assumes that the input is a trivial path in finite 
##  dimensional a path algebra, and it finds the position of this 
##  trivial path/vertex in the list of vertices for the quiver used 
##  to define the algebra.
##
InstallOtherMethod ( VertexPosition, 
   "for an element in a PathAlgebra",
   true,
   [ IsElementOfMagmaRingModuloRelations ],
   0,
   function( elm )

   if "pathRing" in NamesOfComponents(FamilyObj(elm)) then 
      return elm![2][1]!.gen_pos;
   else
      return false;
   fi;
end
);

#######################################################################
##
#O  HomFromProjective(<m>,<M>)
##
##  This function checks if the element  <m>  is an element of  <M>  and
##  if the element  <m>  is uniform. If so, the function constructs the 
##  map from the indecomposable projective generated by the vertex of 
##  support to the module  <M>  given by sending the vertex to the 
##  element  <m>. 
##
InstallMethod ( HomFromProjective, 
   "for an element in a PathAlgebraMatModule and the PathAlgebraMatModule",
   IsElmsColls,
   [ IsRightAlgebraModuleElement, IsPathAlgebraMatModule ],
   0,
   function( m, M )

   local A, K, num_vert, support, pos, P, B, mats, zeros, i, f;
#
# Finding the algebra acting on M and the number of vertices.
#
   A := RightActingAlgebra(M);
   K := LeftActingDomain(M);
   num_vert := Length(m![1]![1]);
#   
# Finding the support of m.
#
   support := SupportModuleElement(m);
#
# If the element  m  is uniform, we proceed to find the map.
# 
   if Length(support) = 1 then 
#
# And we need the "number of" the vertex ...
#
      pos := VertexPosition(support[1]);
#
# And the basis for the correct projective module, and the proj. module itself
#
	  B := BasisOfProjectives(A)[pos];
	  P := IndecProjectiveModules(A)[pos];
#
# Then we calculate the matrices for the homomorphism
#
      mats  := List([1..num_vert],x -> List(B[x], y -> ExtRepOfObj(m^y)![1][x]));
      zeros := Filtered([1..num_vert], x -> DimensionVector(P)[x] = 0);
      for i in zeros do 
         if DimensionVector(M)[i] <> 0 then         
            mats[i] := NullMat(1,DimensionVector(M)[i],K);
         else 
            mats[i] := NullMat(1,1,K);
         fi;
      od;
#
# Construct the homomorphism
#
      f := RightModuleHomOverAlgebra(P,M,mats);
	  return f;
   else
      return fail;
   fi;
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
      PN := DirectSumOfModules(PN);
      projections := DirectSumProjections(PN);

      return projections*maps;;
   fi;
end
);

#######################################################################
##
#O  ExtOverAlgebra(<M>,<N>)
##
##  This function returns the kernel of the projective cover 
##  Omega(<M>) --> P(M) and a basis of Ext^1(<M>,<N>) inside 
##  Hom(Omega(<M>),<N>) if the group is non-zero, otherwise it returns
##  an empty list.
##
InstallMethod( ExtOverAlgebra, 
   "for two PathAlgebraMatModule's",
   true, 
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
   function( M, N )

   local K, f, g, PM, syzygy, G, H, Img1, zero, genssyzygyN, VsyzygyN, 
         Img, gensImg, VImg, pi, ext, preimages, homvecs, dimsyz, dimN, vec, t, 
         l, i, H2;
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
         return [];
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
            return [];
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
            for vec in preimages do # iterate on each homomorphism
               t := 0;
               l := []; # to store the maps for one homomorphism
               for i in [1..Size(dimsyz)] do    # matrix i in the hom. is a (dimsyz[i] x dimN[i])-matrix
                  if ( dimsyz[i] = 0 ) then 
                     if ( dimN[i] = 0 ) then 
                        Add(l,[vec{[t+1]}]);
                        t := t + 1;
                     else
                        Add(l,[vec{[1..dimN[i]+t]}]);
                        t := t + dimN[i];
                     fi;
                  else
                     if ( dimN[i] = 0 ) then 
                        Add(l,List([1..dimsyz[i]], x -> vec{[x+t]}));
                        t := t + dimsyz[i];
                     else
                        Append(l,[List([1..dimsyz[i]],x-> vec{(x-1)*dimN[i]+[1..dimN[i]] + t})]);
	                    t := t + dimsyz[i]*dimN[i];
                     fi;
                  fi;
               od;
               Append(homvecs,[l]);
            od;
#
# Making homomorphisms of the elements
#
            H2 := List(homvecs, x -> RightModuleHomOverAlgebra(syzygy,N,x));

            return [g,H2];
         fi;
      fi;
   fi;
end
);

#######################################################################
##
#O  ExtOverAlgebraAdd(<M>,<N>)
##
##  This function returns the kernel of the projective cover 
##  Omega(<M>) --> P(M) and a basis of Ext^1(<M>,<N>) inside 
##  Hom(Omega(<M>),<N>) if the group is non-zero, otherwise it returns
##  an empty list.
##
InstallMethod( ExtOverAlgebraAdd, 
   "for two PathAlgebraMatModule's",
   true, 
   [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ], 0,
   function( M, N )

   local K, f, g, PM, syzygy, G, H, Img1, zero, genssyzygyN, VsyzygyN, 
         Img, gensImg, VImg, pi, ext, preimages, homvecs, dimsyz, dimN, vec, t, 
         l, i, H2, coefficients;
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
         return [];
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
            return [g,[]];
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
            for vec in preimages do # iterate on each homomorphism
               t := 0;
               l := []; # to store the maps for one homomorphism
               for i in [1..Size(dimsyz)] do    # matrix i in the hom. is a (dimsyz[i] x dimN[i])-matrix
                  if ( dimsyz[i] = 0 ) then 
                     if ( dimN[i] = 0 ) then 
                        Add(l,[vec{[t+1]}]);
                        t := t + 1;
                     else
                        Add(l,[vec{[1..dimN[i]+t]}]);
                        t := t + dimN[i];
                     fi;
                  else
                     if ( dimN[i] = 0 ) then 
                        Add(l,List([1..dimsyz[i]], x -> vec{[x+t]}));
                        t := t + dimsyz[i];
                     else
                        Append(l,[List([1..dimsyz[i]],x -> vec{(x-1)*dimN[i]+[1..dimN[i]] + t})]);
	                    t := t + dimsyz[i]*dimN[i];
                     fi;
                  fi;
               od;
               Append(homvecs,[l]);
            od;
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
#O  AlmostSplitSequence(<M>,<N>)
##
##  This function finds the almost split sequence ending in the module
##  <M>, if the module is not projective. It returns fail otherwise. 
##  The almost split sequence is returned as a pair of maps, the 
##  monomorphism and the epimorphism. The function assumes that the 
##  module  <M>  is indecomposable, and the range of the epimorphism
##  is a module that is isomorphic to the input, not necessarily 
##  identical. 
##
InstallMethod( AlmostSplitSequence, 
   "for a PathAlgebraMatModule",
   true, 
   [ IsPathAlgebraMatModule ], 0,
   function( M )

   local K, DTrM, f, g, PM, syzygy, G, H, Img1, zero, 
         genssyzygyDTrM, VsyzygyDTrM, Img, gensImg, VImg, 
         stop, test, ext, preimages, homvecs, dimsyz, dimDTrM, 
         EndDTrM, radEndDTrM, nonzeroext, temp, L, pos, i;
#
# Add test of input.
#
   K := LeftActingDomain(M);
   if IsProjectiveModule(M) then 
      return fail;
   else 
      DTrM := DTr(M);
#
# creating a short exact sequence 0 -> Syz(M) -> P(M) -> M -> 0
# f: P(M) -> M, g: Syz(M) -> P(M)  
#
      f := ProjectiveCover(M);
      g := KernelInclusion(f);
      PM := Source(f);
      syzygy := Source(g);
#
# using Hom(-,DTrM) on the s.e.s. above
#
      G := HomOverAlgebra(PM,DTrM);
      H := HomOverAlgebra(syzygy,DTrM);
#
# Making a vector space of Hom(Syz(M),DTrM)
# by first rewriting the maps as vectors
#
      genssyzygyDTrM := List(H, x -> Flat(x!.maps));
      VsyzygyDTrM := VectorSpace(K, genssyzygyDTrM);
#
# finding a basis for im(g*)
# first, find a generating set of im(g*)
# 
      Img1 := g*G;
#
# removing 0 maps by comparing to zero = Zeromap(syzygy,DTrM)
#
      zero := ZeroMapping(syzygy,DTrM);
      Img  := Filtered(Img1, x -> x <> zero);
#
# Rewriting the maps as vectors
#
      gensImg := List(Img, x -> Flat(x!.maps));
#
# Making a vector space of <Im g*>
      VImg := Subspace(VsyzygyDTrM, gensImg);  
#
# Finding a non-zero element in Ext1(M,DTrM)
#
      i := 1;
      stop := false;
      repeat 
         test := Flat(H[i]!.maps) in VImg;
         if test then 
            i := i + 1;
         else 
            stop := true;
         fi;
      until stop;
      nonzeroext := H[i];
#
# Finding the radical of End(DTrM)
#
      EndDTrM := EndOverAlgebra(DTrM);
      radEndDTrM := RadicalOfAlgebra(EndDTrM);
      radEndDTrM := List(BasisVectors(Basis(radEndDTrM)), x -> FromEndMToHomMM(DTrM,x));
#
# Finding an element in the socle of Ext^1(M,DTrM)
#
      temp := nonzeroext;

      L := List(temp*radEndDTrM, x -> Flat(x!.maps) in VImg);
      while not ForAll(L, x -> x = true) do
         pos := Position(L,false);
         temp := temp*radEndDTrM[pos];
         L := List(temp*radEndDTrM, x -> Flat(x!.maps) in VImg);
      od;
#
# Constructing the almost split sequence in Ext^1(M,DTrM)
#
      ext := PushOut(g,temp);

      return [ext[1],CoKernelProjection(ext[1])];
   fi;
end
);

