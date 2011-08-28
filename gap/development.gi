InstallMethod( IsSelfinjective, 
   "for a finite dimension quotient of a path algebra",
   [ IsSubalgebraFpPathAlgebra ], 0,
   function( A ) 

   local KQ, rels, I, B, gb, gbb, Inj, T, num_vert, total, i;

   Inj := IndecomposableInjectiveRepresentations(A);
   T := List(Inj, M -> DimensionVector(TopOfRep(M)));
   num_vert := Length(T);
   total := [];
   for i in [1..num_vert] do
      Add(total,0);
   od;
   for i in [1..num_vert] do
      total := total + T[i];
   od;
   
   if ( num_vert = Sum(total) ) and ( ForAll(total, x -> x > 0) )  then
      return true;
   else
      return false;
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
         N := RadicalOfRep(N);
         i := i + 1;
      until
         Dimension(N) = 0;
      return i;
   fi;
end
);

InstallOtherMethod( LoewyLength, 
   "for a SubalgebraFpPathAlgebra",
   [ IsSubalgebraFpPathAlgebra ], 0,
   function( A ) 

   local N, i;

   N := IndecomposableProjectiveRepresentations(A);
   N := List(N, x -> LoewyLength(x));

   return Maximum(N);
end
);

InstallOtherMethod( LoewyLength, 
   "for a PathAlgebra",
   [ IsPathAlgebra ], 0,
   function( A ) 

   local N, i;

   if IsAcyclic(QuiverOfPathAlgebra(A)) then 
      N := IndecomposableProjectiveRepresentations(A);
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

   P := IndecomposableProjectiveRepresentations(A);
   C := [];
   for i in [1..Length(P)] do
      Add(C,DimensionVector(P[i]));
   od;

   return C;
end
);

InstallOtherMethod( CartanMatrix, 
   "for a SubalgebraFpPathAlgebra",
   [ IsSubalgebraFpPathAlgebra ], 0,
   function( A ) 

   local P, C, i;

   P := IndecomposableProjectiveRepresentations(A);
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
   [ IsSubalgebraFpPathAlgebra ], 0,
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
   [ IsSubalgebraFpPathAlgebra ], 0,
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
         Add(radlayers,DimensionVector(TopOfRep(N)));
         N := RadicalOfRep(N);
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
      N := DualOfPathAlgebraMatModule(M);
      repeat     
         Add(series,DimensionVector(TopOfRep(N)));
         N := RadicalOfRep(N);
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

InstallMethod( DimensionMatModule,
   "for a PathAlgebraMatModule",
   [ IsPathAlgebraMatModule ], 0,
   function( M );

   return Sum(DimensionVector(M));
end
);

InstallMethod( Centre,
   "for a path algebra",
   [ IsSubalgebraFpPathAlgebra ], 0,
   function( A ) 

   local Q, K, num_vert, vertices, arrows, B, cycle_list, i, j, b, 
         pos, commutators, center, zero_count, a, c, x, cycles, matrix,
         data, coeffs, fam, elms, solutions, gbb;

   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
   K := LeftActingDomain(A);
   num_vert := OrderOfQuiver(Q);
   vertices := VerticesOfQuiver(Q);
   arrows   := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+SizeOfQuiver(Q)]};
  
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

   matrix := NullMat(cycles,Length(B)*SizeOfQuiver(Q),K);

   for j in [1..cycles] do
      for i in [1..SizeOfQuiver(Q)] do
         matrix[j]{[Length(B)*(i-1)+1..Length(B)*i]} := Coefficients(B,commutators[i+(j-1)*SizeOfQuiver(Q)]); 
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

   top := TopOfRep(M); 
   dimension := 0; 
   for i in [1..Length(DimensionVector(M))] do 
      if DimensionVector(top)[i] <> 0 then 
         P := IndecomposableProjectiveRepresentations(RightActingAlgebra(M),[i]); 
         dimension := dimension + DimensionMatModule(P[1])*DimensionVector(top)[i];
      fi;
   od; 

   if dimension = DimensionMatModule(M) then 
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

   return IsProjectiveModule(DualOfPathAlgebraMatModule(M));
end
);

InstallMethod( IsSimpleModule, 
   "for a module over a quotient of a path algebra",
   [ IsPathAlgebraMatModule ], 0,
   function( M ) ; 

   if DimensionMatModule(M) = 1 then 
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

   if DimensionMatModule(RadicalOfRep(M)) = 0 then 
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

InstallMethod( 1st_Syzygy,
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

   A := RightActingAlgebra(M);
   B := CanonicalBasis(A);
   Q := QuiverOfPathAlgebra(A); 
   K := LeftActingDomain(A);
   num_vert := OrderOfQuiver(Q);
   vertices := VerticesOfQuiver(Q);
   if IsPathAlgebra(A) then 
      verticesinalg := GeneratorsOfAlgebra(A){[1..num_vert]};
      arrows   := GeneratorsOfAlgebra(A){[1+num_vert..num_vert+SizeOfQuiver(Q)]};   
   else 
      verticesinalg := GeneratorsOfAlgebra(A){[2..1+num_vert]};
      arrows   := GeneratorsOfAlgebra(A){[2+num_vert..1+num_vert+SizeOfQuiver(Q)]};
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
   G   := MinimalSetOfGenerators(M);
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
      first_syzygy := RightModuleOverQuotientOfPathAlgebra(A,big_mat); 
   fi;      

   return first_syzygy;
end
);

InstallMethod( nth_Syzygy,
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
         result := 1st_Syzygy(result);
         Print("Top of the ",i,"th syzygy: ",DimensionVector(TopOfRep(result)),"\n");
         if IsProjectiveModule(result) then 
            Print("The module has projective dimension ",i,".\n");
            break;
         fi;
      od;
   fi;

   return result;
end
);

InstallMethod( nth_SyzygyNC,
   "for a path algebra module and a positive integer",
   [ IsPathAlgebraMatModule, IS_INT ], 0,
   function( M, n ) 

   local i, result;
 
   result := ShallowCopy(M);
   for i in [1..n] do
      Print("Computing syzygy number: ",i," ....\n");
      result := 1st_Syzygy(result);
      if Dimension(result) = 0 then 
         Print("The module has projective dimension ",i-1,".\n");
         break;
      fi;
   od;

   return result;
end
);

InstallMethod( DirectSumOfPathAlgebraMatModules,
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
                     MatricesOfPathAlgebraMatModule(L[r])[i];
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
         direct_sum := RightModuleOverQuotientOfPathAlgebra(A,big_mat); 
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
         Add(list_of_projs,RightModuleHomOverPathAlgebra(direct_sum,L[i],maps));
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
         Add(list_of_incls,RightModuleHomOverPathAlgebra(L[i],direct_sum,maps));
      od;

      if Sum(dim_vect) <> 0 then 
         SetIsDirectSum(direct_sum,true);
         SetDirectSumProjections(direct_sum,list_of_projs);
         SetDirectSumInclusions(direct_sum,list_of_incls);
      else 
         SetIsDirectSum(direct_sum,false);         
      fi;
      return direct_sum;
   fi;

   return fail;
end
);

InstallMethod( PushOut,
   "for two homomorphisms starting in a common module over a path algebra",
   [ IsPathAlgebraMatModuleMap, IsPathAlgebraMatModuleMap ], 0,
   function( f, g ) 

   local B, C, BplusC, inclusions, projections, h;

   if Source(f) = Source(g) then 
      B := Range(f);
      C := Range(g);
      BplusC := DirectSumOfPathAlgebraMatModules([B,C]);
      inclusions := DirectSumInclusions(BplusC);
      h := f*inclusions[1]-g*inclusions[2];
      h := CokerProjection(h);

      return [inclusions[1]*h,inclusions[2]*h];
   else
      Print("Error: The two maps entered don't start in the same module.\n");
      return fail;
   fi;
end
);

InstallMethod( PullBack,
   "for two homomorphisms ending in a common module over a path algebra",
   [ IsPathAlgebraMatModuleMap, IsPathAlgebraMatModuleMap ], 0,
   function( f, g ) 

   local A, B, AplusB, projections, h;

   if Range(f) = Range(g) then 
      A := Source(f);
      B := Source(g);
      AplusB := DirectSumOfPathAlgebraMatModules([A,B]);
      projections := DirectSumProjections(AplusB);
      h := projections[1]*f-projections[2]*g;
      h := KerInclusion(h);

      return [h*projections[1],h*projections[2]];
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
      N1 := 1st_Syzygy(N0);
      if ModuleIsomorphismTest(M,N1) then
         return i;
      else
         N0 := N1;
      fi;
   od;
   return false;
end
);