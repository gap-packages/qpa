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