#######################################################################
##
#O  UnitForm( <B> )
##
##  A unitform is identitied with a symmetric square matrix, where the 
##  entries along the diagonal are all 2. The bilinear form associated
##  to the unitform given by a matrix  B  is defined for two vectors  
##  x  and  y  as:
##                    x*B*y^T
##  The quadratic form associated to the unitform given by a matrix  B
##  is defined for a vector  x  as:
##                 (1/2)*x*B*x^T
##  The bilinear and the quadratic forms associated to a unitform  B  are
##  attributes of the unitform and they are given by 
##         BilinearFormOfUnitForm(B)
##  and
##         QuadraticFormOfUnitForm(B).
##  The matrix of a unitform is also given as an attribute by
##         SymmetricMatrixOfUnitForm(B).
##
InstallMethod( UnitForm, 
    "for a matrix",
    [ IsMatrix ],
        
    function( B )
    local dimB, bilinearform, quadraticform, unitform;
    
    #
    # Checking the input.
    #
    dimB := DimensionsMat(B);
    if dimB[1] <> dimB[2] then
        Error("the entered matrix is not a square matrix,\n");
    fi;
    if not ForAll(Flat(B), IS_INT) then
        Error("the entered matrix is not an integer matrix,\n");
    fi;
    if not ForAll([1..dimB[1]], x -> B[x][x] = 2) then
        Error("the matrix doesn't have 2's along the diagonal,\n");
    fi;
    if not B = TransposedMat(B) then 
        Error("the entered matrix is not a symmetric matrix,\n");
    else
        #
        # if the input is OK, define the bilinear and quadratic form, 
        # and set the associated symmetric matrix.
        #
        bilinearform := function( x, y );
            
            if not ( IsVector(x) or IsVector(y) ) then
                Error("the arguments of the bilinear form must be in the category IsVector,\n");
            fi;
            if Length(x) <> Length(y) then 
                Error("the arguments are not vectors of the same length,\n");
            fi;
            if Length(x) <> Length(B[1]) then 
                Error("the arguments don't have correct size,\n");
            fi;
            
            return x*B*y;
        end;
        quadraticform := function( x )
            if not IsVector(x) then
                Error("the argument of the quadratic form must be in the category IsVector,\n");
            fi;
            if Length(x) <> Length(B[1]) then 
                Error("the argument doesn't have the correct size,\n");
            fi;
            
            return (1/2)*x*B*x;
        end;
        
        unitform := Objectify( NewType( ElementsFamily( FamilyObj( B ) ), 
                            IsUnitForm and IsUnitFormRep ), rec( type := "undetermined" ));
        SetSymmetricMatrixOfUnitForm(unitform, B);
        SetBilinearFormOfUnitForm(unitform, bilinearform);
        SetQuadraticFormOfUnitForm(unitform, quadraticform);
        return unitform;
    fi;
end
);

#######################################################################
##
#M  ViewObj( <B> )
##
##  This function defines how View prints a UnitForm  <B>.
##
InstallMethod( ViewObj, 
    "for a UnitForm",
    true,
    [ IsUnitForm ], NICE_FLAGS+1,
    function ( B );

    View(SymmetricMatrixOfUnitForm(B));
end
);

#######################################################################
##
#O  TitsUnitFormOfAlgebra( <A> )
##
##  This function returns the Tits unitform associated to a finite 
##  dimensional quotient of a path algebra given that the underlying
##  quiver has no loops or minimal relations that starts and ends in 
##  the same vertex. That is, then it returns a symmetric matrix  B  
##  such that for  x = (x_1,...,x_n)
##  (1/2)*(x_1,...,x_n)B(x_1,...,x_n)^T = 
##      \sum_{i=1}^n x_i^2 - \sum_{i,j} dim_k Ext^1(S_i,S_j)x_ix_j 
##                + \sum_{i,j} dim_k Ext^2(S_i,S_j)x_ix_j 
##  where  n  is the number of vertices in  Q. 
##
InstallMethod( TitsUnitFormOfAlgebra,
    "for a QuotientOfPathAlgebra",
    [ IsAdmissibleQuotientOfPathAlgebra ],
        
    function( A )
    local fam, I, KQ, Q, arrows, gens, JI, IJ, JI_IJ, Aprime, IAprime,
          famAprime, gb, ImodJI_IJ, gensforImodJI_IJ, ext2matrix, pids, 
          i, j, temp;

    if not IsFiniteDimensional(A) then
        Error("the entered algebra is not finite dimensional,\n");
    fi;
    # 
    # Checking if the quiver has any loops.
    #    
    KQ := OriginalPathAlgebra(A); 
    Q := QuiverOfPathAlgebra(A);
    if not ForAll([1..NumberOfVertices(Q)], i -> AdjacencyMatrixOfQuiver(Q)[i][i] = 0) then
        Error("the quiver of the algebra has at least one loop,\n");
    fi;
    # 
    # Finding elements generating the ideal  JI + IJ and primitive orthogonal 
    # idempotents which sum is the identity. 
    # 
    fam := ElementsFamily(FamilyObj(A)); 
    I := fam!.ideal;
    arrows := ArrowsOfQuiver(Q);
    gens := GeneratorsOfTwoSidedIdeal(I);
    JI := Filtered(Flat(List(gens, g -> List(arrows, a -> a*g))), x -> x <> Zero(x));
    IJ := Filtered(Flat(List(gens, g -> List(arrows, a -> g*a))), x -> x <> Zero(x));
    JI_IJ := Concatenation(JI,IJ); 
    if Length(JI_IJ) = 0 then 
        Aprime := KQ;   
        ImodJI_IJ := gens;
        pids := List(VerticesOfQuiver(Q), v -> v*One(KQ));
    else
        Aprime := KQ/JI_IJ;  
        famAprime := ElementsFamily(FamilyObj(Aprime));
        IAprime := famAprime!.ideal;
        gb := GroebnerBasisOfIdeal(IAprime);
        ImodJI_IJ := Unique(List(gens, x -> 
                             ElementOfQuotientOfPathAlgebra(famAprime, CompletelyReduce(gb,x),true)));
        pids := List(VerticesOfQuiver(Q), v -> 
                     ElementOfQuotientOfPathAlgebra(famAprime,CompletelyReduce(gb,v*One(KQ)),true));
    fi;
    gensforImodJI_IJ := List([1..Length(pids)], x -> List([1..Length(pids)], y -> []));
    ext2matrix := NullMat(Length(pids),Length(pids));
    for i in [1..Length(pids)] do
        for j in [1..Length(pids)] do
            gensforImodJI_IJ[i][j] := Filtered(pids[i]*ImodJI_IJ*pids[j], y -> y <> Zero(y));
            gensforImodJI_IJ[i][j] := BasisVectors(Basis(Subspace(Aprime,gensforImodJI_IJ[i][j])));
            ext2matrix[i][j] := Length(gensforImodJI_IJ[i][j]);
        od; 
    od;
    if not ForAll([1..Length(pids)], i -> ext2matrix[i][i] = 0) then
        Error("the algebra has at least one minimal relation starting and ending in the same vertex,\n");
    fi;    
    temp := IdentityMat(NumberOfVertices(Q)) - AdjacencyMatrixOfQuiver(Q) + ext2matrix;
    
    return UnitForm(temp + TransposedMat(temp)); 
end
);  

InstallMethod( TitsUnitFormOfAlgebra,
    "for a PathAlgebra",
    [ IsPathAlgebra ],
        
    function( A )
    local Q, temp; 
    
    Q := QuiverOfPathAlgebra(A);
    temp := IdentityMat(NumberOfVertices(Q)) - AdjacencyMatrixOfQuiver(Q);
    return UnitForm(temp + TransposedMat(temp)); 
end
);

#######################################################################
##
#O  ReflectionByBilinearForm( <B>, <i>, <z> )
##
##  Computing the reflection of a vector  z  with respect to a bilinear form 
##  given by a matrix  B  in the coordinate  i.
##
InstallMethod( ReflectionByBilinearForm,
     "for bilinear form, integer and a vector",
     [ IsUnitForm, IS_INT, IsVector ],
        
     function( B, i, z )
     local q, n, e_i; 
     
     q := BilinearFormOfUnitForm(B); 
     n := DimensionsMat(SymmetricMatrixOfUnitForm(B))[1];
     e_i := ShallowCopy(Zero(SymmetricMatrixOfUnitForm(B)[1]));
     e_i[i] := 1; 
     if Length(z) = n then
         return z - q(z,e_i)*e_i;
     else
         Error("the entered vector for reflection was not of correct size,\n");
     fi;
end
);

#######################################################################
##
#O  IsWeaklyNonnegativeUnitForm( <B> )
##
##  This function returns true if the unit form is weakly non-negative and false
##  otherwise. It is based on the algorithm given by A. Dean and J. A. de la
##  Pena in "Algorithms for quadratic forms", Linear algebra and its 
##  applications, 235 (1996) 35-46. Comments in the code are using the same
##  notation as in this paper.
##
InstallMethod( IsWeaklyNonnegativeUnitForm,
    "for sets",
    [ IsUnitForm ],
        
    function( B )
    local q, A, Bmatrix, n, simples, i, temp, C, flag, test_one, test_two, v, R, r, j;
    
    q := BilinearFormOfUnitForm(B);
    A := SymmetricMatrixOfUnitForm(B);
    n := DimensionsMat(A)[1];
    simples := IdentityMat(n);
    C := ShallowCopy(simples);  # this is C_1.
    flag := true;
    while flag do
        for v in C do     # performing the test in (2) (a)
            if ForAny(simples, e -> q(v,e) <= -3) then 
                return false;
            fi;
        od;           
        for v in C do     # performing the test in (2) (b)
            if ForAny(simples, e -> ( q(v,e) = -2 ) and not ForAll((v+e)*A, x -> x >= 0 ) ) then 
                return false;
            fi;
        od;          
        #
        #  Now we know that the procedure didn't fail.
        #  Construct  R_s, i.e. step (3) in the algorithm.
        #
        R := Filtered(C, v -> ( ForAll(v, y -> y <= 6) and ForAny(simples, e -> q(v,e) = -1 ) ) ); 
        if Length(R) = 0 then       #  step (4) 
            flag := false;
            return true;  # C_{s+1} is empty, return "successful".
        else
            C := [];
            for r in R do  # construct C_{s+1}, i.e. step (5) in the algorithm. 
                j := Position(simples, First(simples, s -> q(r,s) = -1 ));
                Add(C,ReflectionByBilinearForm(B,j,r));  
            od;
        fi;
    od;
end
);

#######################################################################
##
#O  IsWeaklyPositiveUnitForm( <B> )
##
##  This function returns true if the unit form is weakly positive and false
##  otherwise. It is based on the algorithm given by Hans-Joachim von Hoehne
##  "Edge reduction for unit forms", Arch. Math. vol. 65 (1995) 300-302. 
##  Comments in the code are using the same notation as in this paper. In addition
##  the function computes the roots of the unit form when true is returned.
##  
InstallMethod( IsWeaklyPositiveUnitForm,
    "for a UnitForm",
    [ IsUnitForm ],
        
    function( B )
    local q, n, flag, reductionpairs, m, testI, roots, temp, testII,
          breakflag, j, i, reductionnumber, lastcolumn, lastrow, qq, s; 
    
    q := ShallowCopy(SymmetricMatrixOfUnitForm(B));    
    flag := false;
    reductionpairs := [];
    m := DimensionsMat(q)[1];
    while not flag do
        n := DimensionsMat(q)[1];
        #
        #  Performing test I) from the paper.
        #
        testI := ForAll(Flat(List([1..n], j -> List([1..j-1], i -> q[i][j]))), x -> x >= 0);
        if testI then
            roots := ShallowCopy(IdentityMat(m));
            for i in [Length(reductionpairs), Length(reductionpairs) - 1..1] do
                temp := roots[reductionpairs[i][1]]{[1..m + i - Length(reductionpairs) - 1]} + 
                        roots[reductionpairs[i][2]]{[1..m + i - Length(reductionpairs) - 1]};
                for j in [1..Length(roots)] do
                    roots[j] := roots[j]{[1..m + i - Length(reductionpairs) - 1]} + 
                                    roots[j][m + i - Length(reductionpairs)]*temp; 
                od;
            od;
            SetPositiveRootsOfUnitForm(B,roots);
            B!.type := "weakly_positive";
            return true;
        fi;
        #
        #  Performing test II) from the paper.
        #
        testII := ForAny(Flat(List([1..n], j -> List([1..j-1], i -> q[i][j]))), x -> x <= -2);
        if testII then
            return false;
        fi;
        breakflag := false;        
        for j in [1..n] do
            for i in [1..j] do
                if q[i][j] = -1 then 
                    Add(reductionpairs, [i,j]);
                    m := m + 1;
                    breakflag := true;
                    break;
                fi;
            od;
            if breakflag then 
                break;
            fi;
        od;
        #
        #  Performing reduction of unit form.
        #
        if breakflag then 
            reductionnumber := Length(reductionpairs); 
            i := reductionpairs[reductionnumber][1];
            j := reductionpairs[reductionnumber][2];
            lastcolumn := TransposedMat(q)[i] + TransposedMat(q)[j];
            lastrow := ShallowCopy(q[i] + q[j]);
            qq := List(q, r -> ShallowCopy(r)); 
            for s in [1..n] do
                Add(qq[s],lastcolumn[s]);
            od;
            Add(lastrow,2);
            Add(qq,lastrow);
            qq[i][j] := 0;
            qq[j][i] := 0;
            q := qq; 
        fi;
    od;
    
    return false; 
end
);




#######################################################################
##
#O  EulerBilinearFormOfAlgebra( <A> )
##
##  This function returns the Euler (non-symmetric) bilinear form 
##  associated to a finite dimensional (basic) quotient of a path algebra <A>.
##  That is,  it returns a bilinear form (function) defined by
##  <x,y> = x*CartanMatrix^(-1)y
##  It makes sense only in case <A> is of finite global dimension.
##  (Recall that in QPA the rows of the CartanMatrix are the dimension 
##  vectors of projectives).
##
InstallMethod( EulerBilinearFormOfAlgebra,
    "for a PathAlgebra or a QuotientOfPathAlgebra",
    [ IsQuiverAlgebra ],
        
    function( A )
    local C, bilinearform;
	
	C := CartanMatrix(A);
        # Operation CartanMatrix checks if A is: 
	# FiniteDimensional and (PathAlgebra or AdmissibleQuotientOfPathAlgebra)
        # in particular, if A is basic	
	if C = fail then
	   Print("Unable to determine the Cartan matrix!\n");
           return fail;
	fi;
	
	if not Determinant(C) in [-1,1] then
	    Print("The Cartan matrix is not invertible over Z!\n");
	    return fail;
	fi;
	
	bilinearform := function( x, y );
            
            if not ( IsVector(x) or IsVector(y) ) then
                Error("the arguments of the bilinear form must be in the category IsVector,\n");
            fi;
            if Length(x) <> Length(y) then 
                Error("the arguments are not vectors of the same length,\n");
            fi;
            if Length(x) <> Length(C[1]) then 
                Error("the arguments don't have correct size,\n");
            fi;
            
            return x*Inverse(C)*y;
       end;
		
    return bilinearform; 
end
);  # EulerBilinearFormOfAlgebra 