#############################################################################
## 
#O  PathsOfLengthTwo ( <q> ) Return : 
##
##  This function return a list of all paths of length two in q, sorted by <.
##  Fails with error message if q is not a Quiver object.
##
##  TODO: This could be made more efficient for large quivers with few
##  paths of length 2 by forming the adjacency matrix, squaring it,
##  and then checking just the source vertices indicated by the
##  result.  I have no idea whether that's worth the effort.
##
InstallMethod( PathsOfLengthTwo, 
    "for a quiver",
    [ IsQuiver ], 0,
    function( q ) 
 
    local p2, q0, v1, a1, arrow1, v2, a2, arrow2;

    if (not IsQuiver(q)) then
        Print("PathsOfLengthTwo(): argument is not a quiver\n");
        return fail;
    fi;

    # p2 will hold the list of all paths of length 2
    p2 := [];
  
    # q0 is a list of all vertices of q
    q0 := VerticesOfQuiver(q);
  
    for v1 in q0 do
        # a1 is all the arrows out of v1
        a1 := OutgoingArrowsOfVertex(v1);
        for arrow1 in a1 do
            v2 := TargetOfPath(arrow1);
            # a2 is all the arrows out of v2
            a2 := OutgoingArrowsOfVertex(v2);
            for arrow2 in a2 do
                Add(p2,arrow1*arrow2);
            od;
        od;
    od;

    Sort(p2);
    return p2;
end
); #end: PathsOfLengthTwo()

#############################################################################
##  
#O  IsQuadraticIdeal ( <gens> ) 
##  
##  This function returns true, if I is a quadratic ideal; false otherwise. 
##  Checks whether the generators passed in have the form
##
##              k1*p1 + k2*p2 + ...
##
##  where each ki is in k and each pi is a path of length 2.  The check
##  on the k's is implicit; nothing should work if they aren't right.
##
InstallMethod( IsQuadraticIdeal, 
    "for a list of elements in a path algebra",
    [ IsHomogeneousList ], 0,
    function( gens ) 
 
    local pa, L, terms_occurring, g;

    if Length(gens) = 0 then 
        Print("IsQuadraticIdeal: The entered list of elements is empty, so quadratic ideal.\n");
        return true; 
    else 
        pa := PathAlgebraContainingElement(gens[1]);
        if not IsPathAlgebra(pa) then 
            Print("IsQuadraticIdeal: The entered list of elements is not a list of elements in a path algebra, but in a quotient of a path algebra.\n");
            return fail;
        fi;
        # for each generator g of I
        for g in gens do
            L := CoefficientsAndMagmaElements(g);
            terms_occurring:= L{[1,3..Length(L)-1]};
            if not ForAll(terms_occurring,x->LengthOfPath(x) = 2) then
                return false;
            fi;
        od;   #end: for g in gens
    fi;
    return true;
end
);  #end: IsQuadraticIdeal()

#############################################################################
## 
#O  QuadraticPerpOfPathAlgebraIdeal( <gens> ) 
##
##  This function returns, when  <gens>  generated an ideal called  I and if 
##  all is well, Iperp inside the opposite path algebra in which  I  is an
##  ideal, and kQ_op. In this way the Koszul dual of kQ/I is KQ_op/Iper.
##  Otherwise, fail with error message (to STDERR) is returned.
##
##  The idea:
##  * make sure that inputs are okay and that I is quadratic;
##  * get an ordered k-basis for kQ_2;
##  * write the generators of I in terms of this basis;
##  * solve a linear system to find Iperp;
##  * translate the basis elements of Iperp to kQ_op;
##  * form the quotient.
##
InstallMethod( QuadraticPerpOfPathAlgebraIdeal, 
    "for a list of elements in a path algebra",
    [ IsHomogeneousList ], 0,
    function( gens ) 

    local Q, Q_op, pa, K, pa_op, p2, p2_op, b2, b2_op, p, A, 
          g, terms, gvec, pos, i, I, gb, b, t, 
          basisOfNullspace, basisOfIperp;

    if (not IsQuadraticIdeal(gens)) then
        Print("BasisOfIperp(): ideal is not quadratic\n");
        return fail;
    else 
        pa := PathAlgebraContainingElement(gens[1]);
        if (not IsPathAlgebra(pa)) then
            Print("BasisOfIperp(): arg is not a list of path algebra elements\n");
            return fail;
        else 
            K := LeftActingDomain(pa);
            Q := QuiverOfPathAlgebra(pa);
            Q_op := OppositeQuiver(Q);
            pa_op := PathAlgebra(LeftActingDomain(pa),Q_op);

            if (Length(gens) = 0) then
                Print ("KoszulDualOfPathAlgebraByIdeal(): empty generator list");
                return NthPowerOfArrowIdeal(pa_op,2);
            else 
                p2 := PathsOfLengthTwo(QuiverOfPathAlgebra(pa));
                
    ## It would seem that all this find-a-basis-of-Iperp stuff could
    ## be factored out, but strange things happen: even if you pass
    ## pa into such a function, the vectors you get back aren't in
    ## pa_op.

    # b2 will hold a list of the generators of kQ_2 in the order
    # determined by PathsOfLengthTwo()
                b2 := [];
                for p in p2 do
                    Add(b2,One(pa)*p);
                od;
                
    # compile the coefficient matrix A the rows of which are the
    # coefficients of the generators of I with respect to the
    # ordered k-basis of kQ_2
                A := [];
                for g in gens do
                    
        # get the terms of g
                    
                    terms := CoefficientsAndMagmaElements(g);
                    
        # Get the coefficients of the terms of g into the list gvec.
        #  --Note.  GAP appears to automatically eliminate terms with zero
        #    coefficient from its representation, so we need not worry
        #    about those.
                    gvec := ListWithIdenticalEntries(Length(b2),Zero(K));
                    for t in terms{[1,3..Length(terms)-1]} do
                        pos := Position(p2,t);
                        gvec[pos] := terms[Position(terms,t)+1];;
                    od;
                    Add(A,gvec);                 
                od; #end: for g in gens
            
                basisOfNullspace := NullspaceMat(TransposedMat(A));
            
        ## translate the basis vectors to kQ_2^op
        ##  --this all still works fine if basisOfNullspace is empty
        # first get a basis for kQ_2^op
        #
        # get the paths of length two in kQ2_op
        #  --thanks to opposite.gi!
                p2_op := [];
                for p in p2 do
                    t := List(WalkOfPath(p),String);
                    b := List(Reversed(t), function(n) return Concatenation(n,"_op"); end);
                    Add(p2_op,Q_op.(b[1])*Q_op.(b[2]));
                od;
            
        # b2_op will hold an ordered basis for kQ2_op
                b2_op := [];
                for b in p2_op do
                    Add(b2_op,One(pa_op)*b);
                od;

        # compile the wanted basis
                basisOfIperp := [];
                for b in basisOfNullspace do
                    t := Zero(pa_op);
                    for i in [1..Length(b)] do
                        t := t + b[i]*b2_op[i];
                    od;
                    Add(basisOfIperp,t);
                od;

                return [pa_op,basisOfIperp];
            fi;
        fi;
    fi;
end
); #end: QuadraticPerpOfPathAlgebraIdeal()
