# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gi,v 1.6 2012/08/01 16:01:10 sunnyquiver Exp $
InstallImmediateMethod( IsFiniteTypeAlgebra, IsNakayamaAlgebra, 0, A -> true);
        
#######################################################################
##
#O  NakayamaAlgebra( K, <admiss_seq> )
##
##  Given a field  <K>  and an admissible sequence  <admiss_seq>  this 
##  function constructs a Nakayama algebra with this data. 
##
InstallMethod ( NakayamaAlgebra,
    "for a field and an admissible sequence",
    [IsField, IsList], 0,
    function( K, admiss_seq )
    
    local num_vert, # number of vertices in the quiver we are creating.
          test,     # the result of the test of being an admissible seq.
          i, j,     # integer variables for for-loops.
          q, r,     # for writing a number n: n = q*num_vert + r.
          Q, KQ,    # the quiver and the path algebra for the 
          # Nakayama algebra we are creating.
          rels,     # relations of the Nakayama algebra we are creating. 
          new_rel,  # partial relation. 
          arrows,   # arrows of the quiver we are creating.
          cycle,    # if we are in the A-tilde case and the relations.
                    # are long, this is used to store the cycles involved
                    # in the relations. 
          I,        # the ideal of relations.
          gb, gbb,  # Groebner basis for the ideal I.
          A;
#
# Testing if we have an admissible sequence
#
    num_vert := Length(admiss_seq);
    test := true;
    i := 1;
    while ( test and ( i < num_vert ) ) do
        test := ( ( admiss_seq[i+1] >= admiss_seq[i] - 1 ) and
                  ( admiss_seq[i] - 1 >= 1 ) );  
        i := i + 1;
    od;
    if test then 
        test := ( admiss_seq[1] >= admiss_seq[num_vert] - 1 );
    fi;
#
# If an admissible sequence, the case admiss_seq[num_vert] = 1, is the
# A_n case, and the other case is A_n-tilde case. 
#
    if test then 
        if ( admiss_seq[num_vert] = 1 ) then
#
# This is the A_n case.
#
            arrows := [];
            for i in [1..num_vert - 1] do
                Add(arrows,[i,i+1]);
            od;
            Q := Quiver(num_vert,arrows);
            KQ := PathAlgebra(K,Q); 
            arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert-1]}; 
            rels := [];
            for i in [1..num_vert - 1] do
                if (( admiss_seq[i] < num_vert - i + 1 ) and 
                    ( admiss_seq[i] - 1 <> admiss_seq[i+1] )) then
                    new_rel := One(KQ);
                    for j in [i..i+admiss_seq[i] - 1] do
                        new_rel := new_rel*arrows[j];
                    od;
                    Add(rels,new_rel); 
                fi;
            od;
#
# end of A_n case
#
        else
#
# This is the A_n-tilde case.
            arrows := [];
            for i in [1..num_vert - 1] do
                Add(arrows,[i,i+1]);
            od;
            Add(arrows,[num_vert,1]); 
            Q := Quiver(num_vert,arrows);
            KQ := PathAlgebra(K,Q); 
            arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert]}; 
            rels := [];
#
# Relations starting in the vertices 1..num_vert - 1. 
#
            for i in [1..num_vert - 1] do
                if ( admiss_seq[i] - 1 <> admiss_seq[i+1] ) then
                    new_rel := One(KQ);
                    r := admiss_seq[i] mod num_vert;
                    q := (admiss_seq[i] - r)/num_vert;
                    cycle := One(KQ);
                    for j in [i..num_vert] do
                        cycle := cycle*arrows[j];
                    od;
                    for j in [1..i-1] do
                        cycle := cycle*arrows[j];
                    od;
                    for j in [1..q] do 
                        new_rel := new_rel*cycle;
                    od;
                    if ( i + r - 1 <= num_vert ) then
                        for j in [i..i+r - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                    else
                        for j in [i..num_vert] do
                            new_rel := new_rel*arrows[j];
                        od;
                        for j in [1..r - (num_vert - i) - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                    fi;
                    Add(rels,new_rel); 
                fi;
            od;
#
# Relation starting in the vertex num_vert. 
#
            if ( admiss_seq[num_vert] - 1 <> admiss_seq[1] ) then
                new_rel := One(KQ);
                r := admiss_seq[num_vert] mod num_vert;
                q := (admiss_seq[num_vert] - r)/num_vert;
                cycle := One(KQ)*arrows[num_vert];
                for j in [1..num_vert-1] do
                    cycle := cycle*arrows[j];
                od;
                for j in [1..q] do 
                    new_rel := new_rel*cycle;
                od;
                if ( r <> 0 ) then
                    new_rel := new_rel*arrows[num_vert];
                    for j in [1..r - 1] do
                        new_rel := new_rel*arrows[j];
                    od;
                fi;
                Add(rels,new_rel);
            fi;
#
# End of A_n-tilde case.
#
        fi;
        if Length(rels) = 0 then 
            SetFilterObj(KQ, IsNakayamaAlgebra and IsPathAlgebra );
            
            return KQ;
        else
            I := Ideal(KQ,rels);
            gb := GBNPGroebnerBasis(rels,KQ);
            gbb := GroebnerBasis(I,gb);
            A := KQ/I;
            SetFilterObj(A, IsNakayamaAlgebra and IsAdmissibleQuotientOfPathAlgebra );               
            return A;
        fi; 
    else
#
# A valid admissible sequence was not entered. 
# 
        Print("The admissible sequence is not a valid one.\n");
        return admiss_seq;
    fi; 
    
end
);

InstallOtherMethod ( NakayamaAlgebra,
    "for an admissible sequence and a field",
    [IsList, IsField], 0,
        function( admiss_seq , K)
    
    return NakayamaAlgebra(K,admiss_seq);
end
);

#######################################################################
##
#O  CanonicalAlgebra( <K>, <weights>, <relcoeff> )
##
##  Given a field  K, a sequence of  weights  and a sequence of 
##  coefficients for the relations, this function constructs the 
##  canonical algebra with this data. If the number of weights is
##  two, then the number of coefficients for the relations must be
##  zero, that is, represented by an empty list, [].
##
InstallMethod ( CanonicalAlgebra,
    "for a field, list of lengths of arms, and coefficients for the relations", 
    [IsField, IsList, IsList], 0,
    function( K, weights, relcoeff)
    
    local n, vertices, i, j, arrows, Q, KQ, num_vert, generators, arms,
          start, temp, relations, I, gb, gbb, A; 

        #
        # Checking the input.
        #
    if Length(weights) <> Length(relcoeff) + 2 then
        Error("number of arms in the quiver don't match the number of coeffiecients for the relations in the quiver,");
    fi;
    if not ForAll( weights, x -> x >= 2) then
        Error("not all weights are greater or equal to 2,");
    fi;
    if not ForAll( relcoeff, x -> x in K ) then
        Error("not all coefficients for the relations are in the entered field,");
    fi;
    if not ForAll( relcoeff, x -> x <> Zero(K) ) then
        Error("not all coefficients for the relations are non-zero,");
    fi;
        #
        # When input is good, start constructing the algebra. 
        # First construct the vertices.
        #
    n := Length(weights);
    vertices := [];
    Add(vertices,"v");
    for i in [1..n] do
        for j in [2..weights[i]] do
            Add(vertices,Concatenation("v",String(i),String(j)));
        od;
    od;
    Add(vertices,"w");
        #
        # Next construct the arrows.
        #
    arrows := [];
    for i in [1..n] do
        Add(arrows,["v",Concatenation("v",String(i),String(2)),Concatenation("a",String(i))]);
        for j in [2..weights[i] - 1] do
            Add(arrows,[Concatenation("v",String(i),String(j)), 
                    Concatenation("v",String(i),String(j+1)),
                    Concatenation("a",String(i),String(j))]);
        od;
        Add(arrows,[Concatenation("v",String(i),String(weights[i])),"w",Concatenation("b",String(i))]);
    od;
        # 
        # Constructing the quiver and path algebra.
        #
    Q := Quiver(vertices,arrows);
    KQ := PathAlgebra(K,Q);
        #
        # If n = 2, it is a path algebra, otherwise construct relations 
        # and quotient of the path algbra just constructed.
        #
    if n = 2 then 
        SetFilterObj(KQ, IsCanonicalAlgebra and IsPathAlgebra );

        return KQ;
    else
        # 
        # First finding all the arrows as elements in the path algebra.
        #
        num_vert := NumberOfVertices(Q);
        generators := GeneratorsOfAlgebra(KQ){[num_vert + 1..num_vert + NumberOfArrows(Q)]};
        #
        # Computing the products of all arrows in an arm of the quiver as 
        # elements in the path algebra. 
        # 
        arms := [];
        start := 1;
        for i in [1..n] do
            temp := One(KQ);
            for j in [start..start + weights[i] - 1] do
                temp := temp*generators[j];
            od;
            start := start + weights[i];
            Add(arms,temp);
        od;
        relations := [];
        # 
        # Constructing the relations and the quotient.
        #
        for i in [3..n] do
            Add(relations,arms[i] - arms[2] + relcoeff[i - 2]*arms[1]);
        od;
        I := Ideal(KQ,relations);
        gb := GBNPGroebnerBasis(relations,KQ);
        gbb := GroebnerBasis(I,gb);
        A := KQ/I;
        SetFilterObj(A, IsCanonicalAlgebra and IsAdmissibleQuotientOfPathAlgebra );
        return A;
    fi;
end
);

#######################################################################
##
#O  CanonicalAlgebra( <K>, <weights>)
##
##  CanonicalAlgebra with only two arguments, the field and the 
##  two weights.
##
InstallOtherMethod ( CanonicalAlgebra,
    "for a field and a list of lengths of two arms", 
    [IsField, IsList ], 0,
    function( K, weights);

    if Length(weights) <> 2 then 
        Error("the list of weights is different from two, need to enter coefficients for the relations then,");
    fi;
    
    return CanonicalAlgebra(K,weights,[]);
end
);

#######################################################################
##
#O  KroneckerAlgebra( <K>, <n> )
##
##  Given a field  K and a positive integer  n, this function constructs  
##  the Kronecker algebra with  n  arrows. 
##
InstallMethod ( KroneckerAlgebra,
    "for a field and a positive integer", 
    [IsField, IS_INT], 0,
    function( K, n)
    
    local arrows, i, Q, KQ;

    if n < 2 then
        Error("the number of arrows in the Kronecker quiver must be greater or equal to 2,");
    fi;
    arrows := [];
    for i in [1..n] do
        Add(arrows,[1,2,Concatenation("a",String(i))]);
    od;
    Q := Quiver(2,arrows);
    KQ := PathAlgebra(K,Q);
    SetFilterObj(KQ,IsKroneckerAlgebra and IsPathAlgebra);
    return KQ;
end
);


#################################################################
##
#P  IsSpecialBiserialQuiver( <Q> ) 
## 
##  It tests if every vertex in quiver <Q> is a source (resp. target) 
##  of at most 2 arrows.
##
##  NOTE: e.g. path algebra of one loop IS NOT special biserial, but
##        one loop IS special biserial quiver.
##
InstallMethod( IsSpecialBiserialQuiver,
    "for quivers",
    true,
    [ IsQuiver ], 0,
    function ( Q )
    
    local vertex_list, v;
    
    vertex_list := VerticesOfQuiver(Q);
    for v in vertex_list do
        if OutDegreeOfVertex(v) > 2 then
            return false;
        fi;
        if InDegreeOfVertex(v) > 2 then
            return false;
        fi;
    od;
    
    return true;
end
); # IsSpecialBiserialQuiver


#################################################################
##
#P  IsSpecialBiserialAlgebra( <A> ) 
##  <A> = a path algebra
##
##  It tests if a path algebra <A> is an algebra of  a quiver Q which is
##  (IsSpecialBiserialQuiver and IsAcyclicQuiver) and the ideal 0 satisfies
##  the "special biserial" conditions (cf. comment below).
##  
##  NOTE: e.g. path algebra of one loop IS NOT special biserial, but
##        one loop IS special biserial quiver.
##		
InstallMethod( IsSpecialBiserialAlgebra,
    "for path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    
    local Q, alpha; 
    
    Q := QuiverOfPathAlgebra(A);
    if not IsSpecialBiserialQuiver(Q) then
        return false;
    fi;
    
    for alpha in ArrowsOfQuiver(Q) do
        if OutDegreeOfVertex(TargetOfPath(alpha)) > 1 then
            return false;
        fi;          
        if InDegreeOfVertex(SourceOfPath(alpha)) > 1 then  
            return false;
        fi;
    od;
    
    return IsAcyclicQuiver(Q);  # <=> 0 is an admissible ideal
end
); # IsSpecialBiserialAlgebra for IsPathAlgebra


########################################################################
##
#P  IsSpecialBiserialAlgebra( <A> ) 
##  <A> = a quotient of a path algebra
##
##  It tests if an original path algebra is an algebra of a quiver Q which is
##  IsSpecialBiserialQuiver, I is an admissible ideal and I satisfies 
##  the "special biserial" conditions, i.e.:
##  for any arrow a there exists at most one arrow b such that ab\notin I
##  and there exists at most one arrow c such that ca\notin I.
##
  
InstallMethod( IsSpecialBiserialAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local Q, PA, I, v, alpha, beta, not_in_ideal;
    
    PA := OriginalPathAlgebra(A);
    Q := QuiverOfPathAlgebra(PA);
    if not IsSpecialBiserialQuiver(Q) then
        return false;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    
    if not IsAdmissibleIdeal(I) then
        return false;
    fi;
    
    for alpha in ArrowsOfQuiver(Q) do
        not_in_ideal := 0;
        for beta in OutgoingArrowsOfVertex(TargetOfPath(alpha)) do
            if not ElementOfPathAlgebra(PA, alpha*beta) in I then
                not_in_ideal := not_in_ideal + 1;
            fi;
        od;
        if not_in_ideal > 1 then
            return false;
        fi;          
        not_in_ideal := 0;
        for beta in IncomingArrowsOfVertex(SourceOfPath(alpha)) do
            if not ElementOfPathAlgebra(PA, beta*alpha) in I then
                not_in_ideal := not_in_ideal + 1;
            fi;
        od;
        if not_in_ideal > 1 then
            return false;
        fi;
    od;
    
    return true;
end
); # IsSpecialBiserialAlgebra for IsQuotientOfPathAlgebra


#################################################################
##
#P  IsStringAlgebra( <A> )
##  <A> = a path algebra
##
##  Note: A path algebra is a string algebra <=> it is a special biserial algebra
##
InstallMethod( IsStringAlgebra,
    "for quotients of path algebras",
    true,
    [ IsPathAlgebra ], 0,
    function ( A )
    return IsSpecialBiserialAlgebra(A);
end
); # IsStringAlgebra for IsPathAlgebra                      


#################################################################
##
#P  IsStringAlgebra( <A> )
##  <A> = a quotient of a path algebra
##
##  Note: kQ/I is a string algebra <=> kQ/I is a special biserial algebra
##                                     and I is a monomial ideal.                                         
InstallMethod( IsStringAlgebra,
    "for quotients of path algebras",
    true,
    [ IsQuotientOfPathAlgebra ], 0,
    function ( A )
    
    local I;
    
    if not IsSpecialBiserialAlgebra(A) then
        return false;
    fi;
    
    I := ElementsFamily(FamilyObj(A))!.ideal;
    return IsMonomialIdeal(I);
end
); # IsStringAlgebra for IsQuotientOfPathAlgebra 