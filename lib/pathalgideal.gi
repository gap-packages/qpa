# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)



#############################################################################
##
#O  TwoSidedIdealByGenerators(<A>,<gens>)
##
##  This function creates and returns an ideal I (an object with the property
##  IsIdealInPathAlgebra) in path algebra <A> generated by a collection
##  of elements <gens>.
##  This is an implementation of standard GAP operation for algebra 
##  (FLMLOR) TwoSidedIdealByGenerators (synonym IdealByGenerators)
##  which is called by standard GAP global function TwoSidedIdeal
##  (synonym Ideal). This is exactly the same code as original
##  method + added SetIsIdealInPathAlgebra(I, true); and the trivial case.
##  Now we can distinguish the ideals in path algebras.
##
InstallMethod( TwoSidedIdealByGenerators,
    "for path algebra and collection",
    IsIdenticalObj,
    [ IsPathAlgebra, IsCollection ], 0,
    function( A, gens )
      local I, lad;
      I:= Objectify( NewType( FamilyObj( A ),
                                  IsFLMLOR
                              and IsAttributeStoringRep ),
                    rec() );
      lad:= LeftActingDomain( A );
      SetLeftActingDomain( I, lad );
      SetGeneratorsOfTwoSidedIdeal( I, gens );
      SetLeftActingRingOfIdeal( I, A );
      SetRightActingRingOfIdeal( I, A );
      
      # The only difference
      if ForAll(gens, x -> IsZero(x)) then # to avoid bugs for "boundary" situations for \in etc.!
        SetIsTrivial(I, true);
      fi;
      SetIsIdealInPathAlgebra(I, true);
      
      CheckForHandlingByNiceBasis( lad, gens, I, false );
      return I;
    end 
); #TwoSidedIdealByGenerators

#############################################################################
##
#P  IsIdealInPathAlgebra(<I>)
##
##  This "true method" for the property IsIdealInPathAlgebra returns true 
##  if <I> is a trivial ideal in path algebra, e.g. obtained by invoking 
##  Ideal(A, []).
##
InstallTrueMethod( IsIdealInPathAlgebra, 
                      IsFLMLOR and IsTrivial and HasLeftActingRingOfIdeal);





####################################################################################
##
#O  \in(<elt>, <I>)
##
##  This function performs a membership test for an ideal in path algebra using 
##  Groebner Bases machinery.
##  It returns true if <elt> is a member of an ideal <I>, false otherwise.
##  For the efficiency reasons, it computes Groebner basis
##  for <I> only if it has not been computed. Similarly, it performs
##  CompletelyReduceGroebnerBasis only if it has not been reduced yet.
##  The method can change the existing Groebner basis (cf. the code).
##  It computes Groebner basis only in case <I> is in the arrow 
##  ideal.
##	
InstallMethod(\in,
	"for a path algebra element and an ideal - membership test using Groebner bases",
	IsElmsColls,
	[ IsRingElement, IsIdealInPathAlgebra], 0,
	function( elt, I )
		local GB, rels, A;
		
    if IsZero(elt) then 
      return true;
    fi;
    if HasIsTrivial(I) and IsTrivial(I) then
      return false;
    fi;
    
    if HasLeftActingRingOfIdeal(I) then
      A := LeftActingRingOfIdeal(I);
    else 
      TryNextMethod();
    fi;
    
    
    if HasGroebnerBasisOfIdeal(I) then
      GB := GroebnerBasisOfIdeal(I);
    else
      if HasGeneratorsOfIdeal(I) then
        rels := GeneratorsOfIdeal(I);
      else
        TryNextMethod();
      fi;      
      GB := GBNPGroebnerBasis(rels, A);
      GB := GroebnerBasis(I, GB);
    fi;
		
    if (not HasIsCompletelyReducedGroebnerBasis(GB)) or
       (HasIsCompletelyReducedGroebnerBasis(GB) and 
        not IsCompletelyReducedGroebnerBasis(GB) ) then
      # This modifies the basis GroebnerBasisOfIdeal(I),
      # i.e. it becomes reduced since now.
      # Strange behaviour, don't know how to avoid it.    
      GB := CompletelyReduceGroebnerBasis(GB);
    fi; 
		
		return CompletelyReduce(GB, elt) = Zero(A);

	end 
);	# \in

############################################################################
##
#P  IsAdmissibleIdeal(<I>)
##  
##  This function returns true if <I> is an admissible ideal, i.e. 
##  I\subset R^2 and R^n \subset I, for some n (cf. the code), where R 
##  is the arrow ideal.
##  Note: the second condition is checked by first computing to which power
##  n  of the radical of  A/I  is zero, then if  R^n  is contained in  I. 
##
InstallMethod( IsAdmissibleIdeal,
  "for an ideal in a path algebra",
  true,
  [IsIdealInPathAlgebra],
  0,
  function(I)
    local gb, rels, i, Monomials, A, B, arrows, J, dim, powerofJ;

    if HasIsTrivial(I) and IsTrivial(I) then
      return IsFiniteDimensional(LeftActingRingOfIdeal(I)); # <=> (arrow ideal)^n \subset I, for some n
    fi;
  
    if HasGroebnerBasisOfIdeal(I) then
      gb := GroebnerBasisOfIdeal(I);
    else
      rels := GeneratorsOfIdeal(I);     
      gb := GBNPGroebnerBasis(rels, LeftActingRingOfIdeal(I));
      gb := GroebnerBasis(I, gb);
    fi;
    
    # returns the list of monomials appearing in elt
    Monomials := function(elt)
      local terms;
      terms := CoefficientsAndMagmaElements(elt);
      return terms{[1,3..Length(terms)-1]};
    end;
    
    # Check if all monomials in all relations belong to (arrow ideal)^2
    rels := GeneratorsOfIdeal(I);
    if not ForAll(rels, r -> ForAll(Monomials(r), m -> (LengthOfPath(m) >= 2)) ) then
      return false;
    fi;
    #
    # Let  J  be the ideal generated by the arrows in B.  Start computing  J^2,  J^4,
    # J^8, and so on.  Stop when this chain stablizes, that is, Dimension(J^{2n}) =
    # Dimension(J^{2n+2}).  If this dimension is different from zero, the ideal is not
    # admissible.  If this dimension is zero, then the ideal is admissible.
    # 
    A := LeftActingRingOfIdeal(I);
    B := A/I;
    if IsFiniteDimensional(B) then
        arrows := List(ArrowsOfQuiver(QuiverOfPathAlgebra(B)), a -> One(B)*a);
        J := Ideal(B, arrows);
        dim := Dimension(J);
        powerofJ := ProductSpace( J, J); 
        while dim <> Dimension(powerofJ) do
            dim := Dimension(powerofJ);
            powerofJ := ProductSpace(powerofJ, powerofJ);
        od;
        return dim = 0;
    else
        return false;
    fi;
end
); # IsAdmissibleIdeal

######################################################################
##
#P  IsMonomialIdeal(<I>)
## 
##  This function returns true if <I> is a monomial ideal, 
##  i.e. <I> is generated by a set of monomials (= "zero-relations").
##  Note:  This uses Groebner bases machinery
##  (which sometimes can cause an infinite loop or another bugs!).
##  It uses an observation: <I> is a monomial ideal <=> Groebner basis 
##  of <I> is a set of monomials.
##  It computes G.b. only in case it has not been computed yet and
##  usual generators of <I> are not monomials.
##
InstallMethod( IsMonomialIdeal,
  "for an ideal in a path algebra",
  true,
  [IsIdealInPathAlgebra],
  0,
  function(I)
    local GB, A, rels;

    if HasIsTrivial(I) and IsTrivial(I) then
      return true;
    fi;
    
    # First check if usual generators are just monomials
    if HasGeneratorsOfIdeal(I) then
      rels := GeneratorsOfIdeal(I);  
      if ForAll(rels, r -> (Length(CoefficientsAndMagmaElements(r)) = 2) ) then
        return true;
      fi;
    else return fail;    
    fi;
    
    # Now we have to check if Groebner basis is a set of monomials
    # Compute Groebner basis if necessary
    if HasGroebnerBasisOfIdeal(I) then
      GB := GroebnerBasisOfIdeal(I);
    else
      GB := GBNPGroebnerBasis(rels, LeftActingRingOfIdeal(I));
      GB := GroebnerBasis(I, GB);
    fi;
    
    # Checks if GB is a list of monomials
    return ForAll(GB, r -> (Length(CoefficientsAndMagmaElements(r)) = 2) );

  end
); # IsMonomialIdeal
  
#######################################################################
##
#O  ProductOfIdeals( <I>, <J> )
##
##  Given two ideals  I  and  J  in a path algebra  A, this function  
##  computes the product  I*J  of the ideal  I  and  J, if the ideal 
##  if  J  admits finitely many nontips in  A. 
##
InstallMethod ( ProductOfIdeals, 
    "for two IdealInPathAlgebra",
    true,
    [ IsIdealInPathAlgebra, IsIdealInPathAlgebra ], 
    0,
    function( I, J )

    local A, rgb, gens_I, gens_IJ, a, b;
#
#   Checking the input
#
    A := LeftActingRingOfIdeal(I);
    if LeftActingRingOfIdeal(J) <> A then
        Error("the ideals are not ideals in the same ring,");
    fi;
#
#   If the ideal  J  admits finitely many nontips in A (that is,
#   if A/J is finite dimensional, then do the computations.
#
    if AdmitsFinitelyManyNontips(GroebnerBasisOfIdeal(J)) then 
        rgb := RightGroebnerBasis(J);
        gens_I := GeneratorsOfTwoSidedIdeal(I);
        gens_IJ := [];
        for a in gens_I do
            for b in rgb do
                Add(gens_IJ, a*b);
            od;
        od;
        return Ideal(A,gens_IJ);
    else
        Error("the second argument does not admit finitely many nontips,");
    fi;
end
);

InstallMethod( IsAdmissibleQuotientOfPathAlgebra, 
    "for a quotient of a path algebra",
    [ IsQuotientOfPathAlgebra ], 0,
    function( A )
    
    local fam;
    
    fam := ElementsFamily(FamilyObj(A));
    
    return IsAdmissibleIdeal(fam!.ideal);
end
  );

#######################################################################
##
#O  MinimalGeneratingSetOfIdeal( <I> )
##
##  The argument of this function is an admissible ideal  <I>  in a 
##  path algebra  KQ  a field  K. 
##  
##  
InstallMethod( MinimalGeneratingSetOfIdeal, 
    "for an admissible ideal in a path algebra",
    [ IsAdmissibleIdeal ], 0,
    function( I )

    local generators, KQ, arrows, JIplusIJ, Aprime, fam, Ibar, B;
    #
    #  Finding a generating set for  JI + IJ. 
    #
    generators := GeneratorsOfTwoSidedIdeal(I);
    if Length(generators) = 0 then
        return [];
    fi;
    KQ := LeftActingRingOfIdeal(I);
    arrows := List(ArrowsOfQuiver(QuiverOfPathAlgebra(KQ)), a -> a*One(KQ));
    JIplusIJ := List(arrows, x -> Filtered(generators*x, y -> y <> Zero(y)));
    Append(JIplusIJ,List(arrows, x -> Filtered(x*generators, y -> y <> Zero(y))));
    JIplusIJ := Flat(JIplusIJ);
    #
    # Constructing the factor KQ/JIplusIJ, if necessary.
    #
    if Length( JIplusIJ ) = 0 then
        return BasisVectors( Basis( I ) );
    else
        Aprime := KQ/JIplusIJ;
        #
        # Finding the ideal  I/JIplusIJ  and a basis of this ideal. Pulling this 
        # basis back to KQ and return it.
        # 
        fam := ElementsFamily(FamilyObj(Aprime));
        Ibar := List(generators, x -> ElementOfQuotientOfPathAlgebra(fam, x, true)); 
        Ibar := Ideal(Aprime, Ibar);
        B := BasisVectors(Basis(Ibar));
    
        return List(B, x -> x![1]);
    fi;
end
  );

#######################################################################
##
#P  IsGentleAlgebra( <A> )
##
##  The argument of this function is a quiver algebra  <A>. The function
##  returns true is  <A>  is a gentle algebra, and false otherwise.
##  
InstallMethod( IsGentleAlgebra, 
    "for a quiver algebra",
    [ IsQuiverAlgebra ], 0,
    function( A )

    local fam, I, minimalgenerators, coeffandmagmaelements, Q, KQ, beta, test;
    
    if IsPathAlgebra(A) then 
        return IsSpecialBiserialAlgebra(A);
    fi;

    fam := ElementsFamily(FamilyObj(A));
    I := fam!.ideal;
    #
    # Checking if  A  is a special biserial algebra, and if it is 
    # quotient of a path algebra by an admissible ideal.
    #
    if not IsSpecialBiserialAlgebra(A) then
        return false;
    fi;
    if not IsAdmissibleQuotientOfPathAlgebra(A) then
        return false;
    fi;
    #
    # Checking if the ideal is generated by paths.
    #
    minimalgenerators := MinimalGeneratingSetOfIdeal(I);
    coeffandmagmaelements := List(minimalgenerators, x -> CoefficientsAndMagmaElements(x));
    if not ForAll(coeffandmagmaelements, x -> Length(x) = 2) then
        return false;
    fi;
    # 
    # Checking if the ideal is generated by paths of length  2.
    #    
    if not ForAll(coeffandmagmaelements, x -> LengthOfPath(x[1]) = 2) then
        return false;
    fi;
    #
    # Checking if for any arrow beta, there is at most one arrow gamma with beta*gamma in I, and
    # if for any arrow beta, there is at most one arrow gamma with gamma*beta in I.
    #
    Q := QuiverOfPathAlgebra(A); 
    KQ := OriginalPathAlgebra(A);
    for beta in ArrowsOfQuiver(Q) do
        test := Filtered(OutgoingArrowsOfVertex(TargetOfPath(beta)), gamma -> ElementOfPathAlgebra(KQ, beta*gamma) in I);
        if Length(test) > 1 then 
            return false;
        fi;
        test := Filtered(IncomingArrowsOfVertex(SourceOfPath(beta)), gamma -> ElementOfPathAlgebra(KQ, gamma*beta) in I);
        if Length(test) > 1 then 
            return false;
        fi;
    od;
    #
    # By now, all is good, return true.
    #
    return true;
end
  );
