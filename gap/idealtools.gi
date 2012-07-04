# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)



##########################################################################
# This is an implementation of standard GAP operation for algebra (FLMLOR)
# TwoSidedIdealByGenerators (synonym IdealByGenerators)
# which is called by standard GAP global function TwoSidedIdeal
# (synonym Ideal). This is exactly the same code as original
# method + added SetIsIdealInPathAlgebra(I, true); and the trivial case.
# Now we can distinguish the ideals in path algebras.

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


# In order to include trivial ideals, i.e. obtained
# by invoking Ideal(A, [])
InstallTrueMethod( IsIdealInPathAlgebra, 
                      IsFLMLOR and IsTrivial and HasLeftActingRingOfIdeal);





####################################################################################
# Membership test for an ideal in path algebra using GBNP Groebner Bases machinery
# Arguments: <elt>, <I>
# Returns: true if <elt> is a member of an ideal <I>, false otherwise
# Description: For the efficiency reasons, it computes Groebner basis
#   for <I> only if it has not been computed. Similarly, it performs
#   CompletelyReduceGroebnerBasis only if it has not been reduced yet.
#   The method can change the existing Groebner basis (cf. the code).
# Remark: It computes Groebner basis only in case <I> is in the arrow 
#   ideal. (There are bugs in GBNPGroebnerBasisNC!)
	
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



######################################################################
# Property IsAdmissibleIdeal for an ideal in path algebra
# Argument: <I> - an ideal
# Returns: true, if <I> is an admissible ideal, i.e. I\subset R^2
#   and R^n \subset I, for some n (cf. the code), where R is the arrow ideal.
# Note: the second condition is checked by verifying if A/I
#  is a finite dimensional algebra. This uses Groebner bases machinery,
#  which sometimes can cause an infinite loop or another bugs!


InstallMethod( IsAdmissibleIdeal,
  "for an ideal in a path algebra",
  true,
  [IsIdealInPathAlgebra],
  0,
  function(I)
    local rels, i, Monomials, A;

    if HasIsTrivial(I) and IsTrivial(I) then
      return IsFiniteDimensional(LeftActingRingOfIdeal(I)); # <=> (arrow ideal)^n \subset I, for some n
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
    
    A := LeftActingRingOfIdeal(I);
  
    return IsFiniteDimensional(A/I); # <=> (arrow ideal)^n \subset I, for some n

  end
); # IsAdmissibleIdeal




######################################################################
# Property IsMonomialIdeal for an ideal in path algebra
# Argument: <I> - an ideal
# Returns: true, if <I> is a monomial ideal, 
#   i.e. <I> is generated by a set of monomials (= "zero-relations").
# Note:  This uses Groebner bases machinery
#  (which sometimes can cause an infinite loop or another bugs!).
#  It uses the observation: <I> is a monomial ideal <=> Groebner basis 
#  of <I> is a set of monomials.
#  It computes G.b. only in case it has not been computed yet and
#  usual generators of <I> are not monomials.


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
