# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)



##########################################################################
# This is an implementation of standard GAP operation for algebra (FLMLOR)
# TwoSidedIdealByGenerators (synonym IdealByGenerators)
# which is called by standard GAP global function TwoSidedIdeal
# (synonym Ideal). This is exactly the same code as original
# method + added SetIsIdealInPathAlgebra(I, true);
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
      SetIsIdealInPathAlgebra(I, true);
      
      CheckForHandlingByNiceBasis( lad, gens, I, false );
      return I;
    end 
); #TwoSidedIdealByGenerators


# In order to include trivial ideals, i.e. obtained
# by invoking Ideal(A, [])
InstallTrueMethod( IsIdealInPathAlgebra, 
                      IsRing and IsTrivial and HasLeftActingRingOfIdeal);





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
		
		A := LeftActingRingOfIdeal(I);
    
    if HasGroebnerBasisOfIdeal(I) then
      GB := GroebnerBasisOfIdeal(I);
    else
      rels := GeneratorsOfIdeal(I);  
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