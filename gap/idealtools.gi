# GAP Implementation
# This contains various tools for ideals in path algebras
# (created A. Mroz, 07.06.2012)

####################################################################################
# Membership test for an ideal in path algebra using GBNP Groebner Bases machinery
# Arguments: elt, I
# Returns: true if elt is a member of an ideal I, false otherwise
# Description: Does not check if I is in arrows ideal
	
InstallMethod(\in,
	"for a path algebra element and an ideal - membership test using Groebner Bases",
	IsElmsColls,
	[ IsRingElement, IsFLMLOR], 0,
	function( elt, I )
		local GB, pGB, A, rels, q, creps, parels, numv;
		
		A := LeftActingRingOfIdeal(I);
		rels := GeneratorsOfIdeal(I);
		
		# The following lines were copied from GBNPGroebnerBasisNC
		# This is the same code as for GBNPGroebnerBasisNC(rels, A), 
		# but without Prints
		#########################################################
		creps := [];
		q := QuiverOfPathAlgebra(A);   
		rels := MakeUniform(rels);
		numv := NumberOfVertices(q); 
		if numv > 1 then
			creps := QPA_Path2Cohen(rels);
			parels := QPA_RelationsForPathAlgebra(A);
			Append(creps,parels);
		else
			creps := QPA_Path2CohenFree(rels);
		fi;
		GB := SGrobner(creps);
		pGB :=  QPA_Cohen2Path(GB,A);
		pGB := Filtered(pGB, x -> x <> Zero(A));
		########################################################
		
		pGB := GroebnerBasis(I, pGB);
		pGB := CompletelyReduceGroebnerBasis(pGB);
		
		return CompletelyReduce(pGB, elt) = Zero(A);

	end 
);	# \in