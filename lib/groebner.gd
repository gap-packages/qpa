# GAP Declarations
# This file was generated from
# $Id: groebner.gd,v 1.3 2012/08/01 16:01:10 sunnyquiver Exp $
DeclareCategory( "IsGroebnerBasis", IsCollection );

DeclareRepresentation(
    "IsGroebnerBasisDefaultRep",
    IsComponentObjectRep and IsAttributeStoringRep,
    ["ideal", "relations", "staticDict"]);

DeclareInfoClass("InfoGroebnerBasis");

DeclareOperation("GroebnerBasis", [IsFLMLOR, IsCollection]);
DeclareOperation("CompletelyReduce", [IsGroebnerBasis, IsRingElement]);
DeclareOperation("CompletelyReduceGroebnerBasis", [IsGroebnerBasis]);
DeclareOperation("TipReduce", [IsGroebnerBasis, IsRingElement]);
DeclareOperation("TipReduceGroebnerBasis", [IsGroebnerBasis]);
DeclareProperty("IsTipReducedGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsHomogeneousGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsCompleteGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsCompletelyReducedGroebnerBasis", IsGroebnerBasis);

DeclareRepresentation( "IsGroebnerBasisIteratorRep",
    IsComponentObjectRep, ["relations", "position"] );

DeclareOperation( "Nontips", [IsCompleteGroebnerBasis and IsGroebnerBasisDefaultRep] );
DeclareOperation( "AdmitsFinitelyManyNontips", [ IsCompleteGroebnerBasis ] );
DeclareAttribute( "NontipSize", IsCompleteGroebnerBasis );
DeclareOperation( "IsPrefixOfTipInTipIdeal", 
    [IsCompleteGroebnerBasis, IsRingElement] );

DeclareCategory( "IsRightGroebnerBasis", IsGroebnerBasis );

DeclareOperation( "RightGroebnerBasis", [IsRing, IsCollection] );
