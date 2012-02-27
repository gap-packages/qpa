# GAP Declarations
# This file was generated from
# $Id: groebner.gd,v 1.2 2012/02/27 12:26:34 sunnyquiver Exp $
DeclareCategory( "IsGroebnerBasis", IsCollection );

DeclareRepresentation(
    "IsGroebnerBasisDefaultRep",
    IsComponentObjectRep and IsAttributeStoringRep,
    ["ideal", "relations", "staticDict"]);

DeclareInfoClass("InfoGroebnerBasis");

DeclareOperation("GroebnerBasis", [IsFLMLOR, IsCollection]);
DeclareOperation("CompletelyReduce", [IsGroebnerBasis, IsRingElement]);
DeclareOperation("CompletelyReduceGroebnerBasis", 
                 [IsGroebnerBasis]);

DeclareOperation("TipReduce", [IsGroebnerBasis, IsRingElement]);
DeclareOperation("TipReduceGroebnerBasis", [IsGroebnerBasis]);

DeclareProperty("IsTipReducedGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsHomogeneousGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsCompleteGroebnerBasis", IsGroebnerBasis);
DeclareProperty("IsCompletelyReducedGroebnerBasis", IsGroebnerBasis);

DeclareRepresentation( "IsGroebnerBasisIteratorRep",
    IsComponentObjectRep,
    ["relations", "position"] );

DeclareOperation( "Nontips", [IsCompleteGroebnerBasis and IsGroebnerBasisDefaultRep] );
DeclareOperation( "AdmitsFinitelyManyNontips", [IsCompleteGroebnerBasis] );
DeclareOperation( "NontipSize", [IsCompleteGroebnerBasis] );
DeclareOperation( "IsPrefixOfTipInTipIdeal", 
    [IsCompleteGroebnerBasis, IsRingElement] );

DeclareCategory( "IsRightGroebnerBasis", IsGroebnerBasis );

DeclareOperation( "RightGroebnerBasis", [IsRing, IsCollection] );
