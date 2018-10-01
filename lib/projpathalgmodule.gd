# GAP Declarations
# This file was generated from
# $Id: present.gd,v 1.3 2012/09/27 08:55:07 sunnyquiver Exp $
DeclareInfoClass( "InfoPathAlgebraModule" );

DeclareProperty( "IsPathAlgebraModule", IsAlgebraModule );
DeclareProperty( "IsVertexProjectiveModule", IsPathAlgebraModule );

DeclareCategory( "IsPathAlgebraVector", IsVector );
DeclareCategoryCollections( "IsPathAlgebraVector" );
DeclareCategoryFamily( "IsPathAlgebraVector" );
DeclareRepresentation("IsPathAlgebraVectorDefaultRep",
                      IsPositionalObjectRep, []);

DeclareOperation( "PathAlgebraVector",
    [IsPathAlgebraVectorFamily, IsHomogeneousList] );
DeclareOperation( "PathAlgebraVectorNC",
    [IsPathAlgebraVectorFamily, IsHomogeneousList, IsInt, IsBool] );
DeclareOperation( "Vectorize", [IsVectorCollection, IsHomogeneousList] );

DeclareSynonym("TargetVertex", TargetOfPath);
DeclareSynonym("SourceVertex", SourceOfPath);

DeclareAttribute( "LeadingComponent", IsPathAlgebraVector );
DeclareAttribute( "LeadingPosition", IsPathAlgebraVector );

DeclareOperation( "IsLeftDivisible",
    [IsPathAlgebraVector, IsPathAlgebraVector]);

DeclareCategory( "IsBasisOfPathAlgebraVectorSpace", IsBasis );
DeclareOperation( "RightProjectiveModule", [IsRing, IsObject] );
DeclareAttribute( "UniformGeneratorsOfModule", IsPathAlgebraModule );
DeclareAttribute( "RightGroebnerBasisOfModule", IsPathAlgebraModule );
DeclareOperation( "RightGroebnerBasisOfModule", [IsRing, IsHomogeneousList] );
DeclareCategory( "IsPathAlgebraModuleGroebnerBasis", IsObject );
DeclareProperty( "IsRightPathAlgebraModuleGroebnerBasis",
                 IsPathAlgebraModuleGroebnerBasis );
DeclareProperty( "IsLeftPathAlgebraModuleGroebnerBasis",
                 IsPathAlgebraModuleGroebnerBasis );
DeclareAttribute( "UnderlyingModule", IsPathAlgebraModuleGroebnerBasis );
DeclareRepresentation( "IsPathAlgebraModuleGroebnerBasisDefaultRep",
    IsAttributeStoringRep, ["staticDictionaries", "gbasisElems"] );
DeclareProperty( "IsFpPathAlgebraModule", IsAlgebraModule );

DeclareCategory( "IsFpPathAlgebraVector", IsPathAlgebraVector );
DeclareCategoryCollections( "IsFpPathAlgebraVector" );
DeclareCategoryFamily( "IsFpPathAlgebraVector" );

DeclareOperation( "LiftPathAlgebraModule", [IsPathAlgebraModule] );
DeclareOperation( "VertexProjectivePresentation", 
                  [ IsAlgebra, IsRingElementTable ] );
# DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );

DeclareOperation( "NewBasis",[ IsFreeLeftModule and IsPathAlgebraVectorCollection, IsList ] );
DeclareOperation( "NewBasis",[ IsFreeLeftModule and IsPathAlgebraModule, IsList ] );
#DeclareOperation( "NewBasis",[ IsFreeLeftModule  , IsList ]);

DeclareOperation( "BasisOfDomain", [IsFreeLeftModule and IsPathAlgebraVectorCollection]);
DeclareOperation( "CompletelyReduceGroebnerBasisForModule", [IsPathAlgebraModuleGroebnerBasis ] );
DeclareOperation( "ProjectivePathAlgebraPresentation", [IsPathAlgebraMatModule ] );