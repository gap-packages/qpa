# GAP Declarations
# $Id: homomorphisms.gd,v 1.9 2012/02/27 12:26:34 sunnyquiver Exp $

DeclareCategory("IsPathAlgebraModuleHomomorphism", IsAdditiveElementWithZero and IsAdditiveElementWithInverse and IsGeneralMapping and RespectsAddition and RespectsZero and RespectsScalarMultiplication and IsTotal and IsSingleValued ); 
DeclareCategoryFamily(  "IsPathAlgebraModuleHomomorphism" );
DeclareCategoryCollections( "IsPathAlgebraModuleHomomorphism" );
DeclareRepresentation("IsPathAlgebraModuleHomomorphismRep", IsComponentObjectRep and IsAttributeStoringRep, ["maps"]);

DeclareOperation( "RightModuleHomOverAlgebra", [IsPathAlgebraModule, IsPathAlgebraModule, IsList] ); 
DeclareAttribute( "PathAlgebraOfMatModuleMap", IsPathAlgebraModuleHomomorphism );
DeclareOperation( "HomOverAlgebra", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "EndOverAlgebra", [IsPathAlgebraModule ] ); 
DeclareOperation( "NumberOfNonIsoDirSummands", [IsPathAlgebraModule ] ); 
DeclareOperation( "LiftIdempotents", [IsPathAlgebraModule ] ); 
DeclareAttribute( "KernelInclusion", IsPathAlgebraModuleHomomorphism ); 
DeclareAttribute( "ImageProjectionInclusion", IsPathAlgebraModuleHomomorphism ); 
DeclareAttribute( "ImageProjection", IsPathAlgebraModuleHomomorphism ); 
DeclareAttribute( "ImageInclusion", IsPathAlgebraModuleHomomorphism ); 
DeclareAttribute( "CoKernelProjection", IsPathAlgebraModuleHomomorphism ); 
DeclareAttribute( "KernelOfWhat", IsPathAlgebraModuleHomomorphism );
DeclareAttribute( "CoKernelOfWhat", IsPathAlgebraModuleHomomorphism );
DeclareAttribute( "ImageOfWhat", IsPathAlgebraModuleHomomorphism );
DeclareProperty( "IsIsomorphism", IsPathAlgebraModuleHomomorphism );
DeclareOperation( "SubRepresentation", [IsPathAlgebraModule, IsList]);
DeclareOperation( "SubRepresentationInclusion", [IsPathAlgebraModule, IsList]);
DeclareOperation( "RadicalOfModule", [IsPathAlgebraModule]);
DeclareOperation( "RadicalOfModuleInclusion", [IsPathAlgebraModule]);
DeclareOperation( "TopOfModule", [IsPathAlgebraModule]);
DeclareOperation( "TopOfModuleProjection", [IsPathAlgebraModule]);
DeclareOperation( "RightFacApproximation", [IsPathAlgebraModule, IsPathAlgebraModule]);
DeclareOperation( "DualOfModuleHomomorphism", [IsPathAlgebraModuleHomomorphism]);
DeclareOperation( "SocleOfModuleInclusion", [IsPathAlgebraModule]);
DeclareOperation( "SocleOfModule", [IsPathAlgebraModule]);
DeclareOperation( "CommonDirectSummand", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "MaximalCommonDirectSummand", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "IsomorphicModules", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "IsDirectSummand", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "IsInAdditiveClosure", [IsPathAlgebraModule, IsPathAlgebraModule ] ); 
DeclareOperation( "MorphismOnKernel", [ IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism ] );
DeclareOperation( "MorphismOnImage", [ IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism ] );
DeclareOperation( "MorphismOnCoKernel", [ IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism ] );
DeclareOperation( "LiftingMorphismFromProjective", [ IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism ] );
DeclareOperation( "LiftingInclusionMorphisms", [ IsPathAlgebraModuleHomomorphism, IsPathAlgebraModuleHomomorphism ] );
