# GAP Declarations
# $Id: homomorphisms.gd,v 1.14 2012/09/28 12:57:10 sunnyquiver Exp $

DeclareCategory("IsPathAlgebraMatModuleHomomorphism", IsAdditiveElementWithZero and IsAdditiveElementWithInverse and IsGeneralMapping and RespectsAddition and RespectsZero and RespectsScalarMultiplication and IsTotal and IsSingleValued ); 
DeclareCategoryFamily(  "IsPathAlgebraMatModuleHomomorphism" );
DeclareCategoryCollections( "IsPathAlgebraMatModuleHomomorphism" );
DeclareRepresentation("IsPathAlgebraMatModuleHomomorphismRep", IsComponentObjectRep and IsAttributeStoringRep, ["maps"]);
DeclareOperation( "RightModuleHomOverAlgebra", [IsPathAlgebraMatModule, IsPathAlgebraMatModule, IsList] ); 
DeclareAttribute( "PathAlgebraOfMatModuleMap", IsPathAlgebraMatModuleHomomorphism );
DeclareOperation( "MatricesOfPathAlgebraMatModuleHomomorphism", [IsPathAlgebraMatModuleHomomorphism ] ); 
DeclareOperation( "HomOverAlgebra", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );
DeclareOperation( "HomOverAlgebraWithBasisFunction", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );
DeclareAttribute( "EndOverAlgebra", IsPathAlgebraMatModule ); 
DeclareOperation( "NumberOfNonIsoDirSummands", [IsPathAlgebraMatModule ] ); 
DeclareOperation( "LiftIdempotents", [IsPathAlgebraMatModule ] ); 
DeclareAttribute( "KernelInclusion", IsPathAlgebraMatModuleHomomorphism ); 
DeclareAttribute( "ImageProjectionInclusion", IsPathAlgebraMatModuleHomomorphism ); 
DeclareAttribute( "ImageProjection", IsPathAlgebraMatModuleHomomorphism ); 
DeclareAttribute( "ImageInclusion", IsPathAlgebraMatModuleHomomorphism ); 
DeclareAttribute( "CoKernelProjection", IsPathAlgebraMatModuleHomomorphism ); 
DeclareAttribute( "KernelOfWhat", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "CoKernelOfWhat", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "ImageOfWhat", IsPathAlgebraMatModuleHomomorphism );
DeclareProperty( "IsIsomorphism", IsPathAlgebraMatModuleHomomorphism );
DeclareOperation( "SubRepresentation", [IsPathAlgebraMatModule, IsList]);
DeclareOperation( "SubRepresentationInclusion", [IsPathAlgebraMatModule, IsList]);
DeclareOperation( "RadicalOfModule", [ IsPathAlgebraMatModule ]);
DeclareOperation( "RadicalOfModuleInclusion", [ IsPathAlgebraMatModule ]);
DeclareAttribute( "TopOfModuleProjection", IsPathAlgebraMatModule);
DeclareAttribute( "TopOfModule", IsPathAlgebraMatModule);
DeclareAttribute( "TopOfModule", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "DualOfModuleHomomorphism", IsPathAlgebraMatModuleHomomorphism);
DeclareOperation( "SocleOfModuleInclusion", [ IsPathAlgebraMatModule ]);
DeclareOperation( "SocleOfModule", [ IsPathAlgebraMatModule ]);
DeclareOperation( "CommonDirectSummand", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] ); 
DeclareOperation( "MaximalCommonDirectSummand", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] ); 
DeclareOperation( "IsomorphicModules", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] ); 
DeclareOperation( "IsDirectSummand", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] ); 
DeclareOperation( "IsInAdditiveClosure", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] ); 
DeclareOperation( "MorphismOnKernel", [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "MorphismOnImage", [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "MorphismOnCoKernel", [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "LiftingMorphismFromProjective", [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "LiftingInclusionMorphisms", [ IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "IntersectionOfSubmodules", [ IsDenseList ]);
DeclareOperation( "SumOfSubmodules", [ IsDenseList ]);
DeclareOperation( "HomFactoringThroughProjOverAlgebra", [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ]);
DeclareOperation( "EndModuloProjOverAlgebra", [ IsPathAlgebraMatModule ]);
DeclareOperation( "FromHomMMToEndM", [ IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "FromEndMToHomMM", [ IsPathAlgebraMatModule, IsMatrix ] );
DeclareProperty( "IsRightMinimal", IsPathAlgebraMatModuleHomomorphism );
DeclareProperty( "IsLeftMinimal", IsPathAlgebraMatModuleHomomorphism );
DeclareProperty( "IsSplitMonomorphism", IsPathAlgebraMatModuleHomomorphism );
DeclareProperty( "IsSplitEpimorphism", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "LeftInverseOfHomomorphism", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "RightInverseOfHomomorphism", IsPathAlgebraMatModuleHomomorphism );
DeclareOperation( "MoreRightMinimalVersion", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation( "MoreLeftMinimalVersion", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareAttribute( "RightMinimalVersion", IsPathAlgebraMatModuleHomomorphism );
DeclareAttribute( "LeftMinimalVersion", IsPathAlgebraMatModuleHomomorphism );
DeclareOperation( "HomFromProjective", [ IsRightAlgebraModuleElement, IsPathAlgebraMatModule ] );
DeclareOperation("InverseOfIsomorphism", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation( "HomomorphismFromImages", [ IsPathAlgebraMatModule, IsPathAlgebraMatModule, IsList ] );
DeclareOperation("MultiplyListsOfMaps", [IsList, IsList, IsList]);
DeclareOperation( "IsomorphismOfModules", [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ]);
DeclareOperation( "EndOfModuleAsQuiverAlgebra", [IsPathAlgebraMatModule ] );
DeclareOperation( "TraceOfModule", [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );
DeclareOperation( "RejectOfModule", [ IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );
DeclareOperation( "AllSubmodulesOfModule" , [ IsPathAlgebraMatModule ]);
DeclareOperation( "AllSimpleSubmodulesOfModule" , [ IsPathAlgebraMatModule ]);
DeclareOperation( "RestrictionViaAlgebraHomomorphismMap", [  IsAlgebraHomomorphism, IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "AllModulesOfLengthAtMost" , [ IsQuiverAlgebra, IS_INT ] );
DeclareOperation( "AllModulesOfLengthPlusOne" , [ IsList ] );
DeclareOperation( "AllIndecModulesOfLengthAtMost" , [ IsQuiverAlgebra, IS_INT ] );
DeclareAttribute( "FromIdentityToDoubleStarHomomorphism", IsPathAlgebraMatModule );
DeclareOperation( "MatrixOfHomomorphismBetweenProjectives", [ IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "FromMatrixToHomomorphismOfProjectives", [ IsQuiverAlgebra, IsMatrix, IsHomogeneousList, IsHomogeneousList ] );
DeclareOperation( "UnderlyingLinearMap",  [ IsPathAlgebraMatModuleHomomorphism ] );
DeclareOperation( "EndOfBasicModuleAsQuiverAlgebra", [ IsDenseList ] );
DeclareOperation( "ADRAlgebraOfAlgebra", [ IsQuiverAlgebra ] );
DeclareOperation( "ADRAlgebraOfBasicModule", [ IsList ] );
DeclareOperation( "MakeBasicListOfModules", [ IsList ] );
