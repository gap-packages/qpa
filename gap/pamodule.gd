# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.7 2012/04/16 09:20:23 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );
DeclareGlobalFunction("PathModuleElem");
DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraMatModule", IsAlgebraModule);
DeclareOperation("RightModuleOverPathAlgebra", [IsAlgebra, IsCollection]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);
DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);
DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);
DeclareOperation( "MatricesOfPathAlgebraModule", [IsPathAlgebraMatModule] ); 
DeclareAttribute( "DimensionVector", IsPathAlgebraMatModule ); 
DeclareAttribute( "MinimalGeneratingSetOfModule", IsPathAlgebraMatModule );
DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
DeclareOperation( "DimensionVectorPartialOrder", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );