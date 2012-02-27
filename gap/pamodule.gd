# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.6 2012/02/27 12:26:34 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );
DeclareGlobalFunction("PathModuleElem");
DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraModule", IsAlgebraModule);
DeclareOperation("RightModuleOverPathAlgebra", [IsAlgebra, IsCollection]);
#DeclareOperation( "NewRightModuleOverPathAlgebra", [IsAlgebra, IsList, IsList ]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);
DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);
DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);
DeclareOperation( "MatricesOfPathAlgebraModule", [IsPathAlgebraModule] ); 
DeclareAttribute( "DimensionVector", IsPathAlgebraModule ); 
DeclareAttribute( "MinimalGeneratingSetOfModule", IsPathAlgebraModule );
DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
DeclareOperation( "DimensionVectorPartialOrder", [IsPathAlgebraModule, IsPathAlgebraModule ] );