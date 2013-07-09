# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.8 2012/08/01 16:01:10 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );
DeclareGlobalFunction("PathModuleElem");
DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraMatModule", IsAlgebraModule);
DeclareOperation("RightModuleOverPathAlgebra", [IsQuiverAlgebra, IsCollection]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);
DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);
DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);
DeclareOperation( "MatricesOfPathAlgebraModule", [IsPathAlgebraMatModule] ); 
DeclareAttribute( "DimensionVector", IsPathAlgebraMatModule ); 
DeclareAttribute( "MinimalGeneratingSetOfModule", IsPathAlgebraMatModule );
DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
DeclareOperation( "DimensionVectorPartialOrder", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );
DeclareAttribute( "AnnihilatorOfModule", IsPathAlgebraMatModule );
DeclareAttribute( "LoewyLength", IsPathAlgebraMatModule  ); 
DeclareOperation( "RadicalSeries",  [ IsPathAlgebraMatModule ] ); 
DeclareOperation( "SocleSeries", [ IsPathAlgebraMatModule ] ); 
#DeclareOperation( "Dimension", [ IsPathAlgebraMatModule ] ); 
DeclareProperty( "IsIndecomposableModule",  IsPathAlgebraMatModule ); 
DeclareProperty( "IsProjectiveModule",  IsPathAlgebraMatModule ); 
DeclareProperty( "IsInjectiveModule", IsPathAlgebraMatModule ); 
DeclareProperty( "IsSimpleModule", IsPathAlgebraMatModule ); 
InstallTrueMethod( IsIndecomposableModule, IsSimpleModule );
DeclareProperty( "IsSemisimpleModule", IsPathAlgebraMatModule );
InstallTrueMethod( IsSemisimpleModule, IsSimpleModule );
DeclareOperation( "DirectSumOfModules", [ IsList] );
DeclareAttribute( "IsDirectSumOfModules", IsPathAlgebraMatModule );
DeclareAttribute( "DirectSumProjections", IsPathAlgebraMatModule );
DeclareAttribute( "DirectSumInclusions", IsPathAlgebraMatModule );
DeclareOperation( "SupportModuleElement", [ IsRightAlgebraModuleElement ] );
DeclareOperation( "RightAlgebraModuleToPathAlgebraMatModule", [ IsRightAlgebraModuleElementCollection ]);
DeclareProperty( "IsRigidModule", IsPathAlgebraMatModule );
DeclareProperty( "IsTauRigidModule", IsPathAlgebraMatModule );
DeclareProperty( "IsExceptionalModule", IsPathAlgebraMatModule );
InstallTrueMethod( IsExceptionalModule, IsIndecomposableModule and IsRigidModule );
