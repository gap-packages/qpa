# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.5 2012/02/02 10:07:31 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );

DeclareGlobalFunction("PathModuleElem");

DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraMatModule", IsAlgebraModule);

DeclareOperation("RightModuleOverPathAlgebra", [IsRing, IsCollection]);
DeclareOperation("RightModuleOverQuotientOfPathAlgebra", [IsSubalgebraFpPathAlgebra, IsCollection]);
DeclareOperation( "NewRightModuleOverPathAlgebra", [IsAlgebra, IsList, IsList ]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);

DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);

DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);
DeclareOperation( "MatricesOfPathAlgebraMatModule", [IsPathAlgebraMatModule] ); 
DeclareOperation( "DimensionVector", [ IsPathAlgebraMatModule ] ); 
DeclareAttribute( "MinimalSetOfGenerators", IsPathAlgebraMatModule );
DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
DeclareOperation( "DimensionVectorPartialOrder", [IsPathAlgebraMatModule, IsPathAlgebraMatModule ] );