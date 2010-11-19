# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.3 2010/11/19 13:24:48 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );

DeclareGlobalFunction("PathModuleElem");

DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraMatModule", IsAlgebraModule);

DeclareOperation("RightModuleOverPathAlgebra", [IsRing, IsCollection]);
DeclareOperation("RightModuleOverQuotientOfPathAlgebra", [IsSubalgebraFpPathAlgebra, IsCollection]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);

DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);

DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);
DeclareOperation( "MatricesOfPathAlgebraMatModule", [IsPathAlgebraMatModule] ); 
DeclareOperation( "DimensionVector", [IsPathAlgebraMatModule] ); 

DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
