# GAP Declarations
# This file was generated from
# $Id: pamodule.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareCategory("IsPathModuleElem", IsVector);
DeclareCategoryFamily( "IsPathModuleElem" );
DeclareCategoryCollections( "IsPathModuleElem" );

DeclareGlobalFunction("PathModuleElem");

DeclareCategory("IsBasisOfPathModuleElemVectorSpace", IsBasis);
DeclareProperty("IsPathAlgebraMatModule", IsAlgebraModule);

DeclareOperation("RightModuleOverPathAlgebra", [IsRing, IsCollection]);
DeclareOperation("SubmoduleAsModule", [IsAlgebraModule]);

DeclareOperation("BasisOfDomain", [IsFreeLeftModule and IsPathModuleElemCollection]);

DeclareOperation("NewBasis", [ IsFreeLeftModule and IsPathModuleElemCollection,
      IsPathModuleElemCollection and IsList ]);

DeclareFilter( "IsAlgebraModuleHomomorphism", IsLeftModuleGeneralMapping );
