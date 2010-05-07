# GAP Declarations
# This file was generated from
# $Id: opal.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareInfoClass("InfoOpal");
DeclareOperation( "OpalGroebnerBasis", [IsFLMLOR] );
DeclareOperation( "OpalBoundedGroebnerBasis", [IsFLMLOR,IsPosInt] );
DeclareOperation( "OpalFiniteGroebnerBasis", [IsFLMLOR] );
DeclareOperation("AlgebraToOpal", [IsFLMLOR, IsOutputStream, IsList]);
DeclareOperation("AssumePathAlgebra", [IsField, IsOutputStream]);
DeclareOperation("WriteCoefficientToOpal", [IsScalar, IsOutputStream]);
DeclareOperation("WriteRelatorsToOpal", [IsFLMLOR, IsOutputStream, IsList]);
DeclareGlobalFunction( "WriteFinalOpalCommands" );
DeclareGlobalFunction( "WriteFinalOpalCommandsWithBound" );
DeclareGlobalFunction( "WriteFinalOpalCommandsFinite" );
DeclareOperation("OrderingToOpal", [IsQuiverOrdering]);
DeclareGlobalFunction( "ExecuteOpal" );
