# GAP Declarations
# This file was generated from
# $Id: decomp.gd,v 1.3 2012/08/15 07:33:53 sunnyquiver Exp $
DeclareOperation("ComplementInFullRowSpace", [IsRowSpace]);
DeclareOperation("LiftIdempotentsForDecomposition",[IsAlgebraGeneralMapping,IsList]);
DeclareOperation("IdempotentsForDecomposition",[IsAlgebra]);
DeclareAttribute("DecomposeModule", IsPathAlgebraMatModule );
DeclareAttribute("DecomposeModuleWithInclusions", IsPathAlgebraMatModule );
DeclareAttribute("DecomposeModuleWithMultiplicities", IsPathAlgebraMatModule);
DeclareOperation( "DecomposeModuleProbabilistic", [ IsHomogeneousList, IsPathAlgebraMatModule ] );
DeclareOperation( "LiftIdempotent", [ IsAlgebraGeneralMapping, IsRingElement ] );
DeclareOperation( "LiftTwoOrthogonalIdempotents", [ IsAlgebraGeneralMapping, IsRingElement, IsRingElement ] );
DeclareOperation( "BlockSplittingIdempotents", [IsPathAlgebraMatModule ] ); 
DeclareOperation( "BlockDecompositionOfModule", [IsPathAlgebraMatModule ] ); 
DeclareOperation( "BasicVersionOfModule", [ IsPathAlgebraMatModule ] ); 
DeclareOperation( "DecomposeModuleViaTop", [ IsPathAlgebraMatModule ] );
DeclareOperation( "DecomposeModuleViaCharPoly", [ IsPathAlgebraMatModule ] );