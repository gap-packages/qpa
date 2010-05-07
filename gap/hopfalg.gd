# GAP Declarations
# This file was generated from
# $Id: hopfalg.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareCategory("IsHopfAlgebra", IsAlgebra);
DeclareCategory("IsHopfAlgebraElement", IsRingElement);
DeclareCategoryCollections("IsHopfAlgebraElement");
DeclareCategoryFamily("IsHopfAlgebraElement");
DeclareRepresentation("IsHopfAlgebraBasisDefaultRep",IsComponentObjectRep,["underlyingBasis"]);
DeclareCategory("IsHopfAlgebraBasis",IsBasis and IsHopfAlgebraBasisDefaultRep);
DeclareAttribute("UnderlyingAlgebra", IsHopfAlgebra);
DeclareAttribute("ComultiplicationMap", IsHopfAlgebra);
DeclareAttribute("CounitMap", IsHopfAlgebra);
DeclareAttribute("AntipodeMap", IsHopfAlgebra);
DeclareOperation("Comultiply", [IsRingElement]);
DeclareOperation("Counit", [IsRingElement]);
DeclareOperation("Antipode", [IsRingElement]);
DeclareOperation("HopfAlgebra", [IsGroupRing]);
DeclareOperation("HopfAlgebraNC", [IsGroupRing]);
DeclareOperation("HopfAlgebraNC", [IsAlgebra,IsAlgebraGeneralMapping,IsAlgebraGeneralMapping,IsAlgebraGeneralMapping]);
DeclareCategory("IsAlgebraAntihomomorphism",IsAlgebraGeneralMapping);
DeclareCategory("IsAlgebraWithOneAntihomomorphism",IsAlgebraWithOneGeneralMapping and IsAlgebraAntihomomorphism);
DeclareOperation("AlgebraAntihomomorphismByImages", [IsAlgebra,IsAlgebra,IsDenseList,IsDenseList]);
DeclareOperation("AlgebraAntihomomorphismByImagesNC", [IsAlgebra,IsAlgebra,IsDenseList,IsDenseList]);
DeclareOperation("AlgebraWithOneAntihomomorphismByImages", [IsAlgebraWithOne,IsAlgebraWithOne,IsDenseList,IsDenseList]);
DeclareOperation("AlgebraWithOneAntihomomorphismByImagesNC", [IsAlgebraWithOne,IsAlgebraWithOne,IsDenseList,IsDenseList]);
