# GAP Declarations
# This file was generated from
# $Id: algtens.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareCategory("IsAlgebraTensorElement",IsTensorElement and IsRingElement);
DeclareCategory("IsAlgebraTensor",IsAlgebra);
DeclareCategoryCollections("IsAlgebraTensorElement");
DeclareOperation("ConstituentAlgebras",[IsAlgebraTensor]);
DeclareOperation("ConstituentElements",[IsAlgebraTensorElement]);
DeclareOperation("SimplifyTensorElement",[IsAlgebraTensorElement]);
DeclareOperation("TensorElement",[IsAlgebraTensor,IsDenseList]);
DeclareOperation("TensorElement",[IsAlgebraTensor,IsRingElement,IsRingElement]);
DeclareOperation("TensorProductOfAlgebras",[IsDenseList]);
DeclareOperation("TensorProductOfAlgebras",[IsAlgebra,IsAlgebra]);
DeclareOperation("TensorProductOfAlgebraMaps",[IsDenseList]);
DeclareOperation("TensorProductOfAlgebraMaps",[IsMapping,IsMapping]);
