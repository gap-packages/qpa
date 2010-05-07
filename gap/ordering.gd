# GAP Declarations
# This file was generated from
# $Id: ordering.gd,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
DeclareCategory("IsQuiverOrdering", IsObject);
DeclareCategory("IsLexicographicOrdering", IsQuiverOrdering);
DeclareCategory("IsLengthOrdering", IsQuiverOrdering);
DeclareCategory("IsReverseOrdering", IsQuiverOrdering);
DeclareCategory("IsVectorOrdering", IsQuiverOrdering);
DeclareCategory("IsWeightOrdering",IsQuiverOrdering);
DeclareCategory("IsBlockOrdering",IsQuiverOrdering);
DeclareCategory("IsWreathOrdering",IsQuiverOrdering);
DeclareCategory("IsLeftVectorOrdering", IsVectorOrdering);
DeclareCategory("IsRightVectorOrdering", IsVectorOrdering);
DeclareCategory("IsLeftLexicographicOrdering", IsLexicographicOrdering);
DeclareCategory("IsRightLexicographicOrdering", IsLexicographicOrdering);

DeclareProperty("IsWellOrdering",IsQuiverOrdering);
DeclareProperty("IsWellReversedOrdering",IsQuiverOrdering);
DeclareProperty("IsTotalOrdering",IsQuiverOrdering);
DeclareProperty("IsAdmissibleOrdering",IsQuiverOrdering);

DeclareAttribute("ComparisonFunction", IsQuiverOrdering);
DeclareAttribute("NextOrdering",IsQuiverOrdering);
DeclareAttribute("OrderingName", IsQuiverOrdering );
DeclareAttribute("LexicographicTable", IsQuiverOrdering);

DeclareOperation("LessThanByOrdering", [IsQuiverOrdering, IsObject, IsObject]);
DeclareOperation("IsWellSubset",[IsQuiverOrdering,IsHomogeneousList]);
DeclareOperation("IsWellReversedSubset",[IsQuiverOrdering,IsHomogeneousList]);
DeclareOperation("IsTotalSubset",[IsQuiverOrdering,IsHomogeneousList]);

DeclareGlobalFunction("LengthOrdering");
DeclareGlobalFunction("LeftLexicographicOrdering");
DeclareGlobalFunction("RightLexicographicOrdering" );
DeclareGlobalFunction("ReverseOrdering");
DeclareGlobalFunction("LeftVectorOrdering");
DeclareGlobalFunction("RightVectorOrdering");
DeclareGlobalFunction("WeightOrdering");
DeclareGlobalFunction("BlockOrdering");
DeclareGlobalFunction("WreathOrdering");
