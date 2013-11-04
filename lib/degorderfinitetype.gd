# GAP Implementation
# This contains various tools for computing degenerations of modules in finite type
# (created A. Mroz, 01.07.2013)

DeclareCategory("IsARQuiverNumerical", IsObject ); 
DeclareAttribute("NumberOfIndecomposables", IsARQuiverNumerical );
DeclareAttribute("NumberOfProjectives", IsARQuiverNumerical );

DeclareGlobalFunction("PredefARQuivers");
DeclareGlobalFunction("ARQuiverNumerical");


DeclareOperation("ModulesOfDimVect", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DegOrderPredecessors", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DegOrderDirectPredecessors", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DegOrderPredecessorsWithDirect", [IsARQuiverNumerical, IsObject]);

DeclareOperation("DegOrderSuccessors", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DegOrderDirectSuccessors", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DegOrderSuccessorsWithDirect", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DimensionVector", [IsARQuiverNumerical, IsObject]);
DeclareOperation("DimHom", [IsARQuiverNumerical, IsObject, IsObject]);
DeclareOperation("DimEnd", [IsARQuiverNumerical, IsObject]);
DeclareOperation("OrbitDim", [IsARQuiverNumerical, IsObject]);
DeclareOperation("OrbitCodim", [IsARQuiverNumerical, IsObject, IsObject]);
DeclareOperation("DegOrderLEQ", [IsARQuiverNumerical, IsObject, IsObject]);
DeclareOperation("DegOrderLEQNC", [IsARQuiverNumerical, IsObject, IsObject]);

DeclareOperation("PrintMultiplicityVector", [IsList]);
DeclareOperation("PrintMultiplicityVectors", [IsList]);


DeclareOperation("TestCodimConjecture", [IsARQuiverNumerical, IsBool, IsBool]);