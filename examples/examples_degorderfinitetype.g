a := ARQuiverNumerical(12, 4, [ [0],[1,0],[1,0],[1,0],[2,3,4,1],[5,2],[5,3],[5,4],[6,7,8,5],[9,6],[9,7],[9,8] ]);
b := ARQuiverNumerical(12, 4, ["orbits", [ [2], [2,[1,0]], [2,[1,0]], [2,[1,0]] ] ]);
###############
a := ARQuiverNumerical("D4 subspace");
###############
a := ARQuiverNumerical("BG", 5);
###############
a := ARQuiverNumerical("D4 subspace");
DimensionVector(a, 7);
DimensionVector(a, [0,1,0,0,0,0,2,0,0,0,0,0]);
###############
a := ARQuiverNumerical("R nilp");
DimensionVector(a, 2);  DimensionVector(a, 3);
DegOrderLEQ(a, 2, 3);
DegOrderLEQ(a, 3, 2);
###############
a := ARQuiverNumerical("BG", 5);
preds := DegOrderPredecessors(a, 60);; Length(preds);
DegOrderLEQ(a, preds[7], 60);
dpreds := DegOrderDirectPredecessors(a, 60);; Length(dpreds);
PrintMultiplicityVectors(dpreds);