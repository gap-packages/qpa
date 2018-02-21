bt1 := [[["v1",1], ["v2",1], ["v3",1], ["v4",1], ["v5", 1], ["v6", 2]], [["e1", "v1", "v2"], ["e2", "v2", "v5"], ["e3", "v3", "v5"], ["e4", "v4", "v5"], ["e5", "v5", "v6"]], [["e1"], ["e2", "e1"], ["e3"], ["e4"], ["e5", "e4", "e2", "e3"], ["e5"]]];
bt2 := [[["v1",1], ["v2",1], ["v3", 3]], [["e1", "v1", "v2"], ["e2", "v2", "v3"]], [["e1"], ["e2", "e1"], ["e2"]]];

#bt3 is an example from Representation Theory by Alexander Zimmerman
bt3:=[[["v1",1], ["v2",1], ["v3", 1], ["v4", 1], ["v5", 3], ["v6", 1], ["v7", 1], ["v8", 1], ["v9", 1], ["v10", 1], ["v11",1], ["v12", 1], ["v13", 1], ["v14", 1]], [["e1", "v1", "v2"], ["e2", "v2", "v3"], ["e3", "v3", "v7"], ["e4", "v3", "v8"], ["e5", "v3", "v9"], ["e6", "v3", "v10"], ["e7", "v3", "v11"], ["e8", "v3", "v12"], ["e9", "v3", "v4"], ["e10", "v4", "v5"], ["e11", "v5", "v6"], ["e12", "v5", "v13"], ["e13", "v5", "v14"]], [["e1"], ["e2", "e1"], ["e3", "e4", "e5", "e2", "e6", "e7", "e8", "e9"], ["e10", "e9"], ["e10", "e13", "e11", "e12"], ["e11"], ["e3"], ["e4"], ["e5"], ["e6"], ["e7"], ["e8"], ["e12"], ["e13"]]];

#bgn are examples from Brauer Graph Algebras by Sibylle Schroll
bg1 := [[["a", 2], ["b", 1], ["c", 3], ["d", 1], ["e", 1], ["f", 1], ["g", 1]], [["1", "a", "c"],["2", "b", "c"], ["3", "c", "d"], ["4", "d", "e"], ["5", "d", "f"], ["6", "d", "g"]], [["1"], ["2"], ["1", "2", "3"], ["3", "4", "5", "6"], ["4"], ["5"], ["6"]]];
bg2 := [[["a", 1], ["b", 1], ["c", 1]], [["1", "a", "a"], ["2", "a", "b"], ["3", "b", "a"], ["4", "b", "c"]], [["1", "2", "3"], ["2", "4", "3"], ["4"]]];
bg3 := [[["a", 1], ["b", 1]],[["1", "a", "b"],["2", "a", "b"], ["3", "a", "b"]],[["1", "2", "3"],["1", "3", "2"]]];
bg4 := [[["a", 1], ["b", 1]],[["1", "a", "b"],["2", "a", "b"], ["3", "a", "b"]],[["1", "2", "3"],["1", "2", "3"]]];

#bcn are examples from Brauer Configuration Algebras: A Generalization of Brauer Graph Algebra by Edward L. Green and Sibylle Shroll
bc1 := [[["1", 1], ["2", 1], ["3", 1], ["4", 1], ["5", 1], ["6", 1], ["7", 1], ["8", 1]], [["V1", "1", "2", "3", "4", "7"], ["V2", "1", "2", "3", "8"], ["V3", "4", "5"], ["V4", "4", "6"], ["V5", "1", "4"]], [["V1", "V5", "V2"], ["V1", "V2"], ["V1", "V2"], ["V1", "V4", "V3", "V5"], ["V3"], ["V4"], ["V1"], ["V2"]]];
bc2 := [[["1", 1], ["2", 1], ["3", 1], ["4", 1], ["5", 2], ["6", 1]], [["V1", "1", "2", "3", "4"], ["V2", "1", "2", "3"], ["V3", "4", "5"], ["V4", "4", "6"], ["V5", "1", "4"]], [["V1", "V5", "V2"], ["V1", "V2"], ["V1", "V2"], ["V1", "V4", "V3", "V5"], ["V3"], ["V4"]]];

bt1ba := BrauerConfigurationAlgebra(Rationals, bt1);
bt2ba := BrauerConfigurationAlgebra(Rationals, bt2);

bt3ba := BrauerConfigurationAlgebra(Rationals, bt3);

bg1ba := BrauerConfigurationAlgebra(Rationals, bg1);
bg2ba := BrauerConfigurationAlgebra(Rationals, bg2);
bg3ba := BrauerConfigurationAlgebra(Rationals, bg3);
bg4ba := BrauerConfigurationAlgebra(Rationals, bg4);

bc1ba := BrauerConfigurationAlgebra(Rationals, bc1);
bc2ba := BrauerConfigurationAlgebra(Rationals, bc2);

