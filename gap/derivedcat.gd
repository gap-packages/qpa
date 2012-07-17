# properties for projective or injective complexes
DeclareProperty("IsProjectiveComplex", IsComplex);
DeclareProperty("IsInjectiveComplex", IsComplex);

# functions for computing projective resolutions of a complex
DeclareOperation("ProjectiveResolutionOfComplex",[IsComplex]);
DeclareOperation("CheckForIsomorphism",[IsPathAlgebraMatModuleHomomorphism]);
DeclareOperation("ReduceDifferential",[IsPathAlgebraMatModuleHomomorphism, IsList, IsList]);
DeclareOperation("ReducePreviousDifferential",
                 [IsPathAlgebraMatModuleHomomorphism, IsPathAlgebraMatModuleHomomorphism]);
DeclareOperation("ReduceTMap", [IsPathAlgebraMatModuleHomomorphism,
                                IsPathAlgebraMatModuleHomomorphism, IsList,
                                IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation("ReducePreviousTMap", [IsPathAlgebraMatModuleHomomorphism,
                                        IsPathAlgebraMatModuleHomomorphism]);

# functions for computing Tau of a complex
DeclareOperation("TauOfComplex",[IsComplex]);
DeclareOperation("ProjectiveToInjectiveComplex",[IsComplex]);
DeclareOperation("ProjectiveToInjectiveFiniteComplex",[IsComplex]);
DeclareOperation("StarOfMapBetweenProjectives",
                 [IsPathAlgebraMatModuleHomomorphism, IsList, IsList]);
DeclareOperation("StarOfMapBetweenIndecProjectives",
                 [IsPathAlgebraMatModuleHomomorphism, IsInt, IsList]);
DeclareOperation("StarOfMapBetweenDecompProjectives",
                 [ IsPathAlgebraMatModuleHomomorphism, IsList, IsList ]);

# functions for describing and printing complexes of projectives or injectives
DeclareOperation("DescriptionOfProjectiveModule",[IsPathAlgebraMatModule]);
DeclareOperation("CompareWithIndecProjective",[IsPathAlgebraMatModule]);
DeclareOperation("CompareWithIndecInjective", [IsPathAlgebraMatModule]);

DeclareOperation("DescriptionOfProjOrInjComplexInDegree",
                 [IsComplex, IsInt, IsBool]);
DeclareOperation("DescriptionOfProjComplexInDegree", [IsComplex, IsInt]);
DeclareOperation("DescriptionOfInjComplexInDegree", [IsComplex, IsInt]);

DeclareOperation("DescriptionOfFiniteProjOrInjComplex", [IsComplex, IsBool]);
DeclareOperation("DescriptionOfFiniteProjComplex", [IsComplex]);
DeclareOperation("DescriptionOfFiniteInjComplex", [IsComplex]);

# "help functions"
DeclareOperation("MultiplyListsOfMaps", [IsList, IsList, IsList]);
DeclareOperation("InverseOfIsomorphism", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation("MatricesOfDualMap", [IsPathAlgebraMatModuleHomomorphism]);
DeclareOperation("FindAllMapComponents", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation("DirectSumMinusSummands", [IsPathAlgebraMatModule, IsList]);