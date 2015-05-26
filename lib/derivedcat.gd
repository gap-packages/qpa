# properties for projective or injective complexes
DeclareProperty("IsProjectiveComplex", IsQPAComplex);
DeclareProperty("IsInjectiveComplex", IsQPAComplex);

# functions for computing projective resolutions of a complex
DeclareOperation("ProjectiveResolutionOfComplex",[IsQPAComplex]);
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
DeclareOperation("TauOfComplex",[IsQPAComplex]);
DeclareOperation("ProjectiveToInjectiveComplex",[IsQPAComplex]);
DeclareOperation("ProjectiveToInjectiveFiniteComplex",[IsQPAComplex]);
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
                 [IsQPAComplex, IsInt, IsBool]);
DeclareOperation("DescriptionOfProjComplexInDegree", [IsQPAComplex, IsInt]);
DeclareOperation("DescriptionOfInjComplexInDegree", [IsQPAComplex, IsInt]);

DeclareOperation("DescriptionOfFiniteProjOrInjComplex", [IsQPAComplex, IsBool]);
DeclareOperation("DescriptionOfFiniteProjComplex", [IsQPAComplex]);
DeclareOperation("DescriptionOfFiniteInjComplex", [IsQPAComplex]);

# "help functions"
#DeclareOperation("MultiplyListsOfMaps", [IsList, IsList, IsList]);
DeclareOperation("MatricesOfDualMap", [IsPathAlgebraMatModuleHomomorphism]);
DeclareOperation("FindAllMapComponents", [ IsPathAlgebraMatModuleHomomorphism ]);
DeclareOperation("DirectSumMinusSummands", [IsPathAlgebraMatModule, IsList]);