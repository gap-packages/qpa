InstallMethod( RightModuleCategory,
               [ IsQuiverAlgebra ],
function( A )
  local category;

  category := CreateCapCategory( Concatenation( "Category of right modules over ",
                                                String( A ) ) );
  category!.path_algebra := A;

  SetIsAbelianCategory( category, true );

  AddPreCompose( category, \* );
  AddIsEqualForObjects( category, \= );
  AddIdentityMorphism( category, IdentityMapping );
  AddMonoAsKernelLift( category,
  function( monomorphism, test_morphism )
    local im_proj, im_inc;
    im_proj := ImageProjection( test_morphism );
    im_inc := ImageInclusion( test_morphism );
    return im_proj * LiftingInclusionMorphisms( monomorphism, im_inc );
  end );
  AddKernelEmb( category, KernelInclusion );
#  AddCokernel( category, CoKernel );
  AddCokernelProj( category, CoKernelProjection );
  AddCokernelColiftWithGivenCokernel( category,
  function( morphism, test_morphism, coker )
    local coker_proj, test_object, test_matrices,  
          proj_matrices, proj_inv_matrices, lift_matrices;
    coker_proj := CokernelProj( morphism );
    test_object := Range( test_morphism );
    test_matrices := MatricesOfPathAlgebraMatModuleHomomorphism( test_morphism );
    proj_matrices := MatricesOfPathAlgebraMatModuleHomomorphism( coker_proj );
    proj_inv_matrices := List( proj_matrices,
                               M -> List( IdentityMat( DimensionsMat( M )[ 2 ] ),
                                          e -> SolutionMat( M, e ) ) );
    lift_matrices := ListN( proj_inv_matrices, test_matrices, \* );
    return RightModuleHomOverAlgebra( coker, test_object, lift_matrices );
  end );
  # AddCokernelColift( category,
  # function( morphism, test_morphism )
  #   local ker_emb;
  #   ker_emb := KernelEmb( test_morphism );
  #    MorphismOnCoKernel( morphism, ker_emb,
  #                              KernelLift( morphism, ker_emb ),
  #                              IdentityMorphism( Range( morphism ) ) )
  AddAdditionForMorphisms( category,
  function( f, g )
    local i, num_vert, x, Fam;
    num_vert := Length(f!.maps);
    x := List([1..num_vert], y -> f!.maps[y] + g!.maps[y]);
    return RightModuleHomOverAlgebra(Source(f),Range(f),x);
  end );
  AddAdditiveInverseForMorphisms( category,
  function( f )
    local i, num_vert, x;
    num_vert := Length(f!.maps);
    x := List([1..num_vert], y -> (-1)*f!.maps[y]);
    return RightModuleHomOverAlgebra(Source(f),Range(f),x);
  end );
  AddZeroMorphism( category, ZeroMapping );
  AddIsZeroForMorphisms( category, IsZero );
  AddIsEqualForMorphisms( category, \= );
  AddZeroObject( category, function() return ZeroModule( A ); end );
  AddDirectSum( category, DirectSumOfQPAModules );
  AddInjectionOfCofactorOfCoproductWithGivenCoproduct( category,
  function( summands, i, sum )
    return DirectSumInclusions( sum )[ i ];
  end );
  AddProjectionInFactorOfDirectProductWithGivenDirectProduct( category,
  function( summands, i, sum )
    return DirectSumProjections( sum )[ i ];
  end );
  
  SetIsAbelianCategory( category, true );
  Finalize( category );

  return category;
end );
