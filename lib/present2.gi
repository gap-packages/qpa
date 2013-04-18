# GAP Implementation
# This file was generated from 
# $Id: present2.gi,v 1.1 2010/05/07 13:30:13 sunnyquiver Exp $
#


InstallMethod( ProjectivePresentation,
  "for algebras",
  true,
  [ IsAlgebra, IsMatrix ], 0,
  function(A,mat)
    local P0,P1,P0List,P1List,M,obj,m,n,i,v,terms,nzterm,nzentry;

    # Number of rows in matrix:
    m := Length(mat);
    if (m < 1) then
      Error("Representation matrix cannot be empty.\n");
    fi;

    # Check each row length to insure column consistency in matrix:
    n := Length(mat[1]);
    i := 2;
    while (i <= m) do
      if (Length(mat[i]) <> n) then
        Error("Matrix must have rows of equal length.\n");
      fi;
      i := i + 1;
    od;

    # Initialize list for projective module P0 in presentation:
    # NOTE: this is the starting vertex in each row, and should be
    #       the same for each element in each row.
    P0List := [];    

    # Check that each row of matrix has uniformly same starting vertex:
    i := 1;
    while (i <= m) do

      # Find first nonzero entry in row and get Source vertex:
      nzentry := First(mat[i], x->not IsZero(x));
      nzterm := CoefficientsAndMagmaElements(nzentry);
      v := SourceOfPath(nzterm[1]);

      if IsLeftUniform(mat[i],v) then
        P0List[i] := v;
      else
        Error("All rows of matrix must contain paths left uniform,\n",
              "starting with the same vertex.\n");
      fi;

      i := i + 1;
    od;
    
    # Initialize list of vertices for projective module P1 in presentation:
    P1List := [];    
    

    # Check that each column of matrix has uniformly same starting vertex:

    obj := rec(); 
    ObjectifyWithAttributes( obj,
        NewType(ProjectivePresentationFamily,
                        IsProjectivePresentationDefaultRep),
        ParentAlgebra,A,
        MatrixRepresentation,mat,
	MatrixNumRows,m,
	MatrixNumCols,n,
        P0VertexList,P0List
      );

    return obj;

  end
);
