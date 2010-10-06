# GAP Implementation
# This file was generated from 
# $Id: predefalgs.gi,v 1.2 2010/10/06 05:29:16 sunnyquiver Exp $
InstallMethod ( NakayamaAlgebra,
    "for an admissible sequence and a field",
    [IsList, IsField], 0,
    function( admiss_seq , K)
        local num_vert, # number of vertices in the quiver we are creating.
              test,     # the result of the test of being an admissible seq.
              i, j,     # integer variables for for-loops.
              q, r,     # for writing a number n: n = q*num_vert + r.
              Q, KQ,    # the quiver and the path algebra for the 
                        # Nakayama algebra we are creating.
              rels,     # relations of the Nakayama algebra we are creating. 
              new_rel,  # partial relation. 
              arrows,   # arrows of the quiver we are creating.
              cycle,    # if we are in the A-tilde case and the relations.
                        # are long, this is used to store the cycles involved
                        # in the relations. 
              I;        # the ideal of relations.
#
# Testing if we have an admissible sequence
#
        num_vert := Length(admiss_seq);
        test := true;
        i := 1;
        while ( test and ( i < num_vert ) ) do
            test := ( ( admiss_seq[i+1] >= admiss_seq[i] - 1 ) and
                ( admiss_seq[i] - 1 >= 1 ) );  
            i := i + 1;
        od;
        if test then 
            test := ( admiss_seq[1] >= admiss_seq[num_vert] - 1 );
        fi;
#
# If an admissible sequence, the case admiss_seq[num_vert] = 1, is the
# A_n case, and the other case is A_n-tilde case. 
#
        if test then 
            if ( admiss_seq[num_vert] = 1 ) then
#
# This is the A_n case.
#
                arrows := [];
                for i in [1..num_vert - 1] do
                    Add(arrows,[i,i+1]);
                od;
                Q := Quiver(num_vert,arrows);
                KQ := PathAlgebra(K,Q); 
                arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert-1]}; 
                rels := [];
                for i in [1..num_vert - 1] do
                    if (( admiss_seq[i] < num_vert - i + 1 ) and 
                        ( admiss_seq[i] - 1 <> admiss_seq[i+1] )) then
                        new_rel := One(KQ);
                        for j in [i..i+admiss_seq[i] - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                        Add(rels,new_rel); 
                    fi;
                od;
#
# end of A_n case
#
            else
#
# This is the A_n-tilde case.
                arrows := [];
                for i in [1..num_vert - 1] do
                    Add(arrows,[i,i+1]);
                od;
                Add(arrows,[num_vert,1]); 
                Q := Quiver(num_vert,arrows);
                KQ := PathAlgebra(K,Q); 
                arrows := GeneratorsOfAlgebra(KQ){[num_vert+1..2*num_vert]}; 
                rels := [];
#
# Relations starting in the vertices 1..num_vert - 1. 
#
                for i in [1..num_vert - 1] do
                    if ( admiss_seq[i] - 1 <> admiss_seq[i+1] ) then
                        new_rel := One(KQ);
                        r := admiss_seq[i] mod num_vert;
                        q := (admiss_seq[i] - r)/num_vert;
                        cycle := One(KQ);
                        for j in [i..num_vert] do
                            cycle := cycle*arrows[j];
                        od;
                        for j in [1..i-1] do
                            cycle := cycle*arrows[j];
                        od;
                        for j in [1..q] do 
                            new_rel := new_rel*cycle;
                        od;
                        if ( i + r - 1 <= num_vert ) then
                            for j in [i..i+r - 1] do
                                new_rel := new_rel*arrows[j];
                            od;
                        else
                            for j in [i..num_vert] do
                                new_rel := new_rel*arrows[j];
                            od;
                            for j in [1..r - (num_vert - i) - 1] do
                                new_rel := new_rel*arrows[j];
                            od;
                        fi;
                        Add(rels,new_rel); 
                    fi;
                od;
#
# Relation starting in the vertex num_vert. 
#
                if ( admiss_seq[num_vert] - 1 <> admiss_seq[1] ) then
                    new_rel := One(KQ);
                    r := admiss_seq[num_vert] mod num_vert;
                    q := (admiss_seq[num_vert] - r)/num_vert;
                    cycle := One(KQ)*arrows[num_vert];
                    for j in [1..num_vert-1] do
                        cycle := cycle*arrows[j];
                    od;
                    for j in [1..q] do 
                        new_rel := new_rel*cycle;
                    od;
                    if ( r <> 0 ) then
                        new_rel := new_rel*arrows[num_vert];
                        for j in [1..r - 1] do
                            new_rel := new_rel*arrows[j];
                        od;
                    fi;
                    Add(rels,new_rel);
                fi;
#
# End of A_n-tilde case.
#
            fi;
            I := Ideal(KQ,rels);
            return KQ/I; 
        else
#
# A valid admissible sequence was not entered. 
# 
            Print("The admissible sequence is not a valid one.\n");
            return admiss_seq;
        fi; 

    end
);
