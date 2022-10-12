InstallGlobalFunction( QPAStringQuiverRho,
  function(Q,rho)
    local rho_q, i,j,k,arr,arr1;
    rho_q := [];
    for i in [1..Length(rho)] do
        arr := [];
        for j in [1..Length(rho[i])] do
            for k in [1..NumberOfArrows(Q)] do
                if rho[i][j] = String(ArrowsOfQuiver(Q)[k])[1] then
                    Add(arr,ArrowsOfQuiver(Q)[k]);
                fi;
            od;
        od;
        arr1 := arr[1];
        for k in [2..Length(arr)] do
            arr1 := arr1*arr[k];
        od;
        Add(rho_q, arr1);
    od;
    return rho_q;
end
);

InstallGlobalFunction( QPAStringSigmaEps,
  function( A )
    local Q, rho1, rho, i, onerel, tworel, str, Composable_Matrix, IsComposable, temp_in, temp_out,
    j, k, sigma, eps, x, y;

    Q := QuiverOfPathAlgebra(A);
    rho := RelationsOfAlgebra(A);

    rho1 := [];
    for i in [1..Length(rho)] do
        onerel := CoefficientsAndMagmaElements(rho[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho1, [str]);
    od;

    Composable_Matrix := function(Q,rho)
        local mat, x, str, i, j;
        mat := [];
        for i in [1..NumberOfArrows(Q)] do
            str := [];
            for j in [1..NumberOfArrows(Q)] do
                if SourceOfPath(ArrowsOfQuiver(Q)[i]) =
                  TargetOfPath(ArrowsOfQuiver(Q)[j]) and not
                  \in(ArrowsOfQuiver(Q)[j]*ArrowsOfQuiver(Q)[i],rho) then
                    Append(str, [1]);
                else
                    Append(str, [0]);
                fi;
            od;
            Append(mat, [str]);
        od;
        return mat;
    end;

    IsComposable := function(Q,rho,arrow1,arrow2)
        local i, Mat, pos1, pos2, j;
        Mat := Composable_Matrix(Q,rho);
        for i in [1..NumberOfArrows(Q)] do
            if arrow1 = ArrowsOfQuiver(Q)[i] then
                break;
            fi;
        od;
        for j in [1..NumberOfArrows(Q)] do
            if arrow2 = ArrowsOfQuiver(Q)[j] then
                break;
            fi;
        od;
        return Mat[i][j];
    end;

    temp_in := function(Q,i,j)
        return SIntChar(String(IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1])
        - SIntChar('a') + 1;
    end;

    temp_out := function(Q,i,j)
        return SIntChar(String(OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1])
        - SIntChar('a') + 1;
    end;

    sigma := [];
    eps := [];
    for i in [1..NumberOfArrows(Q)] do
        Append(sigma, [0]);
        Append(eps, [0]);
    od;

    rho := QPAStringQuiverRho(Q,rho1);

    for i in [1..NumberOfVertices(Q)] do
        x := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
        y := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
        if Length(x) = 2 then
            eps[temp_in(Q, i, 1)] := 1;
            eps[temp_in(Q, i, 2)] := -1;
            if Length(y) = 1 then
                if IsComposable(Q, rho, y[1], x[1]) = 1 then
                    sigma[temp_out(Q, i, 1)] := -1;
                elif IsComposable(Q, rho, y[1], x[2]) = 1 then
                    sigma[temp_out(Q, i, 1)] := 1;
                fi;
            elif Length(y) = 2 then
                if IsComposable(Q, rho, y[1], x[1]) = 1 or
                  IsComposable(Q, rho, y[2], x[2]) = 1 then
                    sigma[temp_out(Q, i, 1)] := -1;
                    sigma[temp_out(Q, i, 2)] := 1;
                else
                    sigma[temp_out(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 2)] := -1;
                fi;
            fi;
        elif Length(x) = 1 then
            if Length(y) = 1 then
                eps[temp_in(Q, i, 1)] := 1;
                sigma[temp_out(Q, i, 1)] := -1;
            elif Length(y) = 2 then
                sigma[temp_out(Q, i, 1)] := 1;
                sigma[temp_out(Q, i, 2)] := -1;
                if IsComposable(Q, rho, y[1], x[1]) = 1 then
                    eps[temp_in(Q, i, 1)] := -1;
                else
                    eps[temp_in(Q, i, 1)] := 1;
                fi;
            fi;
        elif Length(x) = 0 then
          sigma[temp_out(Q, i, 1)] := 1;
          sigma[temp_out(Q, i, 2)] := -1;
        fi;
    od;
    return [sigma, eps];
end
);

InstallMethod( IsValidString,
"for stringalgebras and strings",
true,
[ IsAlgebra, IsString ], 0,
function( A, input_str )
    local IsValidStringSA, IsReducedWalk, IsValidWalk, Q, rho1, rho, vertices, arrows, i, x, y, str;
    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho := [];
    for i in [1..Length(rho1)] do
        x := CoefficientsAndMagmaElements(rho1[i]);
        y := WalkOfPath(x[1]);
        str := Concatenation(List(y, a -> String(a)));
        Append(rho, [str]);
    od;

    IsReducedWalk := function(Q, input_str)
        local x,y,i,j;
        x := 0;
        y := ['-','1'];
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            Append(y, String(ArrowsOfQuiver(Q)[i]));
        od;
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            Append(y, [CharInt(SIntChar(String(ArrowsOfQuiver(Q)[i])[1])-32)]);
        od;

        for i in [1..Length(input_str)] do
            if not input_str[i] in y then
                return 0;
            fi;
        od;
        for i in [1..Length(input_str) - 1] do
            if SIntChar(input_str[i]) - SIntChar(input_str[i+1]) = 32 or
              SIntChar(input_str[i + 1]) - SIntChar(input_str[i]) = 32 then
                x := 1;
            else continue;
            fi;
        od;
        if x = 0 then return 1;
        else return 0;
        fi;
    end;

    IsValidWalk := function (Q, input_str)
        local i,j,k,x,s,t,array, y;
        x := 0;
        y := [];
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            Append(y, String(ArrowsOfQuiver(Q)[i]));
        od;
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            Append(y, [CharInt(SIntChar(String(ArrowsOfQuiver(Q)[i])[1]) - 32)]);
        od;
        for i in [1..Length(input_str)] do
            if not input_str[i] in y then
                return false;
            fi;
        od;
        for i in [1..Length(input_str) - 1] do
            if SIntChar(input_str[i]) > 96 then
                if SIntChar(input_str[i + 1]) > 96 then
                    if SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i])
                      - SIntChar('a') + 1]) <>
                      TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i + 1])
                      - SIntChar('a') + 1]) then
                        x := 1;
                    fi;
                else
                    if SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i])
                      - SIntChar('a') + 1]) <>
                      SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i + 1])
                      - SIntChar('A') + 1]) then
                        x := 1;
                    fi;
                fi;
            else
                if SIntChar(input_str[i + 1]) > 96 then
                    if TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i])
                      - SIntChar('A') + 1]) <>
                      TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i + 1])
                      - SIntChar('a') + 1]) then
                        x := 1;
                    fi;
                else
                    if TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i])
                      - SIntChar('A') + 1]) <>
                      SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[i + 1])
                      - SIntChar('A') + 1]) then
                        x := 1;
                    fi;
                fi;
            fi;
        od;
        if x = 0 then
            if SIntChar(input_str[Length(input_str)]) > 96 then
            s := SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[Length(input_str)])
            - SIntChar('a') + 1]);
            else
            s := TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[Length(input_str)])
            - SIntChar('A') + 1]);
            fi;

            if SIntChar(input_str[1]) > 96 then
            t := TargetOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[1])
            - SIntChar('a') + 1]);
            else
            t := SourceOfPath(ArrowsOfQuiver(Q)[SIntChar(input_str[1])
            - SIntChar('A') + 1]);
            fi;
            array := [s, t];
        fi;

        if x = 1 then return false;
        else return true;
        fi;
    end;

    IsValidStringSA := function(A, Q, rho1, input_str)
        local i, j, k, rel, temp, temp_caps, temp1, y, x, l, length, count, num,
        num_digits, lit, sigma, eps, array, rho;
        rho := QPAStringQuiverRho(Q,rho1);
        array := QPAStringSigmaEps(A);
        sigma := array[1];
        eps := array[2];
        num := NumberOfVertices(Q);
        num_digits := 0;
        while num <> 0 do
            num := num - num mod 10;
            num := num/10;
            num_digits := num_digits + 1;
        od;
        y := "0123456789(-,)";
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            y := Concatenation(String(ArrowsOfQuiver(Q)[i]), y);
        od;
        for i in [1..Length(ArrowsOfQuiver(Q))] do
            y := Concatenation([CharInt(SIntChar(String(
            ArrowsOfQuiver(Q)[i])[1]) - 32)], y);
        od;
        for i in [1..Length(input_str)] do
            if not input_str[i] in y then
                return false;
            fi;
        od;
        if input_str[1] = '(' then
            num := 0;
            for j in [2..num_digits + 1] do
                num := num + (SIntChar(input_str[j])
                - SIntChar('0')) * 10^(num_digits - j + 1);
        		od;
            if \in(num,[1..NumberOfVertices(Q)]) and input_str[num_digits + 2] =
              ',' and input_str[num_digits + 3] = '1' and
              input_str[num_digits+4] = ')' then
                return true;
            elif \in(num,[1..NumberOfVertices(Q)]) and input_str[num_digits + 2]
              = ',' and input_str[num_digits + 3] = '-' and
              input_str[num_digits + 4] = '1' and
              input_str[num_digits + 5] = ')' then
                return true;
            fi;
        fi;
        for i in [1..14] do
            Unbind(y[Length(y)]);
        od;

        count := 0;
        if Length(input_str) > 2 then
            for i in [1..Length(input_str) - 2] do
            if input_str[i] = '-' and input_str[i + 1] = '1' and
              \in(input_str[i+2],y) then
                if SIntChar(input_str[i + 2]) > 96 and
                  eps[SIntChar(input_str[i + 2]) - SIntChar('a') + 1] = 1 then
                    return false;
                elif SIntChar(input_str[i + 2]) < 96 and
                  sigma[SIntChar(input_str[i + 2]) - SIntChar('A') + 1] = 1 then
                    return false;
                else
                    for j in [i..Length(input_str) - 2] do
                        input_str[j] := input_str[j + 2];
                    od;
                    count := count + 2;
                fi;
                if input_str[i + 1] = '1' then input_str[i + 1] := '2';
                fi;
            fi;
            od;
        fi;

        if Length(input_str) > 1 then
        for i in [1..Length(input_str) - 1] do
            if input_str[i] = '1' and \in(input_str[i + 1], y) then
                if SIntChar(input_str[i + 1]) > 96 and
                  eps[SIntChar(input_str[i + 1]) - SIntChar('a') + 1] = -1 then
                    return false;
                elif SIntChar(input_str[i + 1]) < 96 and
                  sigma[SIntChar(input_str[i + 1]) - SIntChar('A') + 1]= -1 then
                    return false;
                else
                    for j in [i..Length(input_str) - 1] do
                        input_str[j] := input_str[j + 1];
                    od;
                    count := count + 1;
                fi;
            fi;
        od;
        fi;

        if Length(input_str) > 1 then
            if input_str[Length(input_str) - 1] = '-'
              and input_str[Length(input_str)] = '1' and
              \in(input_str[Length(input_str) - 2], y) then
                if SIntChar(input_str[Length(input_str) - 2]) > 96 and
                  sigma[SIntChar(input_str[Length(input_str) - 2]) -
                  SIntChar('a') + 1] = -1 then
                    return false;
                elif SIntChar(input_str[Length(input_str)-2]) < 96 and
                  eps[SIntChar(input_str[Length(input_str) - 2]) -
                  SIntChar('A') + 1] = -1 then
                    return false;
                else count := count + 2;
                fi;
            elif input_str[Length(input_str)] = '1' and
              \in(input_str[Length(input_str) - 1], y) then
                if SIntChar(input_str[Length(input_str) - 1]) > 96 and
                  sigma[SIntChar(input_str[Length(input_str) - 1]) -
                  SIntChar('a') + 1] = 1 then
                    return false;
                elif SIntChar(input_str[Length(input_str) - 1]) < 96 and
                  eps[SIntChar(input_str[Length(input_str) - 1]) -
                  SIntChar('A') + 1] = 1 then
                    return false;
                else
                    count := count + 1;
                fi;
            fi;
        fi;

        for i in [1..count] do
            Unbind(input_str[Length(input_str)]);
        od;

        for i in [1..Length(input_str)] do
            if not input_str[i] in y then
                return false;
            fi;
        od;

        if IsReducedWalk(Q, input_str) = 0 then return false;
        fi;

        if IsValidWalk(Q, input_str) = false then
            return false;
        fi;

        x := 0;
        for i in [1..Length(rho)] do
            temp := "";
            temp_caps := "";
            rel := WalkOfPath(rho[i]);
            for j in [1..Length(rel)] do
                temp := Concatenation(String(rel[j]), temp);
            od;
            for j in [1..Length(rel)] do
                temp1 := CharSInt(SIntChar(String(rel[j])[1]) - 32);
                temp_caps[j] := temp1;
            od;

            l := Length(input_str)-Length(temp);

            for j in [0..l] do
                for k in [1..Length(temp)] do
                    if input_str[j + k] <> temp[k] then
                        break;
                    fi;
                od;
                if k = Length(temp) and input_str[j + k] = temp[k] then
                    x := 3;
                    break;
                fi;
            od;
            for j in [0..l] do
                for k in [1..Length(temp_caps)] do
                    if input_str[j + k] <> temp_caps[k] then
                        break;
                    fi;
                od;
                if k = Length(temp_caps) and input_str[j +k] = temp_caps[k] then
                    x := 3;
                    break;
                fi;
            od;
        od;

        if x=0 then return true;
        else return false;
        fi;
    end;

    return IsValidStringSA(A, Q, rho, input_str);
end
);

InstallGlobalFunction( QPAStringDirectLeft,
  function(A, input_str, sigma, eps)
    local Q, x, y, z, temp, output, num, num_digits, j, temp3;
    output := IsValidString(A, input_str);
    #if output <> "Valid Positive Length String" then
    #       return output;
    #fi;
    Q := QuiverOfPathAlgebra(A);
    if output = true and input_str[1] = '(' then
        num := NumberOfVertices(Q);
        num_digits := 0;
        while num <> 0 do
            num := num - num mod 10;
            num := num/10;
            num_digits := num_digits + 1;
        od;
        num := 0;
        for j in [2..num_digits + 1] do
            num := num +
            (SIntChar(input_str[j]) - SIntChar('0')) * 10^(num_digits - j + 1);
        od;

        temp := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[num]);
        if input_str[num_digits + 3] = '1' then
            if Length(temp) = 1 and
              sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1 then
                z := String(temp[1]);
            elif Length(temp) = 2 then
                if sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1 then
                    z := String(temp[1]);
                else
                    z := String(temp[2]);
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        elif input_str[num_digits + 3] = '-' then
            if Length(temp) = 1 and
              sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1 then
                z := String(temp[1]);
            elif Length(temp) = 2 then
                if sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1 then
                    z := String(temp[1]);
                else
                    z := String(temp[2]);
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        fi;
    fi;

    if SIntChar(input_str[1]) > 96 then
        x := SIntChar(String(TargetOfPath(ArrowsOfQuiver(Q)[
        SIntChar(input_str[1]) - SIntChar('a') + 1]))[2]) - SIntChar('0');

        temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp3) = 2 then
            z := [String(temp3[1])[1]];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                z := [String(temp3[2])[1]];
                Append(z, input_str);
                if IsValidString(A, z) = false then
                    return "Cannot Perform The Operation";
                fi;
            fi;
        elif Length(temp3) = 1 then
            z := [String(temp3[1])[1]];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                return "Cannot Perform The Operation";
            fi;
        else return "Cannot Perform The Operation";
        fi;
    else
        x := SIntChar(String(SourceOfPath(ArrowsOfQuiver(Q)[
        SIntChar(input_str[1]) - SIntChar('A') + 1]))[2]) - SIntChar('0');

        temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp3) = 2 then
            if SIntChar(String(temp3[1])[1])
              - SIntChar(input_str[1]) = 32 then
                y := String(temp3[2])[1];
            else
                y := String(temp3[1])[1];
            fi;
            z := [y];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                return "Cannot Perform The Operation";
            fi;
        else
            return "Cannot Perform The Operation";
        fi;
    fi;
    return z;
end
);

InstallGlobalFunction( QPAStringInverseLeft,
  function(A, input_str, sigma, eps)
    local Q, x, y, z, temp, output, num, num_digits, j, temp2;
    output := IsValidString(A, input_str);
    #if output = "Valid Positive Length String" then
    #  return output;
    #fi;
    Q := QuiverOfPathAlgebra(A);
    if output = true and input_str[1] = '(' then
        num := NumberOfVertices(Q);
        num_digits := 0;
        while num <> 0 do
            num := num - num mod 10;
            num := num/10;
            num_digits := num_digits + 1;
        od;
        num := 0;
        for j in [2..num_digits + 1] do
            num := num +
            (SIntChar(input_str[j]) - SIntChar('0')) * 10^(num_digits - j + 1);
        od;
        temp := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[num]);
        if input_str[num_digits + 3] = '1' then
            if Length(temp) = 1 and eps[SIntChar(String(temp[1])[1])
              - SIntChar('a') + 1] = -1 then
                z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
            elif Length(temp) = 2 then
                if eps[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1 then
                    z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
                else
                    z := [CharInt(SIntChar(String(temp[2])[1]) - 32)];
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        elif input_str[num_digits + 3] = '-' then
            if Length(temp) = 1 and
              eps[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1 then
                z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
            elif Length(temp) = 2 then
                if eps[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1 then
                    z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
                else
                    z := [CharInt(SIntChar(String(temp[2])[1]) - 32)];
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        fi;
    fi;

    if SIntChar(input_str[1]) > 96 then
        x := SIntChar(String(TargetOfPath(ArrowsOfQuiver(Q)[
        SIntChar(input_str[1]) - SIntChar('a') + 1]))[2]) - SIntChar('0');

        temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp2) = 2 then
            if String(temp2[1])[1] = input_str[1] then
                y := String(temp2[2])[1];
            else
                y := String(temp2[1])[1];
            fi;
            z := [CharInt(SIntChar(y) - 32)];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                return "Cannot Perform The Operation";
            fi;
        else
            return "Cannot Perform The Operation";
        fi;
    else
        x := SIntChar(String(SourceOfPath(ArrowsOfQuiver(Q)[
        SIntChar(input_str[1]) - SIntChar('A') + 1]))[2]) - SIntChar('0');

        temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp2) = 2 then
            z := [CharInt(SIntChar(String(temp2[1])[1]) - 32)];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                z := [CharInt(SIntChar(String(temp2[2])[1]) - 32)];
                Append(z, input_str);
                if IsValidString(A, z) = false then
                    return "Cannot Perform The Operation";
                fi;
            fi;
        elif Length(temp2) = 1 then
            z := [CharInt(SIntChar(String(temp2[1])[1]) - 32)];
            Append(z, input_str);
            if IsValidString(A, z) = false then
                return "Cannot Perform The Operation";
            fi;
        else return "Cannot Perform The Operation";
        fi;
    fi;
    return z;
end
);

InstallGlobalFunction( QPAStringDirectRight,
  function(A, input_str, sigma, eps)
    local Q, x, y, z, temp, output, i, num, num_digits, j, temp2;
    #if IsValidString(Q,rho,input_str) = false then
    #return false;
    #fi;
    output := IsValidString(A, input_str);
    Q := QuiverOfPathAlgebra(A);

    if output = true and input_str[1] = '(' then
        num := NumberOfVertices(Q);
        num_digits := 0;
        while num <> 0 do
            num := num - num mod 10;
            num := num/10;
            num_digits := num_digits + 1;
        od;
        num := 0;
        for j in [2..num_digits+1] do
            num := num +
            (SIntChar(input_str[j])-SIntChar('0')) * 10^(num_digits -j + 1);
        od;
        temp := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[num]);
        if input_str[num_digits + 3] = '1' then
            if Length(temp) = 1 and eps[SIntChar(String(temp[1])[1])
              - SIntChar('a')+1] = 1 then
                z := temp;
            elif Length(temp) = 2 then
                if eps[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1
                  then z := String(temp[1]);
                else
                    z := String(temp[2]);
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        elif input_str[num_digits + 3] = '-' then
            if Length(temp) = 1 and eps[SIntChar(String(temp[1])[1])
              - SIntChar('a') + 1] = -1 then
                z := String(temp[1]);
            elif Length(temp) = 2 then
                if eps[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1
                  then z := String(temp[1]);
                else
                    z := String(temp[2]);
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        fi;
    fi;
    #temp := input_str;
    temp := [];
    for i in [1..Length(input_str)] do
        Add(temp, input_str[i]);
    od;
    if SIntChar(temp[Length(temp)]) > 96 then
        x := SIntChar(String(SourceOfPath(ArrowsOfQuiver(Q)[
        SIntChar(temp[Length(temp)]) - SIntChar('a') + 1]))[2])
        - SIntChar('0');
        temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[x]);
        if Length(temp2) = 2 then
            z := [String(temp2[1])[1]];
            temp := Concatenation(temp, z);
            if IsValidString(A, temp) = false then
                z := [String(temp2[2])[1]];
                temp := input_str;
                temp := Concatenation(temp, z);
                if IsValidString(A, temp) = false then
                    return "Cannot Perform The Operation";
                fi;
            fi;
        elif Length(temp2) = 1 then
            z := [String(temp2[1])[1]];
            temp := Concatenation(temp, z);
            if IsValidString(A, temp) = false then
                return "Cannot Perform The Operation";
            fi;
        else return "Cannot Perform The Operation";
        fi;
    else
        x := SIntChar(String(TargetOfPath(ArrowsOfQuiver(Q)[
        SIntChar(temp[Length(temp)]) - SIntChar('A') + 1]))[2])
        - SIntChar('0');
        temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[x]);
        if Length(temp2) = 2 then
            if SIntChar(String(temp2[1])[1])
              - SIntChar(temp[Length(temp)]) = 32 then
                y := String(temp2[2])[1];
            else
                y := String(temp2[1])[1];
            fi;
            z := [y];
            temp := Concatenation(temp, z);
            if IsValidString(A, temp) = false then
                return "Cannot Perform The Operation";
            fi;
        else
            return "Cannot Perform The Operation";
        fi;
    fi;
    return temp;
end
);

InstallGlobalFunction( QPAStringInverseRight,
  function(A, input_str, sigma, eps)
    local Q, x, y, z, temp, output, i, num_digits, num, j, temp3;
    #if IsValidString(Q,rho,input_str) = false then
    #  return false;
    #fi;
    output := IsValidString(A, input_str);
    Q := QuiverOfPathAlgebra(A);

    if output = true and input_str[1] = '(' then
        num := NumberOfVertices(Q);
        num_digits := 0;
        while num <> 0 do
            num := num - num mod 10;
            num := num/10;
            num_digits := num_digits + 1;
        od;
        num := 0;
        for j in [2..num_digits + 1] do
            num := num +
            (SIntChar(input_str[j]) - SIntChar('0')) * 10^(num_digits - j + 1);
        od;
        temp := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[num]);
        if input_str[num_digits + 3] = '1' then
            if Length(temp) = 1 and sigma[SIntChar(String(temp[1])[1])
              - SIntChar('a') + 1] = 1 then
                z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
            elif Length(temp) = 2 then
                if sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = 1 then
                    z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
                else
                    z := [CharInt(SIntChar(String(temp[2])[1]) - 32)];
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        elif input_str[num_digits + 3] = '-' then
            if Length(temp) = 1 and
              sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1 then
                z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
            elif Length(temp) = 2 then
                if sigma[SIntChar(String(temp[1])[1]) - SIntChar('a') + 1] = -1 then
                    z := [CharInt(SIntChar(String(temp[1])[1]) - 32)];
                else
                    z := [CharInt(SIntChar(String(temp[2])[1]) - 32)];
                fi;
            else z := "Cannot Perform The Operation";
            fi;
            return z;
        fi;
    fi;

    #temp1 := input_str;
    temp := [];
    for i in [1..Length(input_str)] do
        Add(temp,input_str[i]);
    od;
    #temp := input_str;
    #input_str := temp1;
    if SIntChar(temp[Length(temp)]) > 96 then
        x := SIntChar(String(SourceOfPath(ArrowsOfQuiver(Q)[
        SIntChar(temp[Length(temp)]) - SIntChar('a') + 1]))[2]) - SIntChar('0');

        temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp3) = 2 then
            if String(temp3[1])[1] = temp[Length(temp)] then
                y := String(temp3[2])[1];
            else
                y := String(temp3[1])[1];
            fi;
            z := [CharInt(SIntChar(y) - 32)];
            Append(temp, z);
            if IsValidString(A, temp) = false then
                return "Cannot Perform The Operation";
            fi;
        else
            return "Cannot Perform The Operation";
        fi;
    else
        x := SIntChar(String(TargetOfPath(ArrowsOfQuiver(Q)[
        SIntChar(temp[Length(temp)]) - SIntChar('A') + 1]))[2]) - SIntChar('0');

        temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[x]);

        if Length(temp3) = 2 then
            z := [CharInt(SIntChar(String(temp3[1])[1]) - 32)];
            temp := Concatenation(temp, z);
            if IsValidString(A, temp) = false then
                temp := input_str;
                z := [CharInt(SIntChar(String(temp3[2])[1]) - 32)];
                Append(temp, z);
                if IsValidString(A, temp) = false then
                    return "Cannot Perform The Operation";
                fi;
            fi;
        elif Length(temp3) = 1 then
            z := [CharInt(SIntChar(String(temp3[1])[1]) - 32)];
            Append(temp, z);
            if IsValidString(A, temp) = false then
                return "Cannot Perform The Operation";
            fi;
        else return "Cannot Perform The Operation";
        fi;
    fi;
    return temp;
end
);

InstallMethod( LocalARQuiver,
"for stringalgebras and strings",
true,
[ IsAlgebra, IsString ], 0,
function( A, input_str )
    local Q, rho1, rho, vertices, arrows, onerel, tworel, str, arr, temp, temp1,
    temp2, temp3, La, Lb, Ra, Rb, array,
    sigma, eps, k, kQ, rel, arr1, gens, i, j;

    if IsStringAlgebra(A) = false then
      Error("The input algebra is not a string algebra");
      return 0;
    fi;

    Q := QuiverOfPathAlgebra(A);
    rho1 := RelationsOfAlgebra(A);

    vertices := VerticesOfQuiver(Q);
    arrows := ArrowsOfQuiver(Q);

    for i in [1..Length(vertices)] do
      if Length(String(vertices[i])) <> 2 or String(vertices[i])[1] <> 'v' or
        String(vertices[i])[2] <> String(i)[1] then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    for i in [1..Length(arrows)] do
      if Length(String(arrows[i])) <> 1 or String(arrows[i])[1] <> CharInt(i+96) then
        Error("The Quiver of A is not acceptable");
        return 0;
      fi;
    od;

    rho := [];
    for i in [1..Length(rho1)] do
        onerel := CoefficientsAndMagmaElements(rho1[i]);
        tworel := WalkOfPath(onerel[1]);
        str := Concatenation(List(tworel, a -> String(a)));
        Append(rho, [str]);
    od;

    temp := IsValidString(A, input_str);
    if temp = false then
        Error("The input string is invalid");
    fi;

    array := QPAStringSigmaEps(A);
    sigma := array[1];
    eps := array[2];

    La := function(A, input_str)
        local x,y;
        #output := IsValidString(Q,rho,input_str);
        #if output <> "Valid Positive Length String" and output <> "Valid Zero Length String" then return output;
        #fi;
        x := QPAStringInverseLeft(A, input_str, sigma, eps);
        if x = "Cannot Perform The Operation" then return 0;
        fi;
        while true do
            y := QPAStringDirectLeft(A, x, sigma, eps);
            if y = "Cannot Perform The Operation" then return x;
            else x := y;
            fi;
        od;
    end;

    Lb := function(A, input_str)
        local x, y;
        x := QPAStringDirectLeft(A, input_str, sigma, eps);
        if x = "Cannot Perform The Operation" then return 0;
        fi;
        while true do
            y := QPAStringInverseLeft(A, x, sigma, eps);
            if y = "Cannot Perform The Operation" then return x;
            else x := y;
            fi;
        od;
    end;

    Ra := function(A, input_str)
        local x, y;
        x := QPAStringInverseRight(A, input_str, sigma, eps);
        if x = "Cannot Perform The Operation" then return 0;
        fi;
        while true do
            y := QPAStringDirectRight(A, x, sigma, eps);
            if y = "Cannot Perform The Operation" then return x;
            else x := y;
            fi;
        od;
    end;

    Rb := function(A, input_str)
        local x, y;
        x := QPAStringDirectRight(A, input_str, sigma, eps);
        if x = "Cannot Perform The Operation" then return 0;
        fi;
        while true do
            y := QPAStringInverseRight(A, x, sigma, eps);
            if y = "Cannot Perform The Operation" then return x;
            else x := y;
            fi;
        od;
    end;

    arr := [];
    temp1 := La(A, input_str);
    if temp1 <> 0 then
        Append(arr, [Concatenation
        ("The canonical embedding StringModule(", input_str, ") to StringModule(", temp1, ")")]);
    fi;

    temp1 := Lb(A, input_str);
    if temp1 <> 0 then
        Append(arr, [Concatenation
        ("The canonical projection StringModule(", temp1, ") to StringModule(", input_str, ")")]);
    fi;

    temp1 := Ra(A, input_str);
    if temp1 <> 0 then
        Append(arr, [Concatenation
        ("The canonical projection StringModule(", temp1, ") to StringModule(", input_str, ")")]);
    fi;

    temp1 := Rb(A, input_str);
    if temp1 <> 0 then
        Append(arr, [Concatenation
        ("The canonical embedding StringModule(", input_str, ") to StringModule(", temp1, ")")]);
    fi;
    return arr;
end
);
