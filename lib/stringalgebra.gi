InstallMethod( IsABand,
"for quivers and relations and strings",
true,
[ IsQuiver, IsList, IsString ], 0,
function ( Q, rho, input_str )
    local computeLPSarray, SourceOfArrow, TargetOfArrow, output,
    w1, w2, n, i, j, k, arr, m_arr, small, capital, count, temp;

    computeLPSarray := function(pat)
        local lps, i, len, m;
        m := Length(pat);
        lps := [];
        for i in [1..m] do
            Append(lps, [1]);
        od;
        len := 1;
        i := 2;
        while i <= m do
            if pat[i] = pat[len] then
                len := len + 1;
                lps[i] := len;
                i := i+1;
            else
                if len <> 1 then
                    len := lps[len - 1];
                else
                    lps[i] := 1;
                    i := i + 1;
                fi;
            fi;
        od;
        for i in [1..Length(lps)] do
            lps[i] := lps[i] - 1;
        od;
        return lps;
    end;

    SourceOfArrow := function(Q, ch)
        local i, ch1;
        ch1 := ch;
        if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            ch := CharInt(SIntChar(ch) + 32);
        fi;
        if SIntChar(ch) - SIntChar('a') + 1 > NumberOfArrows(Q) then
            return "Invalid String";
        fi;
        for i in [1..NumberOfArrows(Q)] do
            if String(ArrowsOfQuiver(Q)[i])[1] = ch then
                break;
            fi;
        od;
        if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            return TargetOfPath(ArrowsOfQuiver(Q)[i]);
        else
            return SourceOfPath(ArrowsOfQuiver(Q)[i]);
        fi;
    end;

    TargetOfArrow := function(Q, ch)
        local i, ch1;
        ch1 := ch;
        if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            ch := CharInt(SIntChar(ch) + 32);
        fi;
        if SIntChar(ch) - SIntChar('a') + 1 > NumberOfArrows(Q) then
            return "Invalid String";
        fi;
        for i in [1..NumberOfArrows(Q)] do
            if String(ArrowsOfQuiver(Q)[i])[1] = ch then
                break;
            fi;
        od;
        if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            return SourceOfPath(ArrowsOfQuiver(Q)[i]);
        else
            return TargetOfPath(ArrowsOfQuiver(Q)[i]);
        fi;
    end;

    output := IsValidString(Q, rho, input_str);
    if output <> "Valid Positive Length String" and output <> "Valid Zero Length String" then
        return output;
    fi;

    if output = "Valid Zero Length String" then return 0;
    fi;
    if Length(input_str) = 1 then return 0;
    fi;

    n := Length(input_str);
    w1 := SourceOfArrow(Q, input_str[n]);
    w2 := TargetOfArrow(Q, input_str[1]);

    if w1 <> w2 then return 0;
    fi;

    capital := 0;
    small := 0;

    for i in [1..Length(input_str)] do
        if SIntChar(input_str[i]) >= SIntChar('A') and
          SIntChar(input_str[i]) <= SIntChar('Z') then
            capital := 1;
        elif SIntChar(input_str[i]) >= SIntChar('a') and
          SIntChar(input_str[i]) <= SIntChar('z') then
            small := 1;
        fi;
    od;
    if capital = 0 or small = 0 then return 0;
    fi;
    temp := Concatenation(input_str, input_str);
    if IsValidString(Q, rho, temp) <> "Valid Positive Length String" then
        return 0;
    fi;
    arr := computeLPSarray(input_str);
    if arr[Length(arr)] <> 0 and Length(input_str) mod (Length(input_str) - arr[Length(arr)]) = 0 then
        return 0;
    else
        return 1;
    fi;
end
);

InstallMethod( StringsLessThan,
"for quivers and relations and levels",
true,
[ IsQuiver, IsList, IsInt ], 0,
function ( Q, rho, level )
    local QuiverRho, making_tree, direct_left,
    inverse_left, k, tree, array, i, vertex, sigma_eps;

    direct_left := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, temp3, output, num, num_digits, j;
        output := IsValidString(Q, rho, input_str);
        #if output <> "Valid Positive Length String" then
        #       return output;
        #fi;
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    z := [String(temp3[2])[1]];
                    Append(z, input_str);
                    if IsValidString(Q, rho, z) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp3) = 1 then
                z := [String(temp3[1])[1]];
                Append(z, input_str);
                if IsValidString(Q,rho,z) = "Invalid String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else
                return "Cannot Perform The Operation";
            fi;
        fi;
        return z;
    end;

    inverse_left := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, output, num, num_digits, j, temp2, temp3;
        output := IsValidString(Q, rho, input_str);
        #if output = "Valid Positive Length String" then
        #  return output;
        #fi;
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    z := [CharInt(SIntChar(String(temp2[2])[1]) - 32)];
                    Append(z, input_str);
                    if IsValidString(Q, rho, z) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp2) = 1 then
                z := [CharInt(SIntChar(String(temp2[1])[1]) - 32)];
                Append(z, input_str);
                if IsValidString(Q, rho, z) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else return "Cannot Perform The Operation";
            fi;
        fi;
        return z;
    end;

    sigma_eps := function( Q, rho1 )
        local QuiverRho, Composable_Matrix, IsComposable, temp_in, temp_out, i,
        j, k, sigma, eps, rho, temp2, temp3;

        QuiverRho := function(Q,rho)
            local rho_q, i, j, k, arr, arr1;
            rho_q := [];
            for i in [1..Length(rho)] do
                arr := [];
                for j in [1..Length(rho[i])] do
                    for k in [1..NumberOfArrows(Q)] do
                        if rho[i][j] = String(ArrowsOfQuiver(Q)[k])[1] then
                            Add(arr, ArrowsOfQuiver(Q)[k]);
                        fi;
                    od;
                od;
                arr1 := arr[1];
                for k in [2..Length(arr)] do
                    arr1 := arr1 * arr[k];
                od;
                Add(rho_q, arr1);
            od;
            return rho_q;
        end;

        Composable_Matrix := function(Q,rho)
            local mat, x, str, i, j;
            mat := [];
            for i in [1..NumberOfArrows(Q)] do
                str := [];
                for j in [1..NumberOfArrows(Q)] do
                    if SourceOfPath(ArrowsOfQuiver(Q)[i]) =
                      TargetOfPath(ArrowsOfQuiver(Q)[j]) and not
                      \in(ArrowsOfQuiver(Q)[j] * ArrowsOfQuiver(Q)[i], rho) then
                        Append(str, [1]);
                    else
                        Append(str, [0]);
                    fi;
                od;
                Append(mat, [str]);
            od;
            return mat;
        end;

        IsComposable := function(Q, rho, arrow1, arrow2)
            local i, Mat, pos1, pos2, j;
            Mat := Composable_Matrix(Q, rho);
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

        temp_in := function(Q, i, j)
            return SIntChar(String(IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        temp_out := function(Q, i, j)
            return SIntChar(String(OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        sigma := [];
        eps := [];
        for i in [1..NumberOfArrows(Q)] do
            Append(sigma, [0]);
            Append(eps, [0]);
        od;

        rho := QuiverRho(Q, rho1);

        for i in [1..NumberOfVertices(Q)] do
            temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            if Length(temp2) = 2 then
                eps[temp_in(Q, i, 1)] := 1;
                eps[temp_in(Q, i, 2)] := -1;
                if Length(temp3) = 1 then
                    if IsComposable(Q, rho, temp3[1], temp2[1]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                    elif IsComposable(Q, rho, temp3[1], temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := 1;
                    fi;
                elif Length(temp3) = 2 then
                    if IsComposable(Q, rho, temp3[1], temp2[1]) = 1 or
                      IsComposable(Q, rho, temp3[2], temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                        sigma[temp_out(Q, i, 2)] := 1;
                    else
                        sigma[temp_out(Q, i, 1)] := 1;
                        sigma[temp_out(Q, i, 2)] := -1;
                    fi;
                fi;
            elif Length(temp2) = 1 then
                if Length(temp3) = 1 then
                    eps[temp_in(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 1)] := -1;
                elif Length(temp3) = 2 then
                    sigma[temp_out(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 2)] := -1;
                    if IsComposable(Q, rho, temp3[1],
                      temp2[1]) = 1 then
                        eps[temp_in(Q, i, 1)] := -1;
                    else
                        eps[temp_in(Q, i, 1)] := 1;
                    fi;
                fi;
            elif Length(temp2) = 0 then
              sigma[temp_out(Q, i, 1)] := 1;
              sigma[temp_out(Q, i, 2)] := -1;
            fi;
        od;
        return [sigma, eps];
    end;

    making_tree := function(Q, rho, vertex, level)
        local tree, x, i, count, j, array, sigma,eps;
        array := sigma_eps(Q,rho);
        sigma := array[1];
        eps := array[2];
        tree := [];
        tree[1] := [1, vertex];
        i := 1;
        count := 0;
        x := direct_left(Q, rho, vertex, sigma, eps);
        if x <> "Cannot Perform The Operation" then
            Add(tree, [2 * 1, x]);
            j := 2;
        else
            count := count + 1;
        fi;

        x := inverse_left(Q, rho, vertex, sigma, eps);
        if x <> "Cannot Perform The Operation" then
            Add(tree, [2 * 1 + 1, x]);
            j := 3;
        else
            count := count + 1;
        fi;

        if count = 2 then return [];
        fi;
        i := 2;
        while true do
            if 2*tree[i][1] > 2^(level+1) - 1 then break;
            fi;
            x := direct_left(Q, rho, tree[i][2], sigma, eps);
            if x <> "Cannot Perform The Operation" then
                Add(tree, [2 * tree[i][1], x]);
                j := 2*tree[i][1];
            fi;

            x := inverse_left(Q, rho, tree[i][2], sigma, eps);
            if x <> "Cannot Perform The Operation" then
      			Add(tree, [2 * tree[i][1] + 1, x]);
      			j := 2 * tree[i][1] + 1;
            fi;
            i := i + 1;
        od;
        return tree;
    end;

    QuiverRho := function(Q,rho)
        local kQ, i, j, k, gens, arr, arr1, rel;
        k := Rationals;
        kQ := PathAlgebra(k,Q);
        gens := GeneratorsOfAlgebra(kQ);
        rel := [];
        for i in [1..Length(rho)] do
            arr := [];
            for j in [1..Length(rho[i])] do
                for k in [NumberOfVertices(Q) + 1..NumberOfVertices(Q)+NumberOfArrows(Q)] do
                    if rho[i][j] = String(gens[k])[5] then
                        Add(arr, gens[k]);
                        break;
                    fi;
                od;
            od;
            arr1 := arr[1];
            for k in [2..Length(arr)] do
                arr1 := arr1 * arr[k];
            od;
            Add(rel, arr1);
        od;
        if IsStringAlgebra(kQ/rel) = false then
            return "Not a String Algebra";
        else return 1;
        fi;
    end;

    if QuiverRho(Q,rho) = "Not a String Algebra" then
        return "Not a String Algebra";
    fi;

    array := [];
    if level < 0 then return array;
    elif level = 0 then
        for k in [1..NumberOfVertices(Q)] do
            vertex := Concatenation("(", String(k), ",", String(1), ")");
            Add(array, vertex);
            vertex := Concatenation("(", String(k), ",", String(-1), ")");
            Add(array, vertex);
        od;
        return array;
    fi;
    for k in [1..NumberOfVertices(Q)] do
        vertex := Concatenation("(", String(k), ",", String(1), ")");
        tree := making_tree(Q, rho, vertex, level);
        for i in [1..Length(tree)] do
            Add(array, tree[i][2]);
        od;
        vertex := Concatenation("(", String(k), ",", String(-1), ")");
        tree := making_tree(Q, rho, vertex, level);
        for i in [1..Length(tree)] do
            Add(array, tree[i][2]);
        od;
    od;
    return array;
end
);

InstallMethod( BandsLessThan,
"for quivers and relations and levels",
true,
[ IsQuiver, IsList, IsInt ], 0,
function ( Q, rho, level )
    local q, b, i;
    q := StringsLessThan(Q, rho, level);
    if q = "Not a String Algebra" then return q;
    fi;
    b := [];
    for i in [1..Length(q)] do
        if IsABand(Q, rho, q[i]) = 1 then
            Append(b, [q[i]]);
        fi;
    od;
    return b;
end
);

InstallMethod( BandRepresentativesLessThan,
"for quivers and relations and levels",
true,
[ IsQuiver, IsList, IsInt ], 0,
function ( Q, rho, level )
    local CyclicPermutationOfABand, b, tempo, br, i;
    CyclicPermutationOfABand := function(Q, rho, band)
        local array, temp, i, j;
        array := [];
        for i in [1..Length(band)] do
            Append(array, [0]);
        od;
        for i in [1..Length(band)] do
            temp := band[1];
            for j in [1..Length(band) - 1] do
                band[j] := band[j + 1];
            od;
            band[j + 1] := temp;
            array[i] := Concatenation(band,"");
        od;
        return array;
    end;

    b := BandsLessThan(Q,rho,level);
    if b = "Not a String Algebra" then return b;
    fi;
    tempo := [];
    br := [];
    for i in [1..Length(b)] do
        if \in(b[i], tempo) = false and SIntChar(b[i][1]) >= 97 and
          SIntChar(b[i][Length(b[i])]) < 97 then
            Append(br, [b[i]]);
            Append(tempo, CyclicPermutationOfABand(Q, rho, b[i]));
        fi;
    od;
    return br;
end
);

InstallMethod( IsDomesticStringAlgebra,
"for quivers and relations",
true,
[ IsQuiver, IsList], 0,
function( Q, rho )
    local sigma_eps, direct_left,
    inverse_right, NumberOfJoints, MaximumLengthOFDirectString, ToVerifyStrAlg,
    k, kQ, gens, rel, arr1, arr, i, j, N, L, level, b, count, arr_quiver, len, temp;

    k := Rationals;
    kQ := PathAlgebra(k,Q);
    gens := GeneratorsOfAlgebra(kQ);
    rel := [];
    for i in [1..Length(rho)] do
        arr := [];
        for j in [1..Length(rho[i])] do
            for k in [NumberOfVertices(Q) + 1..NumberOfVertices(Q) +
              NumberOfArrows(Q)] do
                if rho[i][j] = String(gens[k])[5] then
                    Add(arr, gens[k]);
                    break;
                fi;
            od;
        od;
        arr1 := arr[1];
        for k in [2..Length(arr)] do
            arr1 := arr1 * arr[k];
        od;
        Add(rel, arr1);
    od;
    if IsStringAlgebra(kQ/rel) = false then
      return "Not a String Algebra";
    fi;

    direct_left := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, temp3, output, num, num_digits, j;
        output := IsValidString(Q, rho, input_str);
        #if output <> "Valid Positive Length String" then
        #       return output;
        #fi;
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    z := [String(temp3[2])[1]];
                    Append(z, input_str);
                    if IsValidString(Q, rho, z) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp3) = 1 then
                z := [String(temp3[1])[1]];
                Append(z, input_str);
                if IsValidString(Q,rho,z) = "Invalid String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else
                return "Cannot Perform The Operation";
            fi;
        fi;
        return z;
    end;

    inverse_right := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, temp3, output, i, num_digits, num, j;
        #if IsValidString(Q,rho,input_str) = "Invalid String" then
        #  return "Invalid String";
        #fi;
        output := IsValidString(Q, rho, input_str);
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, temp) = "Invalid String" then
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
                if IsValidString(Q, rho, temp) = "Invalid String" then
                    temp := input_str;
                    z := [CharInt(SIntChar(String(temp3[2])[1]) - 32)];
                    Append(temp, z);
                    if IsValidString(Q, rho, temp) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp3) = 1 then
                z := [CharInt(SIntChar(String(temp3[1])[1]) - 32)];
                Append(temp, z);
                if IsValidString(Q, rho, temp) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else return "Cannot Perform The Operation";
            fi;
        fi;
      return temp;
    end;

    sigma_eps := function( Q, rho1 )
        local QuiverRho, Composable_Matrix, IsComposable, temp_in, temp_out, i,
        j, k, sigma, eps, rho, temp2, temp3;

        QuiverRho := function(Q,rho)
            local rho_q, i, j, k, arr, arr1;
            rho_q := [];
            for i in [1..Length(rho)] do
                arr := [];
                for j in [1..Length(rho[i])] do
                    for k in [1..NumberOfArrows(Q)] do
                        if rho[i][j] = String(ArrowsOfQuiver(Q)[k])[1] then
                            Add(arr, ArrowsOfQuiver(Q)[k]);
                        fi;
                    od;
                od;
                arr1 := arr[1];
                for k in [2..Length(arr)] do
                    arr1 := arr1 * arr[k];
                od;
                Add(rho_q, arr1);
            od;
            return rho_q;
        end;

        Composable_Matrix := function(Q,rho)
            local mat, x, str, i, j;
            mat := [];
            for i in [1..NumberOfArrows(Q)] do
                str := [];
                for j in [1..NumberOfArrows(Q)] do
                    if SourceOfPath(ArrowsOfQuiver(Q)[i]) =
                      TargetOfPath(ArrowsOfQuiver(Q)[j]) and not
                      \in(ArrowsOfQuiver(Q)[j] * ArrowsOfQuiver(Q)[i], rho) then
                        Append(str, [1]);
                    else
                        Append(str, [0]);
                    fi;
                od;
                Append(mat, [str]);
            od;
            return mat;
        end;

        IsComposable := function(Q, rho, arrow1, arrow2)
            local i, Mat, pos1, pos2, j;
            Mat := Composable_Matrix(Q, rho);
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

        temp_in := function(Q, i, j)
            return SIntChar(String(IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        temp_out := function(Q, i, j)
            return SIntChar(String(OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        sigma := [];
        eps := [];
        for i in [1..NumberOfArrows(Q)] do
            Append(sigma, [0]);
            Append(eps, [0]);
        od;

        rho := QuiverRho(Q, rho1);

        for i in [1..NumberOfVertices(Q)] do
            temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            if Length(temp2) = 2 then
                eps[temp_in(Q, i, 1)] := 1;
                eps[temp_in(Q, i, 2)] := -1;
                if Length(temp3) = 1 then
                    if IsComposable(Q, rho, temp3[1],
                      temp2[1]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                    elif IsComposable(Q, rho, temp3[1], temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := 1;
                    fi;
                elif Length(temp3) = 2 then
                    if IsComposable(Q, rho, temp3[1], temp2[1]) = 1 or
                      IsComposable(Q, rho, temp3[2],
                      temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                        sigma[temp_out(Q, i, 2)] := 1;
                    else
                        sigma[temp_out(Q, i, 1)] := 1;
                        sigma[temp_out(Q, i, 2)] := -1;
                    fi;
                fi;
            elif Length(temp2) = 1 then
                if Length(temp3) = 1 then
                    eps[temp_in(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 1)] := -1;
                elif Length(temp3) = 2 then
                    sigma[temp_out(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 2)] := -1;
                    if IsComposable(Q, rho, temp3[1],
                      temp2[1]) = 1 then
                        eps[temp_in(Q, i, 1)] := -1;
                    else
                        eps[temp_in(Q, i, 1)] := 1;
                    fi;
                fi;
            elif Length(temp2) = 0 then
              sigma[temp_out(Q, i, 1)] := 1;
              sigma[temp_out(Q, i, 2)] := -1;
            fi;
        od;
        return [sigma, eps];
    end;

    MaximumLengthOFDirectString := function(Q,rho)
        local arr, i, max, x, y, array, sigma, eps;
        array := sigma_eps(Q, rho);
        sigma := array[1];
        eps := array[2];
        arr := [];
        for i in [1..NumberOfArrows(Q)] do
            x := [CharInt(i + 96)];
            while true do
                y := direct_left(Q, rho, x, sigma, eps);
                if y = "Cannot Perform The Operation" then
                    Append(arr, [Length(x)]);
                    break;
                else x := y;
                fi;
            od;
        od;
        max := arr[1];
        for i in [2..Length(arr)] do
            if arr[i] > max then max := arr[i];
            fi;
        od;
        return max;
    end;

    NumberOfJoints := function(Q,rho)
        local arr, i, temp, array, sigma, eps;
        arr := [];
        array := sigma_eps(Q, rho);
        sigma := array[1];
        eps := array[2];
        for i in [1..NumberOfArrows(Q)] do
            temp := inverse_right(Q, rho, [CharInt(i+96)], sigma, eps);
            if temp <> "Cannot Perform The Operation" then
                Append(arr, [temp]);
            fi;
        od;
        return Length(arr);
    end;

    N := NumberOfJoints(Q, rho);
    L := MaximumLengthOFDirectString(Q, rho);
    level := 2 * L * N;
    b := BandsLessThan(Q, rho, level);

    count := 0;
    arr_quiver := [];
    len := Length(b);
    for i in [1..len] do
        for j in [1..len] do
            if i <> j then
                temp := Concatenation(b[i], b[j]);
                if IsValidString(Q, rho, temp) = "Valid Positive Length String" then
                    count := count + 1;
                    Append(arr_quiver, [[j, i, Concatenation("a", String(count))]]);
                fi;
            fi;
        od;
    od;
    return IsAcyclicQuiver(Quiver(len, arr_quiver));
end
);

InstallMethod( BridgeQuiver,
"for quivers and relations",
true,
[ IsQuiver, IsList], 0,
function( Q, rho )
    local sigma_eps, arr, arr1,kQ, gens,rel, SourceOfArrow, TargetOfArrow,
    KMPSearch, computeLPSarray,Q2, CyclicPermutationOfABand, NumberOfJoints,
    MaximumLengthOFDirectString, inverse_right, direct_left, WeakBridgeQuiver,
    BridgeQuiver, BandFreeStrings, IsBandFreeString, wb, wb1, lambda, i, j, k,
    str, temp, array, array1, perm, Q1, num, l, N, L, level, q, b, br, tempo,
    count, arr_quiver, len;

    k := Rationals;
    kQ := PathAlgebra(k, Q);
    gens := GeneratorsOfAlgebra(kQ);
    rel := [];
    for i in [1..Length(rho)] do
        arr := [];
        for j in [1..Length(rho[i])] do
            for k in [NumberOfVertices(Q) + 1..NumberOfVertices(Q) +
              NumberOfArrows(Q)] do
                if rho[i][j] = String(gens[k])[5] then
                    Add(arr,gens[k]);
                    break;
                fi;
            od;
        od;
        arr1 := arr[1];
        for k in [2..Length(arr)] do
            arr1 := arr1*arr[k];
        od;
        Add(rel, arr1);
    od;
    if IsStringAlgebra(kQ/rel) = false then
      return "Not a String Algebra";
    fi;

    IsBandFreeString := function(Q, rho, input_str)
        local i, j, k, temp;
        #var := IsValidString(Q,rho,input_str);
        #if var <> "Valid Positive Length String" then
        #  return var;
        #fi;
        if input_str[1] = '(' then return 1;
        fi;

        for i in [1..Length(input_str)] do
            for j in [i + 1..Length(input_str)] do
                if TargetOfArrow(Q, input_str[i]) = SourceOfArrow(Q, input_str[j]) then
                    temp := [];
                    for k in [i..j] do Add(temp,input_str[k]);
                    od;
                    if IsABand(Q, rho, temp) = 1 then return 0;
                    fi;
                fi;
            od;
        od;

        return 1;
    end;

    BandFreeStrings := function(Q, rho, v1, v2, q)
        local tree, tree1, i, j, k, array, count, vertex, x;
        array := [];
        if v1 > NumberOfVertices(Q) or v2 > NumberOfVertices(Q) then
            return array;
        fi;

        for i in [1..Length(q)] do
            if q[i][1] = '(' then
                if v1 = v2 and q[i][2] = String(v1)[1] then
                    Add(array, q[i]);
                fi;
            else
                if String(TargetOfArrow(Q, q[i][1])) = Concatenation("v", String(v2))
                  and String(SourceOfArrow(Q, q[i][Length(q[i])])) =
                  Concatenation("v", String(v1)) then
                    if IsBandFreeString(Q, rho, q[i]) = 1 then
                        Add(array, q[i]);
                    fi;
                fi;
            fi;
        od;
        return array;
    end;

    direct_left := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, temp3, output, num, num_digits, j;
        output := IsValidString(Q, rho, input_str);
        #if output <> "Valid Positive Length String" then
        #       return output;
        #fi;
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    z := [String(temp3[2])[1]];
                    Append(z, input_str);
                    if IsValidString(Q, rho, z) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp3) = 1 then
                z := [String(temp3[1])[1]];
                Append(z, input_str);
                if IsValidString(Q,rho,z) = "Invalid String" then
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
                if IsValidString(Q, rho, z) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else
                return "Cannot Perform The Operation";
            fi;
        fi;
        return z;
    end;

    inverse_right := function(Q, rho, input_str, sigma, eps)
        local x, y, z, temp, temp3, output, i, num_digits, num, j;
        #if IsValidString(Q,rho,input_str) = "Invalid String" then
        #  return "Invalid String";
        #fi;
        output := IsValidString(Q, rho, input_str);
        if output = "Valid Zero Length String" then
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
                if IsValidString(Q, rho, temp) = "Invalid String" then
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
                if IsValidString(Q, rho, temp) = "Invalid String" then
                    temp := input_str;
                    z := [CharInt(SIntChar(String(temp3[2])[1]) - 32)];
                    Append(temp, z);
                    if IsValidString(Q, rho, temp) = "Invalid String" then
                        return "Cannot Perform The Operation";
                    fi;
                fi;
            elif Length(temp3) = 1 then
                z := [CharInt(SIntChar(String(temp3[1])[1]) - 32)];
                Append(temp, z);
                if IsValidString(Q, rho, temp) = "Invalid String" then
                    return "Cannot Perform The Operation";
                fi;
            else return "Cannot Perform The Operation";
            fi;
        fi;
      return temp;
    end;

    sigma_eps := function( Q, rho1 )
        local QuiverRho, Composable_Matrix, IsComposable, temp_in, temp_out, i,
        j, k, sigma, eps, rho, temp2, temp3;

        QuiverRho := function(Q,rho)
            local rho_q, i, j, k, arr, arr1;
            rho_q := [];
            for i in [1..Length(rho)] do
                arr := [];
                for j in [1..Length(rho[i])] do
                    for k in [1..NumberOfArrows(Q)] do
                        if rho[i][j] = String(ArrowsOfQuiver(Q)[k])[1] then
                            Add(arr, ArrowsOfQuiver(Q)[k]);
                        fi;
                    od;
                od;
                arr1 := arr[1];
                for k in [2..Length(arr)] do
                    arr1 := arr1 * arr[k];
                od;
                Add(rho_q, arr1);
            od;
            return rho_q;
        end;

        Composable_Matrix := function(Q,rho)
            local mat, x, str, i, j;
            mat := [];
            for i in [1..NumberOfArrows(Q)] do
                str := [];
                for j in [1..NumberOfArrows(Q)] do
                    if SourceOfPath(ArrowsOfQuiver(Q)[i]) =
                      TargetOfPath(ArrowsOfQuiver(Q)[j]) and not
                      \in(ArrowsOfQuiver(Q)[j] * ArrowsOfQuiver(Q)[i], rho) then
                        Append(str, [1]);
                    else
                        Append(str, [0]);
                    fi;
                od;
                Append(mat, [str]);
            od;
            return mat;
        end;

        IsComposable := function(Q, rho, arrow1, arrow2)
            local i, Mat, pos1, pos2, j;
            Mat := Composable_Matrix(Q, rho);
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

        temp_in := function(Q, i, j)
            return SIntChar(String(IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        temp_out := function(Q, i, j)
            return SIntChar(String(OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i])[j])[1]) - SIntChar('a') + 1;
        end;

        sigma := [];
        eps := [];
        for i in [1..NumberOfArrows(Q)] do
            Append(sigma, [0]);
            Append(eps, [0]);
        od;

        rho := QuiverRho(Q, rho1);

        for i in [1..NumberOfVertices(Q)] do
            temp2 := IncomingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            temp3 := OutgoingArrowsOfVertex(VerticesOfQuiver(Q)[i]);
            if Length(temp2) = 2 then
                eps[temp_in(Q, i, 1)] := 1;
                eps[temp_in(Q, i, 2)] := -1;
                if Length(temp3) = 1 then
                    if IsComposable(Q, rho, temp3[1], temp2[1]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                    elif IsComposable(Q, rho, temp3[1], temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := 1;
                    fi;
                elif Length(temp3) = 2 then
                    if IsComposable(Q, rho, temp3[1], temp2[1]) = 1 or
                      IsComposable(Q, rho, temp3[2], temp2[2]) = 1 then
                        sigma[temp_out(Q, i, 1)] := -1;
                        sigma[temp_out(Q, i, 2)] := 1;
                    else
                        sigma[temp_out(Q, i, 1)] := 1;
                        sigma[temp_out(Q, i, 2)] := -1;
                    fi;
                fi;
            elif Length(temp2) = 1 then
                if Length(temp3) = 1 then
                    eps[temp_in(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 1)] := -1;
                elif Length(temp3) = 2 then
                    sigma[temp_out(Q, i, 1)] := 1;
                    sigma[temp_out(Q, i, 2)] := -1;
                    if IsComposable(Q, rho, temp3[1],
                      temp2[1]) = 1 then
                        eps[temp_in(Q, i, 1)] := -1;
                    else
                        eps[temp_in(Q, i, 1)] := 1;
                    fi;
                fi;
            elif Length(temp2) = 0 then
              sigma[temp_out(Q, i, 1)] := 1;
              sigma[temp_out(Q, i, 2)] := -1;
            fi;
        od;
        return [sigma, eps];
    end;

    MaximumLengthOFDirectString := function(Q, rho)
        local arr, i, max, x, y, array, sigma, eps;
        array := sigma_eps(Q, rho);
        sigma := array[1];
        eps := array[2];
        arr := [];
        for i in [1..NumberOfArrows(Q)] do
            x := [CharInt(i + 96)];
            while true do
                y := direct_left(Q, rho, x, sigma, eps);
                if y = "Cannot Perform The Operation" then
                    Append(arr, [Length(x)]);
                    break;
                else x := y;
                fi;
            od;
        od;
        max := arr[1];
        for i in [2..Length(arr)] do
            if arr[i] > max then max := arr[i];
            fi;
        od;
        return max;
    end;

    NumberOfJoints := function(Q, rho)
        local arr, i, temp, array, sigma, eps;
        arr := [];
        array := sigma_eps(Q, rho);
        sigma := array[1];
        eps := array[2];
        for i in [1..NumberOfArrows(Q)] do
            temp := inverse_right(Q, rho, [CharInt(i+96)], sigma, eps);
            if temp <> "Cannot Perform The Operation" then
                Append(arr, [temp]);
            fi;
        od;
        return Length(arr);
    end;

    CyclicPermutationOfABand := function(Q, rho, band)
        local array, temp, i, j;
        array := [];
        for i in [1..Length(band)] do
            Append(array, [0]);
        od;
        for i in [1..Length(band)] do
            temp := band[1];
            for j in [1..Length(band) - 1] do
                band[j] := band[j + 1];
            od;
            band[j + 1] := temp;
            array[i] := Concatenation(band,"");
        od;
        return array;
    end;

    WeakBridgeQuiver := function(Q, rho, br, q)
        local i, j, k, arr_quiver, arr_quiver1, len, j1, j2, temp, Q1, wb, wb1;
        arr_quiver := [];
        arr_quiver1 := [];
        len := Length(br);
        for i in [1..len] do
            for j in [1..len] do
                #if i <> j then
                j1 := SIntChar(String(TargetOfArrow(Q,br[j][1]))[2]) - SIntChar('0');
                j2 := SIntChar(String(SourceOfArrow(Q,br[i][Length(br[i])]))[2]) - SIntChar('0');
                temp := BandFreeStrings(Q, rho, j1, j2, q);
                for k in [1..Length(temp)] do
                    if IsValidString(Q, rho, Concatenation(br[i], temp[k], br[j]))
                      = "Valid Positive Length String" then
                        Add(arr_quiver1, [j, i, temp[k]]);
                        Add(arr_quiver, [br[j], br[i], temp[k]]);
                    fi;
                od;
                #fi;
            od;
        od;
        wb := arr_quiver;
        wb1 := arr_quiver1;
        Q1 := Quiver(len, arr_quiver1);
        return [Q1, wb, wb1];
    end;

    computeLPSarray := function(pat)
        local lps, i, len, m;
        m := Length(pat);
        lps := [];
        for i in [1..m] do
            Append(lps, [1]);
        od;
        len := 1;
        i := 2;
        while i <= m do
            if pat[i] = pat[len] then
                len := len + 1;
                lps[i] := len;
                i := i + 1;
            else
                if len <> 1 then
                    len := lps[len - 1];
                else
                    lps[i] := 1;
                    i := i + 1;
                fi;
            fi;
        od;
        for i in [1..Length(lps)] do
            lps[i] := lps[i]-1;
        od;
        return lps;
    end;

    KMPSearch := function(pat, txt)
        local M, N, lps, i, j;
        M := Length(pat);
        N := Length(txt);
        lps := computeLPSarray(pat);
        for i in [1..Length(lps)] do
            lps[i] := lps[i] + 1;
        od;
        i := 1;
        j := 1;
        while i <= N do
            if pat[j] = txt[i] then
                j := j + 1;
                i := i + 1;
            fi;
            if j=M + 1 then
                return i - j + 1;
            elif i <= N and pat[j] <> txt[i] then
                if j <> 1 then
                    j := lps[j - 1];
                else
                    i := i + 1;
                fi;
            fi;
        od;
        return 0;
    end;

    SourceOfArrow := function(Q,ch)
        local i, ch1;
        ch1 := ch;
        if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            ch := CharInt(SIntChar(ch) + 32);
        fi;
        if SIntChar(ch)-SIntChar('a') + 1 > NumberOfArrows(Q) then
            return "Invalid String";
        fi;
        for i in [1..NumberOfArrows(Q)] do
            if String(ArrowsOfQuiver(Q)[i])[1] = ch then
                break;
            fi;
        od;
        if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            return TargetOfPath(ArrowsOfQuiver(Q)[i]);
        else
            return SourceOfPath(ArrowsOfQuiver(Q)[i]);
        fi;
    end;

    TargetOfArrow := function(Q, ch)
        local i, ch1;
        ch1 := ch;
        if SIntChar(ch) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            ch := CharInt(SIntChar(ch) + 32);
        fi;
        if SIntChar(ch)-SIntChar('a') + 1 > NumberOfArrows(Q) then
            return "Invalid String";
        fi;
        for i in [1..NumberOfArrows(Q)] do
            if String(ArrowsOfQuiver(Q)[i])[1] = ch then
                break;
            fi;
        od;
        if SIntChar(ch1) < SIntChar('Z') + 1 and SIntChar(ch) > SIntChar('A') - 1 then
            return SourceOfPath(ArrowsOfQuiver(Q)[i]);
        else
            return TargetOfPath(ArrowsOfQuiver(Q)[i]);
        fi;
    end;

    BridgeQuiver := function(Q,rho,br,wb,wb1)
        local lambda, i, j, k, str, temp, array, array1, perm, Q1, num, l;
        temp := [];
        lambda := [];
        for k in [1..Length(wb)] do
            lambda[k] := 1;
        od;
        for i in [1..Length(wb)] do
            for j in [1..Length(wb)] do
                if wb[i][1] = wb[j][2] then
                    str := Concatenation(wb[i][3], wb[j][3]);
                    if IsValidString(Q, rho, str) = "Valid Positive Length String" then
                        if IsBandFreeString(Q, rho, str) = 1 then
                            for k in [1..Length(wb)] do
                                if str = wb[k][3] then
                                    lambda[k] := 0;
                                    break;
                                fi;
                            od;
                        else
                            perm := CyclicPermutationOfABand(Q, rho, wb[i][1]);
                            for k in [1..Length(perm)] do
                                temp := [];
                                num := KMPSearch(perm[k], str);
                                if num <> 0 then
                                    for l in [1..num - 1] do
                                        Add(temp, str[l]);
                                    od;
                                    for l in [num + Length(perm[k])..Length(str)] do
                                        Add(temp, str[l]);
                                    od;
                                    break;
                                fi;
                            od;

                            for k in [1..Length(wb)] do
                                if temp = wb[k][3] then
                                    lambda[k] := 0;
                                    break;
                                fi;
                            od;
                        fi;
                    fi;
                fi;
            od;
        od;

        array := [];
        array1 := [];
        for k in [1..Length(lambda)] do
            if lambda[k] = 1 then
                Add(array, wb[k]);
                Add(array1, wb1[k]);
            fi;
        od;
        Q1 := Quiver(Length(br), array1);
        return Q1;
    end;

    N := NumberOfJoints(Q, rho);
    L := MaximumLengthOFDirectString(Q, rho);
    level := 2 * L * N;
    q := StringsLessThan(Q, rho, level);

    b := [];
    for i in [1..Length(q)] do
        if IsABand(Q, rho, q[i]) = 1 then
            Append(b, [q[i]]);
        fi;
    od;

    tempo := [];
    br := [];
    for i in [1..Length(b)] do
        if \in(b[i], tempo) = false and SIntChar(b[i][1]) >= 97 and
          SIntChar(b[i][Length(b[i])]) < 97 then
            Append(br, [b[i]]);
            Append(tempo, CyclicPermutationOfABand(Q, rho, b[i]));
        fi;
    od;

    count := 0;
    arr_quiver := [];
    len := Length(b);
    for i in [1..len] do
        for j in [1..len] do
            if i <> j then
                temp := Concatenation(b[i], b[j]);
                if IsValidString(Q, rho, temp) = "Valid Positive Length String" then
                    count := count + 1;
                    Append(arr_quiver, [[j, i, Concatenation("a", String(count))]]);
                fi;
            fi;
        od;
    od;

    if IsAcyclicQuiver(Quiver(len, arr_quiver)) = 0 then return 0;
    fi;

    Q1 := WeakBridgeQuiver(Q, rho, br, q);
    Q2 := BridgeQuiver(Q, rho, br, Q1[2], Q1[3]);
    return Q2;
end
);
