intrinsic computeLPolys2(N::.) -> .
{docstring}
    _<x1,x2,x3> := PolynomialRing(Integers(), 3);
    Picard := x1^3*x3 + x2^2*x3^2 + 3*x2^4 + 7*x2*x3^3;

    C := Curve(ProjectiveSpace(Integers(),2), Picard);

    R<x> := PolynomialRing(Integers());

    LPolys := [R!0 : k in [0..230]];

    for p in [1..N] do
        if IsPrime(p) then
            Cp := BaseChange(C, FiniteField(p));
            if IsNonsingular(Cp) then
                print R!Reverse(Coefficients(LPolynomial(Cp)));
                LPolys[p] := R!Reverse(Coefficients(LPolynomial(Cp)));
            end if;
        end if;
    end for;
    return LPolys;
end intrinsic;


intrinsic computeLPolys3(N::.) -> .
{docstring}
    _<x1,x2,x3> := PolynomialRing(Integers(), 3);
    Picard := x1^3*x3 + x2^2*x3^2 + 3*x2^4 + 2*x2*x3^3;

    C := Curve(ProjectiveSpace(Integers(),2), Picard);

    R<x> := PolynomialRing(Integers());

    LPolys := [R!0 : k in [0..230]];

    for p in [1..N] do
        if IsPrime(p) then
            Cp := BaseChange(C, FiniteField(p));
            if IsNonsingular(Cp) then
                print R!Reverse(Coefficients(LPolynomial(Cp)));
                LPolys[p] := R!Reverse(Coefficients(LPolynomial(Cp)));
            end if;
        end if;
    end for;
    return LPolys;
end intrinsic;


intrinsic computeLPolys1(N::.) -> .
{docstring}
    _<x1,x2,x3> := PolynomialRing(Integers(), 3);
    Q11 :=
    -4169*x1^4 - 956*x1^3*x2 + 7440*x1^3*x3 + 55770*x1^2*x2^2 +
    43486*x1^2*x2*x3 + 42796*x1^2*x3^2 - 38748*x1*x2^3 - 30668*x1*x2^2*x3 +
    79352*x1*x2*x3^2 - 162240*x1*x3^3 + 6095*x2^4 + 19886*x2^3*x3 -
    89869*x2^2*x3^2 - 1079572*x2*x3^3 - 6084*x3^4;

    C := Curve(ProjectiveSpace(Integers(),2), Q11);

    R<x> := PolynomialRing(Integers());

    LPolys := [R!0 : k in [0..230]];

    for p in [1..N] do
        if IsPrime(p) then
            Cp := BaseChange(C, FiniteField(p));
            if IsNonsingular(Cp) then
                print R!Reverse(Coefficients(LPolynomial(Cp)));
                LPolys[p] := R!Reverse(Coefficients(LPolynomial(Cp)));
            end if;
        end if;
    end for;
    return LPolys;
end intrinsic;
