/***
 *  Determining the Sato-Tate group
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

intrinsic SatoTateGroup(EndoStructBase::List, GeoEndoRep::SeqEnum, GalK::List : Shorthand := "") -> MonStgElt
{Given a description EndoStructBase of an endomorphism algebra, a
representation GeoEndoRep of a geometric endomorphism algebra, and a Galois
group GalK, determines the corresponding Sato-Tate group. Via Shorthand, a
description of the geometric endomorphism algebra tensored with RR can be
passed.}

g := #Rows(EndoStructBase[1][1][1]);
if g eq 2 then
    return SatoTateGroupG2QQ(EndoStructBase, GeoEndoRep, GalK : Shorthand := Shorthand);
else
    // TODO: Add other cases when they appear.
    return "undef";
end if;

end intrinsic;


intrinsic SatoTateGroup(EndoStructBase::List, GeoEndoRep::SeqEnum, K::Fld, h::Map : Shorthand := "") -> MonStgElt
{Given a description EndoStructBase of an endomorphism algebra, a
representation GeoEndoRep of a geometric endomorphism algebra, and a field K,
determines the corresponding Sato-Tate group. Via Shorthand, a description of
the geometric endomorphism algebra tensored with RR can be passed.}

/* Apply previous function after finding a corresponding subgroup */
L := BaseRing(EndoStructBase[1][1]);
GalK := SubgroupGeneratorsUpToConjugacy(L, K, h);
return SatoTateGroup(EndoStructBase, GeoEndoRep, GalK : Shorthand := Shorthand);

end intrinsic;


intrinsic SatoTateGroupG2QQ(EndoStructBase::List, GeoEndoRep::SeqEnum, GalK::List : Shorthand := "") -> MonStgElt
{Given a description EndoStructBase of an endomorphism algebra, a
representation GeoEndoRep of a geometric endomorphism algebra, and a Galois
group GalK, determines the corresponding Sato-Tate group. Via Shorthand, a
description of the geometric endomorphism algebra tensored with RR can be
passed. Assumes that the genus equals 2.}

vprint EndoFind: "";
vprint EndoFind: "Calculating Sato-Tate group...";
gensH, Gphi := Explode(GalK);
if #gensH eq 0 then
    Shorthand := SatoTateShorthandG2(EndoStructBase);
    if Shorthand eq "A" then
        return "USp(4)";
    elif Shorthand eq "B" then
        return "G_{3,3}";
    elif Shorthand eq "C" then
        return "G_{1,3}";
    elif Shorthand eq "D" then
        return "F";
    elif Shorthand eq "E" then
        return "E_1";
    elif Shorthand eq "F" then
        return "C_1";
    end if;
end if;
H := sub< Domain(Gphi) | gensH >;
L := BaseRing(GeoEndoRep[1][1]);
GalL := [* sub< Domain(Gphi) |  [ ] >, Gphi *];

/* First case: we are already over the geometric field. It then suffices to
 * calculate the shorthand and we do not have to manipulate the geometric
 * representation. */
if #H eq 1 then
    Shorthand := SatoTateShorthandG2(EndoStructBase);
    if Shorthand eq "A" then
        return "USp(4)";
    elif Shorthand eq "B" then
        return "G_{3,3}";
    elif Shorthand eq "C" then
        return "G_{1,3}";
    elif Shorthand eq "D" then
        return "F";
    elif Shorthand eq "E" then
        return "E_1";
    elif Shorthand eq "F" then
        return "C_1";
    end if;
end if;

/* Determine the Shorthand if it was not passed */
if Shorthand eq "" then
    GeoEndoStructBase := EndomorphismData(GeoEndoRep, GalL);
    Shorthand := SatoTateShorthandG2(GeoEndoStructBase);
end if;
descRR := EndoStructBase[3][3];
K := FixedFieldExtra(L, [ Gphi(gen) : gen in gensH ]);
vprint EndoFind: "Real endomorphism algebra over base field:", descRR;

/* Usually the shorthand and endomorphism structure of the base field determine
 * everything; in the rare cases where they do not we recalculate a bit. */
if Shorthand eq "A" then
    return "USp(4)";

elif Shorthand eq "B" then
    if descRR eq ["RR"] then
        return "N(G_{3,3})";
    elif descRR eq ["RR", "RR"] then
        return "G_{3,3}";
    end if;

elif Shorthand eq "C" then
    if descRR eq ["RR", "RR"] then
        return "N(G_{1,3})";
    elif descRR eq ["RR", "CC"] or descRR eq ["CC", "RR"] then
        return "G_{1,3}";
    end if;

elif Shorthand eq "D" then
    if descRR eq ["RR"] then
        return "F_{ac}";
    elif descRR eq ["RR", "RR"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "F_{ab}";
        elif IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "F_{a,b}";
        end if;
    elif descRR eq ["RR", "CC"] or descRR eq ["CC", "RR"] then
        return "F_a";
    elif descRR eq ["CC", "CC"] then
        return "F";
    end if;

elif Shorthand eq "E" then
    if descRR eq ["RR"] then
        /* Magma being silly with DihedralGroup(2)... */
        if IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "J(E_2)";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "J(E_3)";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "J(E_4)";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "J(E_6)";
        end if;
    elif descRR eq ["RR", "RR"] then
        return "J(E_1)";
    elif descRR eq ["CC"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "E_2";
        elif IsIsomorphic(H, CyclicGroup(3)) then
            return "E_3";
        elif IsIsomorphic(H, CyclicGroup(4)) then
            return "E_4";
        elif IsIsomorphic(H, CyclicGroup(6)) then
            return "E_6";
        end if;
    elif descRR eq ["M_2 (RR)"] then
        return "E_1";
    end if;

elif Shorthand eq "F" then
    if descRR eq ["RR"] then
        if IsIsomorphic(H, DihedralGroup(4)) then
            return "D_{4,1}";

        elif IsIsomorphic(H, DihedralGroup(6)) then
            /* See Fité--Kedlaya--Rotger--Sutherland (4.3) for the next step */
            H_prime := Center(H);
            gensH_prime := Generators(H_prime);
            GalK_prime := [* gensH_prime, Gphi *];
            EndoStruct_prime := EndomorphismData(GeoEndoRep, GalK_prime);
            descRR_prime := EndoStruct_prime[3][3];
            if descRR_prime eq ["M_2 (RR)"] then
                return "D_{6,1}";
            elif descRR_prime eq ["HH"] then
                return "J(D_3)";
            end if;

        elif IsIsomorphic(H, SymmetricGroup(4)) then
            return "O_1";
        elif IsIsomorphic(H, DirectProduct(DirectProduct(CyclicGroup(2), CyclicGroup(2)), CyclicGroup(2))) then
            return "J(D_2)";
        elif IsIsomorphic(H, DirectProduct(DihedralGroup(4), CyclicGroup(2))) then
            return "J(D_4)";
        elif IsIsomorphic(H, DirectProduct(DihedralGroup(6), CyclicGroup(2))) then
            return "J(D_6)";
        elif IsIsomorphic(H, DirectProduct(AlternatingGroup(4), CyclicGroup(2))) then
            return "J(T)";
        elif IsIsomorphic(H, DirectProduct(SymmetricGroup(4), CyclicGroup(2))) then
            return "J(O)";
        end if;

    elif descRR eq ["RR", "RR"] then
        if IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            return "D_{2,1}";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "D_{3,2}";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "D_{4,2}";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "D_{6,2}";
        end if;

    elif descRR eq ["CC"] then
        if IsIsomorphic(H, CyclicGroup(4)) then
            return "C_{4,1}";

        elif IsIsomorphic(H, CyclicGroup(6)) then
            /* Here we take the unique subgroup of order 2 */
            H_prime := Subgroups(H : OrderEqual := 2)[1]`subgroup;
            gensH_prime := Generators(H_prime);
            GalK_prime := [* gensH_prime, Gphi *];
            EndoStruct_prime := EndomorphismData(GeoEndoRep, GalK_prime);
            descRR_prime := EndoStruct_prime[3][3];
            if descRR_prime eq ["M_2 (RR)"] then
                return "C_{6,1}";
            elif descRR_prime eq ["HH"] then
                return "J(C_3)";
            end if;

        elif IsIsomorphic(H, DirectProduct(CyclicGroup(2), CyclicGroup(2))) then
            /* In this case it suffices to check whether the polynomial that
             * defines the center of the geometric endomorphism ring in fact has
             * a root in the ground field */
            struct := EndomorphismDataWithSatoTate(GeoEndoRep, [* [ ], Gphi *]);
            A := struct[2][1];
            A := AlgebraOverCenter(A);
            M := BaseRing(A);
            f := DefiningPolynomial(M);
            if HasRootPari(f, K) then
                return "D_2";
            else
                return "J(C_2)";
            end if;

        elif IsIsomorphic(H, DirectProduct(CyclicGroup(4), CyclicGroup(2))) then
            return "J(C_4)";
        elif IsIsomorphic(H, DirectProduct(CyclicGroup(6), CyclicGroup(2))) then
            return "J(C_6)";
        elif IsIsomorphic(H, DihedralGroup(3)) then
            return "D_3";
        elif IsIsomorphic(H, DihedralGroup(4)) then
            return "D_4";
        elif IsIsomorphic(H, DihedralGroup(6)) then
            return "D_6";
        elif IsIsomorphic(H, AlternatingGroup(4)) then
            return "T";
        elif IsIsomorphic(H, SymmetricGroup(4)) then
            return "O";
        end if;

    elif descRR eq ["CC", "CC"] then
        if IsIsomorphic(H, CyclicGroup(2)) then
            return "C_2";
        elif IsIsomorphic(H, CyclicGroup(3)) then
            return "C_3";
        elif IsIsomorphic(H, CyclicGroup(4)) then
            return "C_4";
        elif IsIsomorphic(H, CyclicGroup(6)) then
            return "C_6";
        end if;

    elif descRR eq ["M_2 (RR)"] then
        return "C_{2,1}";
    elif descRR eq ["HH"] then
        return "J(C_1)";
    elif descRR eq ["M_2 (CC)"] then
        return "C_1";
    end if;
end if;
error "All cases in SatoTateGroupG2QQ fell through";

end intrinsic;


intrinsic SatoTateShorthand(GeoEndoStructBase::List) -> MonStgElt
{Returns the letter describing the neutral connected component of the Sato-Tate
group corresponding to GeoEndoStructBase.}

g := #Rows(GeoEndoStructBase[1][1][1]);
if g eq 2 then
    return SatoTateShorthandG2(GeoEndoStructBase);
else
    /* Add other cases when they appear. */
    return "undef";
end if;

end intrinsic;


intrinsic SatoTateShorthandG2(GeoEndoStructBase::List) -> MonStgElt
{Returns the letter describing the neutral connected component of the Sato-Tate
group corresponding to GeoEndoStructBase. Assumes that the genus equals 2.}

descRR := GeoEndoStructBase[3][3];
case descRR:
    when ["RR"]:       return "A";
    when ["RR", "RR"]: return "B";
    when ["RR", "CC"]: return "C";
    when ["CC", "RR"]: return "C";
    when ["CC", "CC"]: return "D";
    when ["M_2 (RR)"]: return "E";
    when ["M_2 (CC)"]: return "F";
    else: error "Shorthand algorithm obtains contradiction with classification";
end case;

end intrinsic;
