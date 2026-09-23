// Comparing the token multisets that RealRepresentationString produces, and
// escalating B against a known truth. Separate from UpperBounds.m because this
// is comparison vocabulary, not part of the CMSV bound itself.

// The simple real algebra M_k(D) behind a token, as <k, dim_R D>. The whole
// vocabulary is what RealRepresentationString emits: RR, CC, HH and M_k of
// those, k unbounded (a rung with few primes can report M_4(RR)).
function SimpleFactorData(tok)
    dims := [1, 2, 4];
    base := Index(["RR", "CC", "HH"], tok);
    if base ne 0 then
        return true, 1, dims[base];
    end if;
    j := Index(tok, "(");
    if #tok lt 6 or tok[1 .. 2] ne "M_" or j lt 4 or tok[#tok] ne ")" then
        return false, 0, 0;
    end if;
    digits := tok[3 .. j - 1];
    base := Index(["RR", "CC", "HH"], Substring(tok, j + 1, #tok - j - 1));
    if base eq 0 or not forall{c : c in Eltseq(digits) | c in ["0" .. "9"]} then
        return false, 0, 0;
    end if;
    k := StringToInteger(digits);
    if k lt 1 then
        return false, 0, 0;
    end if;
    return true, k, dims[base];
end function;

// RealRepresentationString cannot always separate type II from type III
// (utils.m:263). Containment is reported when it holds for either resolution.
function TokenCandidates(tok)
    if tok cmpeq "M_2(RR) or HH" then
        return ["M_2(RR)", "HH"];
    end if;
    return [tok];
end function;

function IsKnownToken(tok)
    for c in TokenCandidates(tok) do
        if not SimpleFactorData(c) then
            return false;
        end if;
    end for;
    return true;
end function;

// Every way of replacing each ambiguous token by a definite one.
function Resolutions(items)
    out := [items];
    for i in [1 .. #items] do
        cands := TokenCandidates(items[i]);
        if #cands eq 1 then
            continue;
        end if;
        next := [];
        for row in out do
            for c in cands do
                r := row;
                r[i] := c;
                Append(~next, r);
            end for;
        end for;
        out := next;
    end for;
    return out;
end function;

// Total dim_R; dim_R M_k(D) = k^2 dim_R D. Both candidates behind the
// ambiguous token have dimension 4, so the resolution chosen does not matter.
function TotalRealDimension(items)
    total := 0;
    for tok in Resolutions(items)[1] do
        _, k, d := SimpleFactorData(tok);
        total +:= k^2 * d;
    end for;
    return total;
end function;

// Supports, as index sets, of the solutions in non-negative integers of
// sum_i c_i units[i] = target, excluding the all-zero c.
function SupportMasks(units, target)
    walk := function(i, rest, mask)
        if i gt #units then
            if rest eq 0 and #mask gt 0 then
                return {mask};
            end if;
            return {};
        end if;
        out := {};
        for c in [0 .. rest div units[i]] do
            out join:= $$(i + 1, rest - c * units[i],
                          (c gt 0 select (mask join {i}) else mask));
        end for;
        return out;
    end function;
    return walk(1, target, {});
end function;

// Per factor M_n(E) of bound: sum_i c_i k_i max(dim D_i, dim E) = n dim E, with
// no factor of bound left unhit (unitality) and none of target unused
// (injectivity). The unit is the real dimension of the simple M_k(D) (x) E
// module, so the c_i are the multiplicities of the isotypic pieces.
function EmbedsResolved(target, bound)
    shape := [];
    for tok in target do
        _, k, d := SimpleFactorData(tok);
        Append(~shape, <k, d>);
    end for;
    reach := {{}};
    for tok in bound do
        _, n, e := SimpleFactorData(tok);
        masks := SupportMasks([s[1] * Max(s[2], e) : s in shape], n * e);
        if #masks eq 0 then
            return false;
        end if;
        reach := {r join m : r in reach, m in masks};
    end for;
    return {1 .. #target} in reach;
end function;

intrinsic RealRepresentationEmbeds(target::SeqEnum, bound::SeqEnum) -> BoolElt
{True if prod(target) admits a unital injective R-algebra homomorphism into
 prod(bound), for some resolution of the ambiguous token. Both arguments are
 multisets of the tokens produced by RealRepresentationString}
    for tok in target cat bound do
        require Type(tok) eq MonStgElt and IsKnownToken(tok):
            Sprintf("unknown algebra token %o", tok);
    end for;
    for t in Resolutions(target) do
        for b in Resolutions(bound) do
            if EmbedsResolved(t, b) then
                return true;
            end if;
        end for;
    end for;
    return false;
end intrinsic;

intrinsic RealRepresentationBound(C::Crv, target::SeqEnum : Bmax := 1024)
    -> SeqEnum, RngIntElt, MonStgElt
{Searches increasing prime bounds, capped at Bmax (default 1024), until the
 flattened bound for C equals the known truth target. Returns <bound, B,
 status>: "sharp" on a match, "exhausted" with the tightest bound containing
 target, "unsound" only if the last rung with a bound missed it}
    sorted := Sort(target);
    best := [];
    bestB := 0;
    bestdim := -1;
    bad := [];
    badB := 0;
    lastbad := false;
    frobs := [];
    prev := 2;
    B := 8;
    while B le Bmax do
        rung := B;
        B *:= 2;
        // Refolds eta/centre per rung, not accumulating as CLV: a known cost.
        frobs cat:= [pair[2] : pair in LPolynomials(C, prev, rung)];
        prev := rung;
        if #frobs eq 0 then
            continue;
        end if;
        raw := RealRepresentationBound(frobs);
        // An empty bound means no prime below B was usable, not a failure.
        if #raw eq 0 then
            continue;
        end if;
        bound := Sort(&cat raw);
        if #bound eq 0 then
            continue;
        end if;
        if not RealRepresentationEmbeds(sorted, bound) then
            // The factor and centre bounds hold only if eta and t are right,
            // which a handful of primes does not give, so a low rung may
            // legitimately miss the truth. Keep climbing; only the last rung
            // with a bound decides, and never record it as the tightest one.
            bad := bound;
            badB := rung;
            lastbad := true;
            continue;
        end if;
        lastbad := false;
        if bound eq sorted then
            return bound, rung, "sharp";
        end if;
        // Ties go to the larger B, hence le.
        dim := TotalRealDimension(bound);
        if bestdim lt 0 or dim le bestdim then
            best := bound;
            bestB := rung;
            bestdim := dim;
        end if;
    end while;
    if lastbad then
        return bad, badB, "unsound";
    end if;
    return best, bestB, "exhausted";
end intrinsic;
