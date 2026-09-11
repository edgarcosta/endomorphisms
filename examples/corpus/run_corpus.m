// Magma driver: run RealRepresentationBound over one chunk of a prepared corpus.
//
// Everything comes from the environment, because Magma's command-line argument
// handling is awkward:
//
//   CORPUS_INPUT   prepared TSV chunk (id, kind, data, expected[, meta])
//   CORPUS_OUTPUT  results TSV, appended to, one line flushed per curve
//   CORPUS_B       prime bound handed to RealRepresentationBound
//   ENDO_SPEC      path to <repo>/endomorphisms/magma/spec
//   POLRED_SPEC    path to CHIMP/MagmaPolred/spec (optional if already attached)
//
// Output columns: id, kind, B, expected, got, cputime, status.
// 'got' is the comma-joined flattened bound; on status 'error' it carries the
// error message instead. Lines are written with PrintFile, which opens/closes
// per call and therefore flushes: a chunk killed mid-run keeps everything it
// had already produced, which is what makes run_corpus.sh resumable. Magma
// block-buffers redirected stdout, so printing would not have that property.
//
// No asserts and no SetQuitOnError: a single bad curve must not take out the
// chunk. The one hard failure is the gp guard below.
//
// Run directly with:
//   CORPUS_INPUT=... CORPUS_OUTPUT=... CORPUS_B=50 magma -b run_corpus.m
// but normally you want run_corpus.sh, which shards and parallelizes.

// ---------------------------------------------------------------- utilities

// Magma's Split drops empty fields, which would silently shift columns.
function SplitOn(s, sep)
    out := [];
    rest := s;
    while true do
        i := Index(rest, sep);
        if i eq 0 then
            Append(~out, rest);
            return out;
        end if;
        Append(~out, rest[1 .. i - 1]);
        rest := rest[i + 1 .. #rest];
    end while;
end function;

// Keep tabs and newlines out of a field (error messages are the usual source).
function Flatten1Line(s)
    if #s eq 0 then
        return s;
    end if;
    return &cat [(c eq "\t" or c eq "\n" or c eq "\r") select " " else c
                 : c in Eltseq(s)];
end function;

function ErrorText(e)
    txt := "unprintable error";
    try
        txt := Sprint(e`Object);
    catch inner
        txt := "unprintable error";
    end try;
    return Flatten1Line(txt);
end function;

function TryAttachSpec(path)
    if path eq "" then
        return false;
    end if;
    ok := true;
    try
        AttachSpec(path);
    catch e
        ok := false;
    end try;
    return ok;
end function;

// -------------------------------------------------------------- environment

input_path := GetEnv("CORPUS_INPUT");
output_path := GetEnv("CORPUS_OUTPUT");
B_text := GetEnv("CORPUS_B");
endo_spec := GetEnv("ENDO_SPEC");
polred_spec := GetEnv("POLRED_SPEC");

if input_path eq "" or output_path eq "" then
    print "run_corpus: CORPUS_INPUT and CORPUS_OUTPUT must be set.";
    exit 2;
end if;

if B_text eq "" then
    B_text := "50";
end if;
B := StringToInteger(B_text);
if B lt 3 then
    printf "run_corpus: CORPUS_B = %o is too small.\n", B;
    exit 2;
end if;

// ENDO_SPEC should be set by run_corpus.sh, which derives it from its own
// location. The fallbacks only cover running this file by hand from the
// repository root or from examples/corpus.
if not TryAttachSpec(endo_spec) then
    attached := false;
    for candidate in ["../../endomorphisms/magma/spec",
                      "endomorphisms/magma/spec",
                      "../endomorphisms/magma/spec"] do
        if TryAttachSpec(candidate) then
            attached := true;
            break;
        end if;
    end for;
    if not attached then
        printf "run_corpus: cannot attach the endomorphisms spec (ENDO_SPEC = %o).\n",
               endo_spec;
        exit 2;
    end if;
end if;

// POLRED_SPEC is optional: Polredabs may already come from MAGMA_USER_SPEC.
// Either way the guard below decides whether we are allowed to proceed.
if polred_spec ne "" then
    if not TryAttachSpec(polred_spec) then
        printf "run_corpus: cannot attach POLRED_SPEC = %o.\n", polred_spec;
        exit 2;
    end if;
end if;

SetVerbose("EndoFind", 0);

// The g3 models are literal strings in x, y, z, so P2 must own those
// names and the univariate ring must not.
R<t> := PolynomialRing(Rationals());
P2<x, y, z> := ProjectiveSpace(Rationals(), 2);
P2CR := CoordinateRing(P2);

// ------------------------------------------------------------------ gp guard
//
// FieldIntersectionMatrix compares centers by set intersection on polredabs'd
// defining polynomials, so polredabs is the canonical key that makes field
// equality testable at all. Polredabs shells out to PARI/gp; when gp is
// missing, Polred catches the failed Pipe, prints a warning and returns a
// NON-canonical fallback. Results are then silently wrong: the meet loses
// elements, centers collapse to Q, and e.g. ["RR","CC"] degrades to
// ["RR","RR"]. So test by round-trip rather than by looking for gp on PATH,
// and refuse to run at all if it does not hold.

gp_ok := false;
try
    gp_ok := (Polredabs(t^2 - 5) eq t^2 - t - 1);
catch e
    gp_ok := false;
end try;

if not gp_ok then
    print "";
    print "run_corpus: FATAL - PARI/gp is not usable from Magma.";
    print "  Polredabs(x^2 - 5) must return x^2 - x - 1.";
    print "  Without a working gp, Polred falls back to a non-canonical";
    print "  polynomial, center detection collapses to Q and the computed";
    print "  bounds are silently wrong (e.g. [RR, CC] degrades to [RR, RR]).";
    print "  Put a working gp on PATH and re-run; Sage ships one, e.g.";
    print "    PATH=\"$(dirname \"$(sage -sh -c 'command -v gp')\"):$PATH\"";
    exit 2;
end if;

// ------------------------------------------------------------------ the work

function BuildCurve(kind, data)
    if kind eq "g2" then
        coeffs := eval("return " cat data cat ";");
        return HyperellipticCurve(R ! coeffs[1], R ! coeffs[2]);
    elif kind eq "g3" then
        return Curve(P2, P2CR ! eval("return " cat data cat ";"));
    end if;
    error "unknown kind " cat kind;
end function;

lines := Split(Read(input_path), "\n");
ndone := 0;

for line in lines do
    if #line eq 0 then
        continue;
    end if;
    fields := SplitOn(line, "\t");
    if #fields lt 4 then
        continue;
    end if;
    id := fields[1];
    kind := fields[2];
    data := fields[3];
    expected := fields[4];
    if id eq "id" and kind eq "kind" then  // header emitted by prepare_input.py
        continue;
    end if;

    t0 := Cputime();
    status := "ok";
    got := "";
    try
        C := BuildCurve(kind, data);
        bound := RealRepresentationBound(C, B);
        if #bound eq 0 then
            // The upper bound could not be established; report it as an empty
            // multiset rather than an error, and let report.py classify it.
            flat := [Strings() | ];
        else
            flat := &cat bound;
        end if;
        got := Join([Strings() | Flatten1Line(s) : s in flat], ",");
    catch e
        status := "error";
        got := ErrorText(e);
    end try;
    elapsed := Cputime(t0);

    PrintFile(output_path, Sprintf("%o\t%o\t%o\t%o\t%o\t%.3o\t%o",
                                   id, kind, B, expected, got, elapsed, status));
    ndone +:= 1;
end for;

printf "run_corpus: %o curves at B = %o -> %o\n", ndone, B, output_path;
exit 0;
