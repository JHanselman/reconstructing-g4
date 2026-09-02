/*
 * acb_theta_cli -- command-line front end to FLINT's acb_theta module.
 *
 * Computes theta functions theta_{a,b}(z, tau) for ALL 2^(2g) characteristics
 * (a,b), together with all partial derivatives in z of total order <= ord,
 * using rigorous ball arithmetic (FLINT >= 3.1; FLINT >= 3.2 preferred).
 *
 * This is the external solver behind ThetaFlint() in magma/FlintWrapper.m,
 * in the same spirit as Belyi/Cext/powser_arnoldi.c: Magma writes a small
 * text file, runs this binary, and `eval`s the output.  No Python involved.
 *
 * Usage:
 *     acb_theta_cli INFILE OUTFILE
 *     acb_theta_cli --selftest
 *
 * INFILE (plain text):
 *     line 1:            g prec_bits ord
 *     next 2g lines:     Re z_1, Im z_1, ..., Re z_g, Im z_g
 *     next 2g^2 lines:   Re tau_11, Im tau_11, Re tau_12, ... (row-major)
 * Each number is anything arb_set_str() accepts: a decimal such as
 * "-1.2345e-7", or a ball "1.2345e-7 +/- 3e-40".  The file is tokenised on
 * whitespace (backslash-newline continuations removed first), so Magma's
 * line wrapping of long output is harmless; upper-case E is accepted.
 *
 * OUTFILE is a single Magma expression:
 *     [ [ ComplexField(D) | [re, im], [re, im], ... ], ... ]
 * an outer sequence of length 2^(2g) indexed by the FLINT characteristic
 * number ab = 2^g * a + b (bits of a and b most significant first, i.e.
 * exactly CHIMP's CharacteristicToInteger), each inner sequence holding
 * theta_{a,b}(z,tau) followed by its partial derivatives:
 *     entry 1        : theta_{a,b}(z, tau)
 *     entries 2..g+1 : d/dz_1 theta, ..., d/dz_g theta          (ord >= 1)
 *     then           : higher-order partials in FLINT's tuple order
 *                      (by total order, then reverse-lexicographic; see
 *                      acb_theta_jet_tuples).  Unlike FLINT's raw output,
 *                      these are TRUE partial derivatives, not Taylor
 *                      coefficients: the 1/k! normalisation is undone.
 * D is the number of decimal digits certified for every entry (minimum
 * over all balls, relative to max(1,|x|) so exact zeros do not spoil it).
 *
 * Exit status: 0 on success; 1 on usage/parse error; 2 if the result balls
 * are too wide to be useful (fewer than 5 accurate digits).
 *
 * Build: see Makefile (make / make mac / make server).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_theta.h>

#if __FLINT_RELEASE < 30100
#error "acb_theta_cli needs FLINT >= 3.1 (acb_theta module)"
#endif

#define MAXLINE 65536

static void
die(const char *msg)
{
    fprintf(stderr, "acb_theta_cli: %s\n", msg);
    exit(1);
}

/* ---- input tokenizer -------------------------------------------------------
 * Magma's Sprintf wraps long output at the terminal width: long tokens are
 * split with a trailing backslash, and at whitespace it simply inserts a
 * newline.  So we do not trust line structure at all: the whole file is
 * read, backslash-newline pairs are removed, and it is split on whitespace.
 * A number is then the token MID optionally followed by the tokens "+/-" RAD.
 */
static char *tokbuf = NULL;
static char *tokpos = NULL;

static void
load_input(const char *infile)
{
    FILE *in = fopen(infile, "r");
    long n;
    char *r, *w;
    if (in == NULL) die("cannot open input file");
    fseek(in, 0, SEEK_END);
    n = ftell(in);
    fseek(in, 0, SEEK_SET);
    if (n < 0) die("cannot read input file");
    tokbuf = flint_malloc((size_t) n + 1);
    if (fread(tokbuf, 1, (size_t) n, in) != (size_t) n) die("cannot read input file");
    tokbuf[n] = 0;
    fclose(in);
    /* remove backslash-newline (and backslash-CR-newline) continuations */
    for (r = w = tokbuf; *r; r++)
    {
        if (*r == '\\' && (r[1] == '\n' || (r[1] == '\r' && r[2] == '\n')))
        {
            r += (r[1] == '\r') ? 2 : 1;
            continue;
        }
        *w++ = *r;
    }
    *w = 0;
    tokpos = tokbuf;
}

static char *
next_token(void)
{
    char *start;
    while (*tokpos == ' ' || *tokpos == '\t' || *tokpos == '\n' || *tokpos == '\r')
        tokpos++;
    if (*tokpos == 0)
        return NULL;
    start = tokpos;
    while (*tokpos && *tokpos != ' ' && *tokpos != '\t' && *tokpos != '\n' && *tokpos != '\r')
        tokpos++;
    if (*tokpos)
        *tokpos++ = 0;
    return start;
}

static char *
peek_token(void)
{
    char *save = tokpos, *t;
    while (*tokpos == ' ' || *tokpos == '\t' || *tokpos == '\n' || *tokpos == '\r')
        tokpos++;
    t = tokpos;
    tokpos = save;
    return *t ? t : NULL;
}

/* Read one real number (MID, or MID +/- RAD) into an arb. */
static void
read_arb_tokens(arb_t x, slong prec, const char *what)
{
    char *mid, *rad = NULL, *pk, *q;
    char fixed[MAXLINE];
    char *r;
    mid = next_token();
    if (mid == NULL)
    {
        fprintf(stderr, "acb_theta_cli: unexpected end of input while reading %s\n", what);
        exit(1);
    }
    pk = peek_token();
    if (pk != NULL && strncmp(pk, "+/-", 3) == 0)
    {
        next_token();          /* the "+/-" */
        rad = next_token();
        if (rad == NULL)
        {
            fprintf(stderr, "acb_theta_cli: missing radius after +/- while reading %s\n", what);
            exit(1);
        }
    }
    /* assemble "mid +/- rad" with Magma's upper-case E lowered and a bare
       "." before an exponent (e.g. "0.E-30") padded to "0.0E-30" */
    r = fixed;
    for (q = mid; *q && r < fixed + MAXLINE - 8; q++)
    {
        *r++ = (*q == 'E') ? 'e' : *q;
        if (*q == '.' && (q[1] == 'E' || q[1] == 'e' || q[1] == 0))
            *r++ = '0';
    }
    if (rad != NULL)
    {
        const char *pm = " +/- ";
        for (q = (char *) pm; *q && r < fixed + MAXLINE - 8; q++) *r++ = *q;
        for (q = rad; *q && r < fixed + MAXLINE - 8; q++)
        {
            *r++ = (*q == 'E') ? 'e' : *q;
            if (*q == '.' && (q[1] == 'E' || q[1] == 'e' || q[1] == 0))
                *r++ = '0';
        }
    }
    *r = 0;
    if (arb_set_str(x, fixed, prec) != 0)
    {
        fprintf(stderr, "acb_theta_cli: cannot parse %s from: '%s'\n", what, fixed);
        exit(1);
    }
}

/* Decimal digits certified for a ball, relative to max(1, |x|). */
static slong
acb_digits(const acb_t x)
{
    slong bits;
    if (!acb_is_finite(x))
        return 0;
    bits = acb_rel_one_accuracy_bits(x);
    if (bits > 100000) bits = 100000;   /* exact values report "infinite" */
    if (bits < 0) return 0;
    return (slong) floor(bits * 0.30102999566398120);
}

static void
print_arb_magma(FILE *out, const arb_t x, slong digits)
{
    char *s = arb_get_str(x, digits, ARB_STR_MORE | ARB_STR_NO_RADIUS);
    char *q;
    for (q = s; *q; q++)
        if (*q == 'e') *q = 'E';
    /* arb_get_str prints "0" for an exact zero; give Magma a real literal */
    if (strchr(s, '.') == NULL && strchr(s, 'E') == NULL)
        fprintf(out, "%s.0p%ld", s, (long) digits);
    else
        fprintf(out, "%sp%ld", s, (long) digits);
    flint_free(s);
}

/* Undo the 1/k! Taylor normalisation of FLINT's jets, in place. */
static void
jet_denormalise(acb_ptr th, slong ord, slong g, slong prec)
{
    slong nb = acb_theta_jet_nb(ord, g);
    slong *tups = flint_malloc(nb * g * sizeof(slong));
    slong j, k;
    fmpz_t f, t;
    fmpz_init(f);
    fmpz_init(t);
    acb_theta_jet_tuples(tups, ord, g);
    for (j = 0; j < nb; j++)
    {
        fmpz_one(f);
        for (k = 0; k < g; k++)
        {
            fmpz_fac_ui(t, tups[j * g + k]);
            fmpz_mul(f, f, t);
        }
        if (!fmpz_is_one(f))
            acb_mul_fmpz(th + j, th + j, f, prec);
    }
    fmpz_clear(f);
    fmpz_clear(t);
    flint_free(tups);
}

/* Compute all thetas + jets: th has length 2^(2g) * jet_nb. */
static void
compute(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    slong n2 = 1 << (2 * g);
    slong ab;
#if __FLINT_RELEASE >= 30200
    acb_theta_jet(th, z, 1, tau, ord, 0, 1, 0, prec);
#else
    if (ord == 0)
        acb_theta_all(th, z, tau, 0, prec);
    else
        acb_theta_jet_all(th, z, tau, ord, prec);
#endif
    if (ord >= 2)
        for (ab = 0; ab < n2; ab++)
            jet_denormalise(th + ab * nb, ord, g, prec);
}

static int
run(const char *infile, const char *outfile)
{
    FILE *out;
    slong g, prec, ord, i, j, nb, n2, digits, d;
    acb_ptr z, th;
    acb_mat_t tau;
    long lg, lprec, lord;

    load_input(infile);
    {
        char *t1 = next_token(), *t2 = next_token(), *t3 = next_token();
        if (t1 == NULL || t2 == NULL || t3 == NULL ||
            sscanf(t1, "%ld", &lg) != 1 || sscanf(t2, "%ld", &lprec) != 1 || sscanf(t3, "%ld", &lord) != 1)
            die("header must be: g prec_bits ord");
    }
    g = lg; prec = lprec; ord = lord;
    if (g < 1 || g > 16) die("g out of range (1..16)");
    if (prec < 16) die("prec_bits must be >= 16");
    if (ord < 0 || ord > 64) die("ord out of range (0..64)");

    z = _acb_vec_init(g);
    acb_mat_init(tau, g, g);
    for (i = 0; i < g; i++)
    {
        read_arb_tokens(acb_realref(z + i), prec, "Re z");
        read_arb_tokens(acb_imagref(z + i), prec, "Im z");
    }
    for (i = 0; i < g; i++)
        for (j = 0; j < g; j++)
        {
            read_arb_tokens(acb_realref(acb_mat_entry(tau, i, j)), prec, "Re tau");
            read_arb_tokens(acb_imagref(acb_mat_entry(tau, i, j)), prec, "Im tau");
        }
    flint_free(tokbuf);

    /* force exact symmetry (Magma's tau is only symmetric to working precision) */
    for (i = 0; i < g; i++)
        for (j = i + 1; j < g; j++)
        {
            acb_add(acb_mat_entry(tau, i, j), acb_mat_entry(tau, i, j),
                    acb_mat_entry(tau, j, i), prec);
            acb_mul_2exp_si(acb_mat_entry(tau, i, j), acb_mat_entry(tau, i, j), -1);
            acb_set(acb_mat_entry(tau, j, i), acb_mat_entry(tau, i, j));
        }

    nb = acb_theta_jet_nb(ord, g);
    n2 = 1 << (2 * g);
    th = _acb_vec_init(n2 * nb);

    compute(th, z, tau, ord, prec);

    digits = 100000;
    for (i = 0; i < n2 * nb; i++)
    {
        d = acb_digits(th + i);
        if (d < digits) digits = d;
    }
    if (digits < 5)
    {
        fprintf(stderr, "acb_theta_cli: result has only %ld certified digits; "
                "increase precision or check that Im(tau) is positive definite\n",
                (long) digits);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n2 * nb);
        acb_mat_clear(tau);
        return 2;
    }

    out = fopen(outfile, "w");
    if (out == NULL) die("cannot open output file");
    fprintf(out, "[\n");
    for (i = 0; i < n2; i++)
    {
        fprintf(out, "[ ComplexField(%ld) | ", (long) digits);
        for (j = 0; j < nb; j++)
        {
            if (j) fprintf(out, ", ");
            fprintf(out, "[");
            print_arb_magma(out, acb_realref(th + i * nb + j), digits);
            fprintf(out, ", ");
            print_arb_magma(out, acb_imagref(th + i * nb + j), digits);
            fprintf(out, "]");
        }
        fprintf(out, " ]%s\n", (i + 1 < n2) ? "," : "");
    }
    fprintf(out, "]\n");
    fclose(out);

    _acb_vec_clear(z, g);
    _acb_vec_clear(th, n2 * nb);
    acb_mat_clear(tau);
    return 0;
}

/* ---- self test: genus 2, compare first derivatives to central differences,
 *      and genus 1 theta_00(0, i) = pi^(1/4)/Gamma(3/4).                   */
static int
selftest(void)
{
    slong prec = 300, g = 2, nb, n2, ab, k, fails = 0;
    acb_mat_t tau;
    acb_ptr z, zp, zm, th, thp, thm;
    acb_t h, diff, tmp;
    arb_t tol, absv;

    /* genus 1 sanity value */
    {
        acb_mat_t t1;
        acb_ptr z1, th1;
        arb_t ref, pi, gam;
        acb_mat_init(t1, 1, 1);
        acb_onei(acb_mat_entry(t1, 0, 0));
        z1 = _acb_vec_init(1);
        th1 = _acb_vec_init(4);
        compute(th1, z1, t1, 0, prec);
        arb_init(ref); arb_init(pi); arb_init(gam);
        arb_const_pi(pi, prec);
        arb_root_ui(pi, pi, 4, prec);
        arb_set_d(gam, 0.75);
        arb_gamma(gam, gam, prec);
        arb_div(ref, pi, gam, prec);
        if (!arb_overlaps(acb_realref(th1 + 0), ref) || !arb_contains_zero(acb_imagref(th1 + 0)))
        {
            printf("FAIL: genus-1 theta_00(0,i)\n"); fails++;
        }
        else
            printf("ok: genus-1 theta_00(0,i) = pi^(1/4)/Gamma(3/4)\n");
        arb_clear(ref); arb_clear(pi); arb_clear(gam);
        _acb_vec_clear(z1, 1); _acb_vec_clear(th1, 4); acb_mat_clear(t1);
    }

    acb_mat_init(tau, g, g);
    acb_set_d_d(acb_mat_entry(tau, 0, 0), 0.31, 1.07);
    acb_set_d_d(acb_mat_entry(tau, 0, 1), -0.12, 0.23);
    acb_set(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 0, 1));
    acb_set_d_d(acb_mat_entry(tau, 1, 1), 0.44, 0.91);
    z = _acb_vec_init(g); zp = _acb_vec_init(g); zm = _acb_vec_init(g);
    acb_set_d_d(z + 0, 0.17, 0.05);
    acb_set_d_d(z + 1, -0.21, 0.13);
    nb = acb_theta_jet_nb(2, g);
    n2 = 1 << (2 * g);
    th = _acb_vec_init(n2 * nb);
    thp = _acb_vec_init(n2);
    thm = _acb_vec_init(n2);
    acb_init(h); acb_init(diff); acb_init(tmp);
    arb_init(tol); arb_init(absv);
    arb_set_d(tol, 1e-30);

    compute(th, z, tau, 2, prec);
    acb_set_d(h, 1e-20);
    for (k = 0; k < g; k++)
    {
        _acb_vec_set(zp, z, g); _acb_vec_set(zm, z, g);
        acb_add(zp + k, zp + k, h, prec);
        acb_sub(zm + k, zm + k, h, prec);
        compute(thp, zp, tau, 0, prec);
        compute(thm, zm, tau, 0, prec);
        for (ab = 0; ab < n2; ab++)
        {
            acb_sub(diff, thp + ab, thm + ab, prec);
            acb_div(diff, diff, h, prec);
            acb_mul_2exp_si(diff, diff, -1);           /* central difference */
            acb_sub(tmp, diff, th + ab * nb + 1 + k, prec);
            acb_abs(absv, tmp, prec);
            if (arb_gt(absv, tol))
            {
                printf("FAIL: d/dz_%ld theta_%ld mismatch\n", (long) k, (long) ab);
                fails++;
            }
        }
        /* second derivative d^2/dz_k^2 via (f(z+h) - 2 f(z) + f(z-h)) / h^2 */
        for (ab = 0; ab < n2; ab++)
        {
            slong tup[2] = {0, 0};
            slong idx;
            tup[k] = 2;
            idx = acb_theta_jet_index(tup, g);
            acb_add(diff, thp + ab, thm + ab, prec);
            acb_submul_ui(diff, th + ab * nb, 2, prec);
            acb_div(diff, diff, h, prec);
            acb_div(diff, diff, h, prec);
            acb_sub(tmp, diff, th + ab * nb + idx, prec);
            acb_abs(absv, tmp, prec);
            if (arb_gt(absv, tol))
            {
                printf("FAIL: d^2/dz_%ld^2 theta_%ld mismatch (denormalisation?)\n", (long) k, (long) ab);
                fails++;
            }
        }
    }
    if (!fails) printf("ok: genus-2 first and second derivatives agree with finite differences\n");

    acb_clear(h); acb_clear(diff); acb_clear(tmp); arb_clear(tol); arb_clear(absv);
    _acb_vec_clear(z, g); _acb_vec_clear(zp, g); _acb_vec_clear(zm, g);
    _acb_vec_clear(th, n2 * nb); _acb_vec_clear(thp, n2); _acb_vec_clear(thm, n2);
    acb_mat_clear(tau);
    printf(fails ? "SELFTEST FAILED\n" : "SELFTEST PASSED\n");
    return fails ? 1 : 0;
}

int
main(int argc, char **argv)
{
    int ret;
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        ret = selftest();
    else if (argc == 3)
        ret = run(argv[1], argv[2]);
    else
    {
        fprintf(stderr, "usage: acb_theta_cli INFILE OUTFILE | --selftest\n");
        ret = 1;
    }
    flint_cleanup();
    return ret;
}
