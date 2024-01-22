/*
 * Functions to compute the last minors of a matrix developed on its final row.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2023 Guilhem Niot
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Guilhem Niot <guilhem@gniot.fr>
 */

#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "perm.h"
#include <stddef.h>
#include <stdint.h>

/** CODE OF FLINT THAT WAS NOT COMMITED YET.
 * Permute rows of a matrix `mat` according to `perm_act`, and propagate the
 * action on `perm_store`.
 * That is, performs for each appropriate index `i`, the operations
 * `perm_store[i] <- perm_store[perm_act[i]]`
 * `rows[i] <- rows[perm_act[i]]` */
static void local_nmod_mat_permute_rows(nmod_mat_t mat, const slong *perm_act,
                                        slong *perm_store)
{
    slong i;
    mp_limb_t **mat_tmp = flint_malloc(mat->r * sizeof(mp_limb_t *));

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (i = 0; i < mat->r; i++)
        mat_tmp[i] = mat->rows[perm_act[i]];
    for (i = 0; i < mat->r; i++)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}

/**
 * Code inspired by the double_det function of Flint and C. Pernet, W. Stein,
 * Fast computation of Hermite normal forms of random integer matrices, J.
 * Number Theory (2010), doi:10.1016/j.jnt.2010.01.017.
 *
 * Takes a n x (n+1) matrix B. Compute the nb_minors last minors of the matrix
 * developped on the (n+1)-th row. Assumes that nb_minors is even.
 */
void compute_minors(fmpz *minors, slong nb_minors, const fmpz_mat_t B)
{
    slong i, j, l, n, l0;
    slong *P;
    mp_limb_t p;
    mp_ptr uimod, vimod;
    fmpz_t bound, prod, s1, s2, t;
    fmpz *ui, *vi;
    fmpz_mat_t dt, Bt;
    fmpq_t tmpq;
    fmpq_mat_t x;
    nmod_mat_t Btmod, dtmod;

    n = B->r;

    fmpz_mat_init(dt, n, 1);
    fmpz_mat_init(Bt, n, n);
    fmpq_mat_init(x, n, 1);

    for (i = 0; i < n; i++)
    {
        // Copy all B but the last column
        for (j = 0; j < n; j++)
        {
            fmpz_set(fmpz_mat_entry(Bt, i, j), fmpz_mat_entry(B, i, j));
        }

        // Last column of B
        fmpz_set(fmpz_mat_entry(dt, i, 0), fmpz_mat_entry(B, i, n));
    }

    /* solve B_{1,..,n}^Tx = B_{n+1}^T */
    fmpq_mat_solve_fmpz_mat(x, Bt, dt);

    fmpz_init(bound);
    fmpz_init(prod);
    fmpz_init(t);
    fmpz_init(s1);
    fmpz_init(s2);
    ui = _fmpz_vec_init(nb_minors); // Will contain dividers of each minor
    vi = _fmpz_vec_init(nb_minors);
    uimod = _nmod_vec_init(nb_minors);
    vimod = _nmod_vec_init(nb_minors);

    /* compute lcm of denominators of vectors x and ys */
    fmpq_init(tmpq);
    for (l = 0; l < nb_minors; l++)
    {
        fmpz_one(&ui[l]);
    }
    for (i = 0; i < n; i++)
    {
        fmpz_lcm(&ui[0], &ui[0], fmpq_mat_entry_den(x, i, 0));
        for (l = 1; l < nb_minors; l++)
        {
            if (fmpq_is_zero(fmpq_mat_entry(x, n - l, 0)))
                continue;

            fmpq_div(tmpq, fmpq_mat_entry(x, i, 0), fmpq_mat_entry(x, n - l, 0));
            fmpz_lcm(&ui[l], &ui[l], fmpq_denref(tmpq));
        }
    }
    for (l = 1; l < nb_minors; l++)
    {
        if (fmpq_is_zero(fmpq_mat_entry(x, n - l, 0)))
            continue;

        fmpq_inv(tmpq, fmpq_mat_entry(x, n - l, 0));
        fmpz_lcm(&ui[l], &ui[l], fmpq_denref(tmpq));
    }
    fmpq_clear(tmpq);

    /* compute Hadamard bounds */
    /* Slight overestimate as we could remove one of the columns to bound each
        * minor */
    fmpz_one(s2);
    for (j = 0; j < n + 1; j++)
    {
        fmpz_zero(s1);
        for (i = 0; i < n; i++)
            fmpz_addmul(s1, fmpz_mat_entry(Bt, i, j), fmpz_mat_entry(Bt, i, j));
        fmpz_sqrtrem(s1, t, s1);
        if (!fmpz_is_zero(t))
            fmpz_add_ui(s1, s1, UWORD(1));
        fmpz_mul(s2, s2, s1);
    }

    for (l = 0; l < nb_minors; l++)
    {
        fmpz_cdiv_q(t, s2, &ui[l]);
        if (fmpz_cmp(bound, t) < 0)
        {
            fmpz_set(bound, t);
        }
    }
    // We need a bigger bound to recover the sign
    fmpz_mul_ui(bound, bound, UWORD(2));

    fmpz_one(prod);
    P = _perm_init(n);
    nmod_mat_init(Btmod, n, n, 2);
    nmod_mat_init(dtmod, n, 1, 2);
    p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
    /* compute determinants divided by u1 and u2 */
    while (fmpz_cmp(prod, bound) <= 0)
    {
        p = n_nextprime(p, 0);

        // If p divides one the divisors, skip it
        int cont = 1;
        for (l = 0; l < nb_minors; l++)
        {
            uimod[l] = fmpz_fdiv_ui(&ui[l], p);
            if (!uimod[l])
            {
                cont = 0;
                break;
            }
        }
        if (!cont) continue;

        _nmod_mat_set_mod(Btmod, p);
        _nmod_mat_set_mod(dtmod, p);

        for (l0 = 0; l0 < nb_minors; l0++) {

            // We start by computing the LU decomposition of the matrix without its n-l0 column
            // Copy the entire matrix except last column
            for (i = 0; i < n; i++) {
                for (j = 0; j < n-l0; j++)
                    nmod_mat_entry(Btmod, i, j) = fmpz_fdiv_ui(fmpz_mat_entry(B, i, j), p);
                for (j = n-l0+1; j <= n; j++)
                    nmod_mat_entry(Btmod, i, j) = fmpz_fdiv_ui(fmpz_mat_entry(B, i, j), p);
            }

            nmod_mat_lu(P, Btmod, 0);

            // We can already compute minors[l0]
            vimod[l0] = UWORD(1);
            for (i = 0; i < n; i++)
                vimod[l0] = n_mulmod2_preinv(vimod[l0], nmod_mat_entry(Btmod, i, i), p,
                                                Btmod->mod.ninv);
            if (_perm_parity(P, n) == 1)
                vimod[l0] = nmod_neg(vimod[l0], Btmod->mod);

            // In case A is rank deficient mod p, we cannot directly compute all the minors
            if (vimod[l0] == 0) {
                continue;
            }

            break;
        }

        // Then, compute the other minors
        // Compute the coordinates of the last vector in the basis
        for (i = 0; i < n; i++)
        {
            nmod_mat_entry(dtmod, i, 0) =
                fmpz_fdiv_ui(fmpz_mat_entry(B, i, n - l0), p); // Column n-l
        }
        local_nmod_mat_permute_rows(dtmod, P, NULL);
        nmod_mat_solve_tril(dtmod, Btmod, dtmod, 1);
        nmod_mat_solve_triu(dtmod, Btmod, dtmod, 0);

        for (int l = l0 + 1; l < nb_minors; l++)
        {
            vimod[l] = n_mulmod2_preinv(vimod[l0], nmod_mat_entry(dtmod, n - l, 0), p,
                                        Btmod->mod.ninv);
            if ((l - l0) % 2 == 0)
            {
                vimod[l] = nmod_neg(vimod[l], Btmod->mod);
            }
        }

        // Update CRT
        for (int l = 0; l < nb_minors; l++) {
            vimod[l] = n_mulmod2_preinv(vimod[l], n_invmod(uimod[l], p),
                                            p, Btmod->mod.ninv);
            fmpz_CRT_ui(&vi[l], &vi[l], prod, vimod[l], p, 1);
        }

        fmpz_mul_ui(prod, prod, p);
    }

    for (l = 0; l < nb_minors; l++)
    {
        fmpz_mul(&minors[l], &ui[l], &vi[l]);
    }

    fmpz_clear(bound);
    fmpz_clear(prod);
    fmpz_clear(s1);
    fmpz_clear(s2);
    fmpz_clear(t);
    _fmpz_vec_clear(ui, nb_minors);
    _fmpz_vec_clear(vi, nb_minors);
    _nmod_vec_clear(uimod);
    _nmod_vec_clear(vimod);
    _perm_clear(P);
    nmod_mat_clear(Btmod);
    nmod_mat_clear(dtmod);

    fmpz_mat_clear(dt);
    fmpz_mat_clear(Bt);
    fmpq_mat_clear(x);
}

int main()
{
    int n = 1000;

    flint_rand_t rstate;
    flint_randinit(rstate);

    fmpz_mat_t B;
    fmpz_mat_init(B, n, n+1);

    fmpz_mat_randbits(B, rstate, 7);


    int nbminors = 100;
    fmpz* minors = _fmpz_vec_init(nbminors);

    printf("execute new opti\n\n");
    compute_minors(minors, nbminors, B);

    // for (int i = 0; i < nbminors; i++) {
    //     fmpz_print(&minors[i]);
    // }
}