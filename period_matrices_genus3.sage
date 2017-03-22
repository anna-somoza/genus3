"""
Period_matrices_of_genus_3_CM_curves -- Sage package for computing 
period matrices of genus 3 CM curves

#*****************************************************************************
# Copyright (C) 2016 Marco Streng <marco.streng@gmail.com> and 
# Pinar Kilicer <pinarkilicer@gmail.com>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*****************************************************************************

"""

"""
 
load_attach_mode(True,True)
load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

prec_bits = 200000
prec_digits = floor((log(2)/log(10))*prec_bits)

This code is a modification of "recip" (which is GPL)
for sextic CM fields

"""

def period_matrix(bas, embs):
    g = ZZ(bas.degree() / 2)
    M = Matrix([[phi(b) for b in bas] for phi in embs])
    Omega1 = M[:,:g]
    Omega2 = M[:,g:]
    return Omega2.inverse()*Omega1


def split(bas):
    g = ZZ(bas.degree()/2)
    return vector(bas[:g]), vector(bas[g:])

def join(bas1,bas2):
    return vector(list(bas1)+list(bas2))

def B_to_M(B):
    g = B.nrows()
    assert B.ncols() == g
    return ABCD_to_M(identity_matrix(g), B, zero_matrix(g, g), identity_matrix(g))

def ABCD_to_M(A,B,C,D):
    g = A.ncols()
    for E in [A,B,C,D]:
        assert E.nrows() == E.ncols() == g
    return Matrix([join(A[i],B[i]) for i in range(A.nrows())] +
                  [join(C[i],D[i]) for i in range(A.nrows())])


def reduce_Siegel_epsilon(bas, embs, verbose=False):
    """
    Given a symplectic basis and the embeddings from a CM-type,
    reduces to the domain F_S^epsilon of the notes:
    
      * Re(tau) has coeffs in [-1/2,1/2]
      * Im(tau) is LLL-reduced
      * Im(tau[0,0]) >= sqrt(3/4) - 0.01

    returns a new basis and a basis transformation
    """
    input_bas = bas
    g = ZZ(bas.degree()/2)
    M = identity_matrix(2*g)
    
    while True:
        # the imaginary part:
        bas, N = reduce_imag_part_LLL(bas, embs)
        M = N * M
        assert act_by_M_on_basis(M, input_bas) == bas
        if verbose:
            print "imag reduction:"
            print N
        # the real part:
        bas, N = reduce_real_part(bas, embs)
        M = N * M
        assert act_by_M_on_basis(M, input_bas) == bas
        if verbose:
            print "real reduction:"
            print N
        # the upper left imaginary entry
        tau = period_matrix(bas, embs)
        if abs(tau[0,0]) < 0.99:
            N = M_for_reduce_using_z11(g)
            M = N * M
            bas = act_by_M_on_basis(N, bas)
            assert act_by_M_on_basis(M, input_bas) == bas
            if verbose:
                print "y11 reduction:"
                print N

        else:
            assert act_by_M_on_basis(M, input_bas) == bas

            return bas, M



def reduce_real_part(bas, embs):
    """
    Given a symplectic basis and the embeddings from a CM-type,
    change the basis such that
    
      * Re(tau) has coeffs in [-1/2,1/2]
    
    returns a new basis and a basis transformation
    """
    tau = period_matrix(bas, embs)
    g = tau.ncols()
    B = Matrix(ZZ, [[0 for i in range(g)] for j in range(g)])
    for i in range(g):
        for j in range(g):
            c = - ZZ(tau[i,j].real().round())
            B[i,j] = c
            B[j,i] = c
            # TODO: maybe test whether tau is "sufficiently symmetric"
            #(tau should be symmetric, but rounding may make it non-symmetric)
    M = B_to_M(B)
    bas1, bas2 = split(bas)
    bas1 = bas1 + bas2*B
    bas_new = join(bas1, bas2)
    assert bas_new == bas * M.transpose()
    return bas_new, M



def reduce_imag_part_LLL(bas, embs, prec=None):
    """
    Given a symplectic basis and the embeddings from a CM-type,
    change the basis such that
    
      * Im(tau) is LLL-reduced

    returns a new basis and a basis transformation
    """
    g = ZZ(bas.degree() / 2)
    if prec is None:
        prec = embs[0].codomain().precision()
    tau = period_matrix(bas, embs)
    Y = Matrix([[ZZ((2^prec*t.imag()).round()) for t in u] for u in tau])
    U = Y.LLL_gram()
    # U^T * Y * U is reduced
    # So want: M = (U^T  0 )
    #              ( 0  U^-1 )
    M = ABCD_to_M(U.transpose(), zero_matrix(g, g), zero_matrix(g, g), U.inverse())
    return bas*M.transpose(), M



def M_for_reduce_using_z11(g):
    M = zero_matrix(2*g,2*g)
    M[0,g] = -1
    M[g,0] = 1
    for i in range(1,g):
        M[i,i] = 1
        M[g+i,g+i] = 1
    return M


def act_by_M_on_basis(M, bas):
    return bas*M.transpose()


def find_xi_g3(K, B):
    """
    OUTPUT: True, xi or False, None
    Returns a totally positive imaginary xi generating B if it exists.
    
    Assumes \OO*_K = W_K\OO*_F and assumes that \OO*_F takes all 8 signs.
    """
    assert B.is_principal()
    xi = B.gens_reduced()[0]
    u = K.unit_group()
    mu = K(u.torsion_generator())
    k = mu.multiplicative_order()
    for l in range(k/2):
        # test whether mu^l is totally imaginary:
        K0b, mb = K.subfield((mu**l*xi)**2, 'b')
        if K0b.degree() == len(K0b.embeddings(AA)) and (-K0b.gen()).is_totally_positive():
            xi = xi*mu**l
            return True, xi

    return False, None



def principally_polarized_ideal_classes(K):
    """
    Returns a set of principally polarized idealss (A, xi)
    that covers all principally polarizable ideal classes [A] exactly once.

        Note that Lemma 3.2.2 in Pinar's thesis says \OO*_K = W_K\OO*_F.
        So it is enough to multiply xi with a generator mu of W_K to get a
        totally imaginary generator.
    """
    lst = []
    c = K.complex_conjugation()
    D = K.different()
    for A in K.class_group_iter():
        # totally_real_cubic_subfield_and_embeddings(K)[0]
        #A = a.ideal()
        ideal_xi = (A*c(A)*D)**-1
        if ideal_xi.is_principal():
            b, xi = find_xi_g3(K, ideal_xi)
            if b:
                lst.append((xi, A))


    return lst


def imaginary_subfield_and_embedding(K):
    """
    Input: K = CM Field of degree 6
    Output: if K has an imaginary quadratic subfield k, then it returns
            a pair (k, m), where m is an embedding k --> K.
            Otherwise, returns None.

    """
    if K.degree() != 6:
        raise ValueError
    subfields = K.subfields()
    for s in subfields:
        if s[0].degree() == 2:
            return s[0], s[1]
            
            

def primitive_CM_type(K, xi=None, prec=None):
    """
    Returns a primitive CM-type Phi of K.
    If xi is specified, then the CM-type has to satisfy phi(xi) / i > 0 in RR
    for all phi in Phi.
    
    If prec is not None, then do an extra numerical sanity check to double-check
    primitivity of the CM-type.
    """
    cc = K.complex_conjugation()
    k, embedding_k = imaginary_subfield_and_embedding(K)
    x0 = embedding_k(k.gen())
    x0 = x0 - cc(x0)


    if xi is None:
        a = K.gen()
        xi = a - cc(a)
        if K.is_totally_positive_real_element(xi/x0) or K.is_totally_positive_real_element(-xi/x0):
            # This xi does not define a primitive CM-type.
            # To fix this: multiply xi by a totally real element with not all signs equal.
            u = K(K.real_generator_in_self())
            u = u - u.trace() / 6
            assert u.trace() == 0
            # Now u has positive and negative conjugates.
            xi = u * xi
            assert not (K.is_totally_positive_real_element(xi/x0) or K.is_totally_positive_real_element(-xi/x0))

    else:


        if not K.is_totally_real_element(xi/x0):
            raise ValueError, "Input xi is not totally positive imaginary"

        if K.is_totally_positive_real_element(xi/x0) or K.is_totally_positive_real_element(-xi/x0):
            raise ValueError, "Input xi does not define a primitive CM-type"

    Phi = CM_Type(xi)

    if not prec is None:
        # here are some old tests, which are numerical only
        # and did not work
        embeds = Phi.embeddings(prec)
        if embeds[0](x0) == embeds[1](x0) == embeds[2](x0):
            raise RuntimeError, "CM type is not primitive"

    return Phi



def period_matrices(K, prec, extra_info=False):
    """
    Chooses one primitive CM type, then computes one representative
    period matrix for every a.v. with CM by O_K of that CM type.
    
    IMPORTANT NOTE: 
    Uses g=3. May yield duplicates if N_{K/K_0}(O_K^*) is not O_{K_0}^2.
    May miss something if O_K^* is not O_{K_0}^*?
    See  Lemma 3.2.2 in Pinar's thesis: should be ok in the cases we're looking at.
    """
    lst = []
    Phi = primitive_CM_type(K, prec=prec)
    K0 = K.real_field()
    [u0,u1,u2] = [K.real_to_self(K0(u)) for u in K0.unit_group().gens()]
    units = [1, u0, u1, u0*u1, u2, u0*u2, u1*u2, u0*u1*u2]
    embs = Phi.embeddings(prec)

    for (xi0, aaa) in principally_polarized_ideal_classes(K):
        for u in units:
            xi = u*xi0
            if Phi.is_totally_positive_imaginary(xi):
                bas = vector(_symplectic_basis(aaa, xi, K.complex_conjugation()))
                (redbas, M) = reduce_Siegel_epsilon(bas, embs, True)
                bigmatrix = Matrix([[phi(b) for b in redbas] for phi in embs])
                Omega1 = bigmatrix[:,:3]
                Omega2 = bigmatrix[:,3:]
                tau = Omega2.inverse()*Omega1
                if extra_info:
                    lst.append((tau, aaa, xi0, u, xi, bas, redbas, M, bigmatrix))
                else:
                    lst.append(tau)
    return lst



    
def compute_period_matrices_of_the_list(list, prec=100):
    list_of_period_matrices = []
    R.<x> = QQ[]
    for poly in list:
        K = CM_Field(poly)
        list_of_period_matrices.append((poly, period_matrices(K, prec)))
    return list_of_period_matrices


