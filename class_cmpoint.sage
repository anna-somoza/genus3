load("https://bitbucket.org/pkilicer/period-matrices-for-genus-3-cm-curves/raw/master/period_matrices_genus3.sage")

def CMPoint(K, CM_type, ideal, xi, prec=None, is_hyperelliptic=False, is_Picard=False):
    """
    Return a CMPoint object, which can be hyperelliptic, Picard or generic depending on input.

    TODO: Add functionality to automatically recognize to which category it belongs.

    INPUT:
        K: a CMField (see class CMFieldfromPoly)
        CM_type: can be generated from the CMField, see documentation.
        ideal: can be generated from the CMField, see documentation.
        xi: can be generated from the CMField, see documentation.
        prec : bits of precision for complex-valued calculations, if None, inherited from CMField.
        is_hyperelliptic (default: False): The CMPoint is hyperelliptic. Returns a CMPoint_hyperelliptic object.
        is_Picard (default: False): The CMPoint is Picard. Returns a CMPoint_picard object.
    """
    if is_hyperelliptic:
        return CMPoint_hyperelliptic(K, CM_type, ideal, xi, prec)
    if is_Picard:
        return CMPoint_picard(K, CM_type, ideal, xi, prec)
    return CMPoint_generic(K, CM_type, ideal, xi, prec)

class CMPoint_generic:
    def __init__(self, K, CM_type, ideal, xi, prec=None):
        """
        K: a CMField (see class CMFieldfromPoly)
        CM_type: can be generated from the CMField, see documentation
        ideal: can be generated from the CMField, see documentation
        xi: can be generated from the CMField, see documentation
        prec : bits of precision for complex-valued calculations, if None, inherited from CMField
        """
        self._K = K
        self._CM_type = CM_type
        self._ideal = ideal
        self._xi = xi
        if prec == None:
            self._prec = self._K._prec
        else:
            self._prec = prec

    def __repr__(self):
        return "CM point for CM field %s \n given by the ideal %s \n and the CM-type %s,\n with polarization given by %s"%(self._K, self._ideal, self._CM_type, self._xi)

    def cubic_fld(self):
        """
        Gives CMpoint its cubic field K0, which is a subfield of the optimal sextic
        """
        self._K0 = self._K._K0
        return self._K0

    def sextic_fld(self):
        """
        Gives CMpoint its sextic field
        """
        return self._K

    def period_matrix(self, test=False, as_tuple=True):
        """
        Computes the period matrix associated with the polarized abelian variety.
        INPUT:
            - Test (default: False) Set to True to see the eigenvalues of the imaginary part of the period matrix
        NOTE: Modified version of the original code by BILV to include LLL reduction method and keep track of the algebraic expression of the matrix if possible.
        """
        try:
            Z = self._period_matrix
            if test:
                CMvalimag = Z.apply_map(lambda i: i.imag())
                print CMvalimag.eigenvalues()
            if as_tuple:
                try:
                    assert len(Z)== 6
                    return Z
                except TypeError:
                    A= Z[0][0]
                    B= Z[0][1]
                    C= Z[0][2]
                    D= Z[1][1]
                    E= Z[1][2]
                    F= Z[2][2]
                    return [A, B, C, D, E, F]
            else:
                try:
                    assert len(Z) == 6
                    return matrix([[Z[0],Z[1],Z[2]],[Z[1],Z[3],Z[4]],[Z[2],Z[4],Z[5]]])
                except TypeError:
                    return Z
        except AttributeError:
            K = self._K
            Phi = self._CM_type
            ideal = self._ideal
            xi = self._xi
            basis = ideal.basis()
            riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
            symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))

            red_basis, red_relation = reduce_Siegel_epsilon(vector(symplectic_basis), Phi, False)
            assert is_symplectic(red_relation), 'The matrix of the basis change is not symplectic'
            big_period_matrix = Matrix([[phi(b) for b in red_basis] for phi in Phi])
            big_period_matrix.subdivide(3,3)
            Omega1 = big_period_matrix.subdivision(0,0)
            Omega2 = big_period_matrix.subdivision(0,1)
            self._A2 = Omega2
            Z = Omega2.adjoint()*Omega1/Omega2.det()
            if test:
                print Z
                CMvalimag = Z.apply_map(lambda i: i.imag())
                print CMvalimag.eigenvalues()

            if as_tuple:
                A= Z[0][0]
                B= Z[0][1]
                C= Z[0][2]
                D= Z[1][1]
                E= Z[1][2]
                F= Z[2][2]

                M =[A, B, C, D, E, F]
                self._period_matrix = M
                return self._period_matrix
            else:
                self._period_matrix = Z
                return self._period_matrix

    def A2(self):
        try:
            return self._A2
        except AttributeError:
            self.period_matrix()
            return self._A2

    def acc_period_matrix(self, test = False):
        """
        This computes the period matrix several times, refining the precision of the CM type until two consecutive period matrices have entries that agree up to the precision of the CM point
        """
        K = self._K
        ideal = self._ideal
        xi = self._xi
        prec = self._prec
        from sage.rings.number_field.number_field import refine_embedding

        basis = ideal.basis()
        riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
        symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
        phis = self._CM_type
        big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
        big_period_matrix.subdivide(3,3)
        Omega1 = big_period_matrix.subdivision(0,0)
        Omega2 = big_period_matrix.subdivision(0,1)
        Zs = [Omega2.adjoint()*Omega1/Omega2.det()]

        equality = False
        iterates = 0

        while equality == False:
            iterates += 1
            phis = [refine_embedding(phi,prec + iterates*20) for phi in phis]
            big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
            big_period_matrix.subdivide(3,3)
            Omega1 = big_period_matrix.subdivision(0,0)
            Omega2 = big_period_matrix.subdivision(0,1)
            Zs.append(Omega2.adjoint()*Omega1/Omega2.det())

            if test:
                print "computing iteration number {0}".format(iterates)
                print Zs[iterates]
                CMvalimag=matrix(RR,3,3)
                for i in range(3):
                    for j in range(3):
                        CMvalimag[i,j] = Zs[iterates][i,j].imag()
                print CMvalimag.eigenvalues()

            if all([compare(Zs[iterates][i,j],Zs[iterates-1][i,j],prec+10) for i in range(3) for j in range(3)]):
                equality = True

        Z = Zs[iterates]

        CC = ComplexField(prec)
        A= CC(Z[0][0])
        B= CC(Z[0][1])
        C= CC(Z[0][2])
        D= CC(Z[1][1])
        E= CC(Z[1][2])
        F= CC(Z[2][2])

        M =[A, B, C, D, E, F]
        self._period_matrix = M
        return self._period_matrix

class CMPoint_hyperelliptic(CMPoint_generic):
    def __repr__(self):
        return "Hyperelliptic CM point for CM field %s \n given by the ideal %s \n and the CM-type %s,\n with polarization given by %s"%(self._K, self._ideal, self._CM_type, self._xi)
    def all_thetas(self, start_bound = 20, prec = None, bound = False):
        try:
            period_matrix = self._period_matrix
        except:
            period_matrix = self.acc_period_matrix()
        if prec == None:
            prec = self._prec

        all_evens = [[[0,0,0],[0,0,0]],[[1,0,0],[0,0,0]],[[0,1,0],[0,0,0]],[[0,0,1],[0,0,0]],[[0,0,1],[1,0,0]],[[0,0,1],[0,1,0]],[[1,1,0],[0,0,0]],[[1,0,1],[0,0,0]],[[0,1,1],[0,0,0]],[[0,0,0],[1,0,1]],[[0,0,0],[0,1,1]],[[0,0,0],[1,1,0]],[[1,1,1],[0,0,0]],[[0,0,0],[1,1,1]],[[0,1,1],[1,0,0]],[[1,0,1],[0,1,0]],[[1,1,0],[0,0,1]],[[0,1,0],[1,0,1]],[[0,0,1],[1,1,0]],[[1,0,0],[0,1,1]],[[1,0,1],[1,0,1]],[[1,1,0],[1,1,0]],[[0,1,1],[0,1,1]],[[1,0,1],[1,1,1]],[[1,1,0],[1,1,1]],[[1,1,1],[0,1,1]],[[1,1,1],[1,0,1]],[[1,1,1],[1,1,0]],[[0,1,1],[1,1,1]],[[0,0,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[0,0,1]],[[1,0,0],[0,1,0]],[[1,0,0],[0,0,1]],[[0,1,0],[0,0,1]],[[0,1,0],[1,0,0]]]

        all_values = [[[val[0][0][1],val[0][0][2]],val[1]] for val in sorted(list(theta_function([(period_matrix,even[0],even[1],2,prec) for even in all_evens])))]

        self._all_thetas = all_values
        return self._all_thetas

    def counting(self, epsilon = 10.^(-2),bound = False):
        """
        this is meant to be used on the output of cmpoint.all_thetas() to count how many theta values are zero
        epsilon can be changed depending on what value one wants to consider to be "zero"
        """
        count = 0
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)
        for value in all_values:
            if value[1].abs() < epsilon:
                count += 1
                return count

    def vanishing_char(self, bound = False, epsilon = 10.^(-2)):
        """
        inputs:
        bound is passed to all_thetas
        epsilon can be changed depending on what value one wants to consider to be "zero"
        outputs:
        if the period matrix is plane quartic, returns None
        if the period matrix is hyperelliptic, returns the theta characteristic "delta" such that theta[delta](Z) = 0
        if there are more than one vanishing characteristics, raises an error
        """
        count = 0
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(bound = bound)

        for value in all_values:
            if value[1].abs() < epsilon:
                count += 1
        if count == 0:
            return None
        elif count > 1:
            raise TypeError('The entries of this period matrix are too large, the theta functions don\'t converge well')
        elif count == 1:
            for value in all_values:
                if value[1].abs() < epsilon:
                    self._vanishing_char = value[0]
                    return self._vanishing_char
    def eta_dict(self, bound = False, epsilon = 10.^(-2)):
        """
        bound is passed to all_thetas (ultimately) True is theta_with_bound and False is theta_without_bound
        returns a dictionary giving values eta_1, eta_2, ... eta_7 for an eta-map associated to the period matrix computed for this cm point
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound, epsilon = epsilon)
        if vanishing_char == None:
            raise TypeError('This is a plane quartic Jacobian')
        else:
            delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])

        if delta == matrix(GF(2)):
            self._eta_dict = eta_bar
            return self._eta_dict
        else:
            delta.set_immutable()
            M = pairs[delta]
            self._eta_dict = {1: M*mumford_eta[1], 2: M*mumford_eta[2], 3: M*mumford_eta[3], 4: M*mumford_eta[4], 5:M*mumford_eta[5], 6: M*mumford_eta[6], 7:M*mumford_eta[7]}
            return self._eta_dict

    def U_set(self, bound = False, epsilon = 10.^(-2)):
        """
        returns U = {2, 4, 6} if delta is non zero and U = {1, 2, 3, 4, 5, 6, 7} if delta is zero (infinity is implicitly in both sets)
        """
        try:
            vanishing_char = self._vanishing_char
        except:
            vanishing_char = self.vanishing_char(bound = bound, epsilon = epsilon)
        if vanishing_char == None:
            raise TypeError('This is a plane quartic Jacobian')
        else:
            delta = matrix(GF(2),[[vanishing_char[0][0]],[vanishing_char[0][1]],[vanishing_char[0][2]],[vanishing_char[1][0]],[vanishing_char[1][1]],[vanishing_char[1][2]]])

        if delta == matrix(GF(2)):
            self._U_set = Set([1,2,3,4,5,6,7])
            return self._U_set
        else:
            self._U_set = Set([2,4,6])
            return self._U_set

    def vareps(self, j, bound = False, epsilon = 10.^(-2)):
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound, epsilon = epsilon)
        v1 = eta_dict[1]-eta_dict[2]
        v2 = eta_dict[2]
        J = matrix(GF(2),[[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
        value = v1.transpose()*J*v2
        if value == 0:
            return 1
        elif value == 1:
            return -1
        

    def one_rosenhain_coeff(self, j, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        """
        This function computes the Rosenhain coefficient j, with 3 <= j <= 7. We assume a1 = 0, a2 = 1. Three values are computed, for the three ways in which you can split the set in Takase's paper. This serves as an additional check that the period matrix is truly hyperelliptic.
        """
        try:
            eta_dict = self._eta_dict
        except:
            eta_dict = self.eta_dict(bound = bound, epsilon = epsilon)
        try:
            U = self._U_set
        except:
            U = self.U_set(bound = bound, epsilon = epsilon)
        if prec == None:
            prec = self._prec
        try:
            all_values = self._all_thetas
        except:
            all_values = self.all_thetas(start_bound, prec, bound)


        T = Set([1,2,3,4,5,6,7])
        S = Set([1,2,j])
        Y = T.difference(S)
        X = Y.subsets(2)

        ajvec = []
        # we introduce an auxiliary variable a, this only makes it so we compute the value for the decomposition given by V and W once (rather than twice, for V and W interchanged)
        if (j >= 3) and (j <= 6):
            a = j + 1
        elif j == 7:
            a = 3
        else:
            raise ValueError('j is not between 3 and 7 inclusively')

        for V in X:
            if a in V:
                W = Y.difference(V)
                setA = U.symmetric_difference(V.union(Set([1,2])))
                setB = U.symmetric_difference(W.union(Set([1,2])))
                setC = U.symmetric_difference(V.union(Set([1,j])))
                setD = U.symmetric_difference(W.union(Set([1,j])))
                A = theta_from_char_and_list(all_values, eta_dict, compute_characteristic_sum_from_set_and_etas(setA,eta_dict))
                B = theta_from_char_and_list(all_values, eta_dict, compute_characteristic_sum_from_set_and_etas(setB,eta_dict))
                C = theta_from_char_and_list(all_values, eta_dict, compute_characteristic_sum_from_set_and_etas(setC,eta_dict))
                D = theta_from_char_and_list(all_values, eta_dict, compute_characteristic_sum_from_set_and_etas(setD,eta_dict))
                aj = self.vareps(j,bound,epsilon)*((A*B)^2)/((C*D)^2)
                ajvec.append(aj)
        return ajvec

    def all_rosenhain_coeffs(self, prec = None, start_bound = 20, bound = False, epsilon = 10.^(-2)):
        """
        returns all of the vectors of Rosenhain coefficients
        """
        if prec == None:
            prec = self._prec

        all_coeffs = []
        for j in range(3,8):
            all_coeffs.append(self.one_rosenhain_coeff(j,prec, start_bound, bound, epsilon))
        return all_coeffs

class CMPoint_picard(CMPoint_generic):
    def __repr__(self):
        return "Picard CM point for CM field %s \n given by the ideal %s \n and the CM-type %s,\n with polarization given by %s"%(self._K, self._ideal, self._CM_type, self._xi)

    def find_endomorphism(self):
        K = self._K
        CM_type = self._CM_type
        Psi = K.complex_embeddings(self._prec)[0]
        roots = K.roots_of_unity()
        L = [a.multiplicative_order() for a in roots]
        for i in range(len(L)):
            if L[i] == 3 and Psi(roots[i]).imag() > 0:
                zeta = roots[i]
                break

        Omega = self.period_matrix(as_tuple=False)
        big_period_matrix = block_matrix([[Omega, identity_matrix(3)], [Omega.conjugate(), identity_matrix(3)]], subdivide=False)
        big_A2 = block_diagonal_matrix(self.A2(), (self.A2()).conjugate(), subdivide=False)
        big_zeta = diagonal_matrix([phi(zeta) for phi in CM_type]+[phi(zeta).conjugate() for phi in CM_type])

        M = ((big_period_matrix.inverse())*(big_A2.inverse())*big_zeta*big_A2*big_period_matrix).apply_map(lambda i: round(i.real()))
        assert compare(norm(M-((big_period_matrix.inverse())*(big_A2.inverse())*big_zeta*big_A2*big_period_matrix)), 0.0), 'M is not integer'

        M = M.transpose()

        self._M = M

        MA = M[0:3,0:3];
        MB = M[0:3,3:6];
        MC = M[3:6,0:3];
        MD = M[3:6,3:6];

        assert compare(norm(Omega - (MA*Omega + MB)*((MC*Omega + MD).inverse())), 0.0), 'Omega does not meet M·Om = Om'


    def Riemann_constant(self):
        """
        Computes the Riemann constant associated to the rational representation of the morphism z(x,y) = (x, xi*y), where xi is the third root of unity.
        """
        try:
            return self._Riemann_constant
        except AttributeError:
            M = self._M

            MA = M[0:3,0:3];
            MB = M[0:3,3:6];
            MC = M[3:6,0:3];
            MD = M[3:6,3:6];

            MDelta = matrix(FiniteField(2), identity_matrix(6) - M.transpose().inverse())
            v = vector(FiniteField(2),(MC*(MD.transpose())).diagonal()+(MA*(MB.transpose())).diagonal())
            Delta = 1/2*vector(QQ,MDelta.solve_right(-v).list())
            self._Riemann_constant = Delta
            return self._Riemann_constant

    def all_taus(self):
        """
        Computes the characteristics that correspond to the 3-torsion points of the curve with jacobian the principally polarized abelian variety
        """
        try:
            return self._taus
        except AttributeError:
            M = self._M
            MTaus = matrix(FiniteField(3), identity_matrix(6) - M.transpose().inverse())
            W = MTaus.right_kernel().list()
            Taus = [1/3*vector(QQ, t.list()) for t in W]
            assert len(Taus) <= 27, 'Too many taus'
            assert len(Taus) >= 27, 'Too few taus'
            self._taus = Taus
            return self._taus

    def zero_locus(self, tope = 3, tol = 1e-5):
        """
        Returns a list of indexes of the Tau elements that are in the zero locus of the theta function.
        """
        try:
            return self._TauZeros
        except AttributeError:
            try:
                Delta = self._Riemann_constant
            except AttributeError:
                Delta = self.Riemann_constant()
            try:
                Taus = self._taus
            except AttributeError:
                Taus = self.all_taus()
            period_matrix = self.period_matrix(as_tuple=True)

            all_chars = [[6*(tau + Delta)[:3],6*(tau + Delta)[3:]] for tau in Taus]
            TauZeros = []
            for val in sorted(list(theta_function([(period_matrix,characteristic[0],characteristic[1],6,50) for characteristic in all_chars]))):
                if abs(val[1]) < tol:
                    i = all_chars.index([val[0][0][1],val[0][0][2]])
                    if i != 0:
                        TauZeros.append(i)
            

            assert len(TauZeros) >= 14, 'Too few zeros'
            assert len(TauZeros) <= 14, 'Too many zeros'
            self._zero_locus= TauZeros
            return self._zero_locus

    def find_ds(self):
        """
        Finds a set of 4 taus such that each three are lineraly independent but their sum is 0.
        """
        try:
            Tau = self._taus
        except AttributeError:
            Tau = self.all_taus()
        try:
            TauZeros = self._zero_locus
        except AttributeError:
            TauZeros = self.zero_locus()

        DZeros = [Tau[aux] for aux in TauZeros]
        for [a,b,c] in Combinations(TauZeros,3).list():
            forth = ((-Tau[a]-Tau[b]-Tau[c]).apply_map(lambda i : i-numerator(i)//denominator(i)))
            LI = are_li(Tau[a],Tau[b],Tau[c]);
            if forth in DZeros and LI:
                d = TauZeros[DZeros.index(forth)]
                if (d > c):
                    assert all((Tau[a]+Tau[b]+Tau[c]+Tau[d]).apply_map(lambda i : 1*(i  in ZZ))), 'sum not zero'
                    return (a,b,c,d)
        return None

    def j_invariants(self, tope = 20):
        """
        Finds the j-invariants of the curve with jacobian the principally polarized abelian variety.
        INPUT:
            - tope: bound on the theta function sum
        """
        try:
            TauZeros = self._zero_locus
        except AttributeError:
            TauZeros = self.zero_locus()
        Taus = self._taus
        Delta = self._Riemann_constant
        period_matrix = self.period_matrix(as_tuple=True)
        DZeros = [Taus[aux] for aux in TauZeros]

        prec = self._prec

        d_list = self.find_ds()
        D = [Taus[aux] for aux in d_list]
        zTC = [0,0,0,0]
        TC = [0,0,0,0]

        zTC[0] = (-Delta + D[1] + 2*D[2] - D[0])
        zTC[1] = (-Delta + 2*D[1] + D[2] - D[0])
        zTC[2] = (-Delta + D[1] + 2*D[3] - D[0])
        zTC[3] = (-Delta + 2*D[1] + D[3] - D[0])

        for val in sorted(list(theta_function([(period_matrix,6*characteristic[:3],6*characteristic[3:],6,prec) for characteristic in zTC]))):
            i = zTC.index(1/6*vector(list(val[0][0][1])+list(val[0][0][2])))
            TC[i] = val[1]
            
        assert all(TC), 'There is one theta constant missing'

        LAM = (TC[0]/TC[1])^3*exp(6*pi*I*((zTC[1]-zTC[0])[:3]*D[0][3:]+(-Delta + 2*D[1] + D[2])[:3]*(2*Delta - 3*D[1] - 3*D[2])[3:])).n(prec)
        MU = (TC[2]/TC[3])^3*exp(6*pi*I*((zTC[2]-zTC[0])[:3]*D[0][3:]+(-Delta + 2*D[1] + D[3])[:3]*(2*Delta - 3*D[1] - 3*D[3])[3:])).n(prec)

        C.<X> = PolynomialRing(ComplexField(prec))
        P = X*(X-1)*(X-LAM)*(X-MU);
        P2=P.subs(X = X-P[3]/4);
        J1 = P2[1]^2/P2[2]^3;
        J2 = P2[0]/P2[2]^2;

        return (J1,J2)

def are_li(x,y,z):
    """
    Given three 3-torsion points, checks whether they are linearly independents or not.
    TODO: Change to matrix rank mod 3.
    """
    for a in xrange(3):
        for b in xrange(3):
            for c in xrange(3):
                if (a,b,c) != (0,0,0) and all((a*x + b*y + c*z).apply_map(lambda i : 1*(i in ZZ)).list()):
                    return False
    return True

def is_symplectic(M):
    JJ = block_matrix([[0,identity_matrix(3)],[-identity_matrix(3),0]])
    return M.transpose()*JJ*M == JJ
