load('class_cmfield.sage')
load('class_cmpoint.sage')
load('thetafunctions.sage')

def mult_k(k, g2, g3, g4):
    assert k.is_square(), 'k must be a square'
    return (g2*k, g3*sqrt(k^3), g4*k^2)
        

def find_gs(j1, j2):
    g2 = j1
    g3 = j1^2
    g4 = j1^2*j2
    for g in [g2,g3,g4]:
        (g2, g3, g4) = mult_k(g.denom()*g.denom().squarefree_part(),g2,g3,g4)
        if all([g2 in ZZ, g3 in ZZ, g4 in ZZ]):
            break
    for fac in g2.factor():
        if fac[1] > 1:
            while all([a in ZZ for a in mult_k(fac[0]^(-2),g2,g3,g4)]):
                (g2, g3, g4) = mult_k(fac[0]^(-2),g2,g3,g4)
    return g2, g3, g4


R.<x> = QQ[]
poly = raw_input('Enter the cubic polynomial:')
K = CMField(R(poly), [3,0,0], 500)

CMtype = K.one_CMtype()
Js = [[],[]]
for pp in K.princ_polarized(CMtype):
    CMpoint = CMPoint(K, CMtype, pp[0], pp[1], is_Picard=True)
    CMpoint.find_endomorphism()
    j1,j2 = CMpoint.j_invariants()
    Js[0].append(j1)
    Js[1].append(j2)
    
if K.is_isomorphic(CyclotomicField(9)):
    print '(g2, g3, g4) = (0, -1, 0)'
else:     
    n = len(Js[0])
    CC = ComplexField(K._prec)
    Coeffs_Hj = [[CC((-1)^k*sum([prod(y) for y in Combinations(Js[i],k)])) for k in range(n+1)] for i in range(2)]
    HJs = map(lambda i : flatten(i.factor())[0:2*len(i.factor()):2], [sum([Coeffs_Hj[i][j].algdep(1).change_ring(QQ).roots()[0][0]*x^(n-j) for j in range(n+1)]) for i in range(2)])
    for case in range(len(HJs[0])):
        if HJs[0][case].degree() == 1:
            (j1, j2) = (HJs[0][case].roots()[0][0], HJs[1][case].roots()[0][0])
            (g2,g3,g4) = find_gs(j1,j2)
            print 'Case #'+str(case)+':\n (g2, g3, g4) = ('+', '.join([str(g2.factor()),str(g3.factor()),str(g4.factor())])+')'
        else:
            KJs = NumberField(HJs[0][case],'a')
            assert KJs.is_isomorphic(NumberField(HJs[1][case],'b'))
            if HJs[0][case].change_ring(K).is_irreducible():
                Fs = K.extension(HJs[0][case],'b').absolute_field('b').optimized_representation()[0].subfields(3)
                F.<g> = Fs[[fld[0].is_isomorphic(KJs) for fld in Fs].index(True)][0].change_names()
            else:
                F.<g> = K.subfields(3)[0][0].change_names()
            Js_alg = [HJs[i][case].change_ring(F).roots() for i in range(2)]
            (j1, j2) = (Js_alg[0][0][0], Js_alg[1][0][0])
            print 'Case #'+str(case)+':\n (j1, j2) = ('+str(j1)+', '+str(j2)+') where '+str(g.minpoly().subs(x = var('g')))
