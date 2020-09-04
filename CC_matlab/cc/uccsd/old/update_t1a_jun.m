import numpy as np


def t1a_update(f, v_a, v_b, occ, unocc, t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d):
    o = slice(0, occ)
    u = slice(occ, occ + unocc)
    q11 = np.einsum('me,ei->mi', f[o,u], t1a, optimize=True)
    q13 = np.einsum('mnie,en->im', v_a[o,o,o,u], t1a, optimize=True)
    q15 = np.einsum('maef,fm->ae', v_a[o,u,u,u], t1a, optimize=True)
    q17 = np.einsum('mnie,en->im', v_b[o,o,o,u], t1b, optimize=True)
    q19 = np.einsum('amef,fm->ae', v_b[u,o,u,u], t1b, optimize=True)
    q21 = np.einsum('mnef,aemn->fa', v_a[o,o,u,u], t2a, optimize=True)
    q23 = np.einsum('mnef,efin->mi', v_a[o,o,u,u], t2a, optimize=True)
    q25 = np.einsum('mnef,fn->em', v_a[o,o,u,u], t1a, optimize=True)
    q37 = np.einsum('em,ei->mi', q25, t1a, optimize=True)
    q29 = np.einsum('mnef,efin->mi', v_b[o,o,u,u], t2b, optimize=True)
    q31 = np.einsum('mnef,em->fn', v_b[o,o,u,u], t1a, optimize=True)
    q33 = np.einsum('mnef,fn->em', v_b[o,o,u,u], t1b, optimize=True)
    q39 = np.einsum('em,ei->mi', q33, t1a, optimize=True)
    q27 = np.einsum('mnef,afmn->ea', v_b[o,o,u,u], t2b, optimize=True)
    q35 = np.einsum('mnef,fn->em', v_a[o,o,u,u], t1b, optimize=True)
    x1 = np.einsum('mi->im', f[o,o]) + np.einsum('mi->im', q11) + np.einsum('im->im', q13) + np.einsum('im->im', q17) + 0.5*np.einsum('mi->im', q23) + np.einsum('mi->im', q37) + np.einsum('mi->im', q29) + np.einsum('mi->im', q39)
    z1 = np.einsum('im,am->ia', x1, t1a, optimize=True)
    x2 = np.einsum('ae->ae', f[u,u]) - np.einsum('ae->ae', q15) + np.einsum('ae->ae', q19) + 0.5*np.einsum('ea->ae', q21) - np.einsum('ea->ae', q27)
    z2 = np.einsum('ae,ei->ai', x2, t1a, optimize=True)
    z3 = np.einsum('maie,em->ia', v_a[o,u,o,u], t1a, optimize=True)
    z4 = np.einsum('ieam,em->ia', v_b[o,u,u,o], t1b, optimize=True)
    x3 = np.einsum('me->em', f[o,u]) + np.einsum('em->em', q25) + np.einsum('em->em', q33)
    z5 = np.einsum('em,aeim->ia', x3, t2a, optimize=True)
    z6 = np.einsum('mnie,aemn->ia', v_a[o,o,o,u], t2a, optimize=True)
    z7 = np.einsum('maef,efim->ai', v_a[o,u,u,u], t2a, optimize=True)
    x4 = np.einsum('me->em', f[o,u]) + np.einsum('em->em', q31) + np.einsum('em->em', q35)
    z8 = np.einsum('em,aeim->ia', x4, t2b, optimize=True)
    z9 = np.einsum('mnie,aemn->ia', v_b[o,o,o,u], t2b, optimize=True)
    z10 = np.einsum('amef,efim->ai', v_b[u,o,u,u], t2b, optimize=True)
    new_t = +np.einsum('ai->ai', f[u,o])
    new_t += -np.einsum('ia->ai', z1)
    new_t += +np.einsum('ai->ai', z2)
    new_t += -np.einsum('ia->ai', z3)
    new_t += +np.einsum('ia->ai', z4)
    new_t += +np.einsum('ia->ai', z5)
    new_t += -0.5*np.einsum('ia->ai', z6)
    new_t += -0.5*np.einsum('ai->ai', z7)
    new_t += +np.einsum('ia->ai', z8)
    new_t += -np.einsum('ia->ai', z9)
    new_t += +np.einsum('ai->ai', z10)

    for i in range(0,occ):
        for a in range(0,unocc):
            denom = +f[occ+a,occ+a]-f[i,i]
            coef = new_t[a,i] / denom
            t1a[a,i] = t1a[a,i] - coef
    return t1a

