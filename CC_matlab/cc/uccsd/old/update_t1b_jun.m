import numpy as np


def t1b_update(f, v_a, v_b, occ, unocc, t1a, t1b, t2a, t2b, t2c, t3a, t3b, t3c, t3d):
    o = slice(0, occ)
    u = slice(occ, occ + unocc)
    q11 = np.einsum('mnei,em->in', v_b[o,o,u,o], t1a, optimize=True)
    q13 = np.einsum('maef,em->af', v_b[o,u,u,u], t1a, optimize=True)
    q15 = np.einsum('mnef,em->fn', v_a[o,o,u,u], t1a, optimize=True)
    q17 = np.einsum('mnef,em->fn', v_b[o,o,u,u], t1a, optimize=True)
    q37 = np.einsum('fn,fi->ni', q17, t1b, optimize=True)
    q21 = np.einsum('mnie,en->im', v_a[o,o,o,u], t1b, optimize=True)
    q23 = np.einsum('maef,fm->ae', v_a[o,u,u,u], t1b, optimize=True)
    q25 = np.einsum('mnef,eamn->fa', v_b[o,o,u,u], t2b, optimize=True)
    q27 = np.einsum('mnef,efmi->ni', v_b[o,o,u,u], t2b, optimize=True)
    q29 = np.einsum('mnef,fn->em', v_b[o,o,u,u], t1b, optimize=True)
    q31 = np.einsum('mnef,aemn->fa', v_a[o,o,u,u], t2c, optimize=True)
    q33 = np.einsum('mnef,efin->mi', v_a[o,o,u,u], t2c, optimize=True)
    q35 = np.einsum('mnef,fn->em', v_a[o,o,u,u], t1b, optimize=True)
    q39 = np.einsum('em,ei->mi', q35, t1b, optimize=True)
    q19 = np.einsum('me,ei->mi', f[o,u], t1b, optimize=True)
    z1 = np.einsum('maei,em->ia', v_b[o,u,u,o], t1a, optimize=True)
    x1 = np.einsum('mi->im', f[o,o]) + np.einsum('im->im', q11) + np.einsum('mi->im', q37) + np.einsum('im->im', q21) + np.einsum('mi->im', q27) + 0.5*np.einsum('mi->im', q33) + np.einsum('mi->im', q39) + np.einsum('mi->im', q19)
    z2 = np.einsum('im,am->ia', x1, t1b, optimize=True)
    x2 = np.einsum('ae->ae', f[u,u]) + np.einsum('ae->ae', q13) - np.einsum('ae->ae', q23) - np.einsum('ea->ae', q25) + 0.5*np.einsum('ea->ae', q31)
    z3 = np.einsum('ae,ei->ai', x2, t1b, optimize=True)
    z4 = np.einsum('maie,em->ia', v_a[o,u,o,u], t1b, optimize=True)
    x3 = np.einsum('me->em', f[o,u]) + np.einsum('em->em', q15) + np.einsum('em->em', q29)
    z5 = np.einsum('em,eami->ia', x3, t2b, optimize=True)
    z6 = np.einsum('mnei,eamn->ia', v_b[o,o,u,o], t2b, optimize=True)
    z7 = np.einsum('maef,efmi->ai', v_b[o,u,u,u], t2b, optimize=True)
    x4 = np.einsum('me->em', f[o,u]) + np.einsum('em->em', q17) + np.einsum('em->em', q35)
    z8 = np.einsum('em,aeim->ia', x4, t2c, optimize=True)
    z9 = np.einsum('mnie,aemn->ia', v_a[o,o,o,u], t2c, optimize=True)
    z10 = np.einsum('maef,efim->ai', v_a[o,u,u,u], t2c, optimize=True)
    new_t = +np.einsum('ai->ai', f[u,o])
    new_t += +np.einsum('ia->ai', z1)
    new_t += -np.einsum('ia->ai', z2)
    new_t += +np.einsum('ai->ai', z3)
    new_t += -np.einsum('ia->ai', z4)
    new_t += +np.einsum('ia->ai', z5)
    new_t += -np.einsum('ia->ai', z6)
    new_t += +np.einsum('ai->ai', z7)
    new_t += +np.einsum('ia->ai', z8)
    new_t += -0.5*np.einsum('ia->ai', z9)
    new_t += -0.5*np.einsum('ai->ai', z10)

    for i in range(0,occ):
        for a in range(0,unocc):
            denom = +f[occ+a,occ+a]-f[i,i]
            coef = new_t[a,i] / denom
            t1b[a,i] = t1b[a,i] - coef
    return t1b

