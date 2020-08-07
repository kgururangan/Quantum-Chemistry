function [X1A] = build_HR_1A(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys)

    Nocc_a = sys.Nocc_alpha; Nunocc_a = sys.Nvir_alpha;

    H1A = HBar_t.H1A; H1B = HBar_t.H1B; H2A = HBar_t.H2A; H2B = HBar_t.H2B;

    X1A = zeros(Nunocc_a, Nocc_a);

    % we don't take out the MP denominator with Zocc/Zumocc, right???

    % (1)
    X1A = X1A - einsum_kg(H1A.oo,r1a,'mi,am->ai');
    % (2)
    X1A = X1A + einsum_kg(H1A.vv,r1a,'ae,ei->ai');
    % (3)
    X1A = X1A + einsum_kg(H2A.voov,r1a,'amie,em->ai');
    % (4)
    X1A = X1A + einsum_kg(H2B.voov,r1b,'amie,em->ai');
    % (5)
    X1A  = X1A - 0.5*einsum_kg(H2A.ooov,r2a,'mnif,afmn->ai');
    % (6)
    X1A = X1A - einsum_kg(H2B.ooov,r2b,'mnif,afmn->ai');
    % (7)
    X1A = X1A + 0.5*einsum_kg(H2A.vovv,r2a,'anef,efin->ai');
    % (8)
    X1A = X1A + einsum_kg(H2B.vovv,r2b,'anef,efin->ai');
    % (9)
    X1A = X1A + einsum_kg(H1A.ov,r2a,'me,aeim->ai');
    % (10)
    X1A = X1A + einsum_kg(H1B.ov,r2b,'me,aeim->ai');

    

end

