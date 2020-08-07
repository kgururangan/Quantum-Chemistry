function [X1B] = build_HR_1B(r1a,r1b,r2a,r2b,r2c,HBar_t,cc_t,sys)

    Nocc_b = sys.Nocc_beta; Nunocc_b = sys.Nvir_beta;

    H1A = HBar_t.H1A; H1B = HBar_t.H1B; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

    X1B = zeros(Nunocc_b,Nocc_b);

    % (1)
    X1B = X1B - einsum_kg(H1B.oo,r1b,'mi,am->ai');
    % (2)
    X1B = X1B + einsum_kg(H1B.vv,r1b,'ae,ei->ai');
    % (3)
    X1B = X1B + einsum_kg(H2B.ovvo,r1a,'maei,em->ai');
    % (4)
    X1B = X1B + einsum_kg(H2C.voov,r1b,'amie,em->ai');
    % (5)
    X1B = X1B - einsum_kg(H2B.oovo,r2b,'nmfi,fanm->ai');
    % (6)
    X1B = X1B - 0.5*einsum_kg(H2C.ooov,r2c,'mnif,afmn->ai');
    % (7)
    X1B = X1B + einsum_kg(H2B.ovvv,r2b,'nafe,feni->ai');
    % (8)
    X1B = X1B + 0.5*einsum_kg(H2C.vovv,r2c,'anef,efin->ai');
    % (9)
    X1B = X1B + einsum_kg(H1A.ov,r2b,'me,eami->ai');
    % (10)
    X1B = X1B + einsum_kg(H1B.ov,r2c,'me,aeim->ai');

end