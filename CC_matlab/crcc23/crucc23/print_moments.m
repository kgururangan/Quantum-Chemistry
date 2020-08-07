function [] = print_moments(MM,L,sys,cc_t,HBar_t,omega,i,j,k,a,b,c,SPIN)

    A = a; B = b; C = c;

    a = a-sys.Nocc_alpha;
    c = c-sys.Nocc_alpha;
    b = b-sys.Nocc_alpha;
    M3 = MM(a,b,c,i,j,k);
    L3 = L(i,j,k,a,b,c);

    H1A = HBar_t.H1A; H1B = HBar_t.H1B;
    H2A = HBar_t.H2A; H2B = HBar_t.H2B; H2C = HBar_t.H2C;

    [D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O] = get_triples_diagonal(cc_t,sys);

    if strcmp(SPIN,'A')

            DMP = sys.fa_vv(a,a) + sys.fa_vv(b,b) + sys.fa_vv(c,c) ...
                  -sys.fa_oo(i,i) - sys.fa_oo(j,j) - sys.fa_oo(k,k); 

            D1 = H1A.vv(a,a) + H1A.vv(b,b) + H1A.vv(c,c) ...
                -H1A.oo(i,i) - H1A.oo(j,j) - H1A.oo(k,k);

            D2 = -H2A.voov(a,i,i,a)-H2A.voov(b,i,i,b)-H2A.voov(c,i,i,c)...
                 -H2A.voov(a,j,j,a)-H2A.voov(b,j,j,b)-H2A.voov(c,j,j,c)...
                 -H2A.voov(a,k,k,a)-H2A.voov(b,k,k,b)-H2A.voov(c,k,k,c)...
                 -H2A.oooo(j,i,j,i)-H2A.oooo(k,i,k,i)-H2A.oooo(k,j,k,j)...
                 -H2A.vvvv(b,a,b,a)-H2A.vvvv(c,a,c,a)-H2A.vvvv(c,b,c,b);
            D2 = -D2;

            D3 = -D3A_O(a,i,j)-D3A_O(a,i,k)-D3A_O(a,j,k)...
                 -D3A_O(b,i,j)-D3A_O(b,i,k)-D3A_O(b,j,k)...
                 -D3A_O(c,i,j)-D3A_O(c,i,k)-D3A_O(c,j,k)...
                 +D3A_V(a,i,b)+D3A_V(a,i,c)+D3A_V(b,i,c)...
                 +D3A_V(a,j,b)+D3A_V(a,j,c)+D3A_V(b,j,c)...
                 +D3A_V(a,k,b)+D3A_V(a,k,c)+D3A_V(b,k,c);
    end

    if strcmp(SPIN,'B')
            DMP = sys.fa_vv(a,a) + sys.fa_vv(b,b) + sys.fb_vv(c,c) ...
                  -sys.fa_oo(i,i) - sys.fa_oo(j,j) - sys.fb_oo(k,k);

            D1 = H1A.vv(a,a) + H1A.vv(b,b) + H1B.vv(c,c) ...
                -H1A.oo(i,i) - H1A.oo(j,j) - H1B.oo(k,k);

            D2 = -H2A.voov(a,i,i,a)-H2A.voov(b,i,i,b)+H2B.vovo(c,i,c,i)...
                 -H2A.voov(a,j,j,a)-H2A.voov(b,j,j,b)+H2B.vovo(c,j,c,j)...
                 +H2B.ovov(k,a,k,a)+H2B.ovov(k,b,k,b)-H2C.voov(c,k,k,c)...
                 -H2A.oooo(j,i,j,i)-H2B.oooo(k,i,k,i)-H2B.oooo(k,j,k,j)...
                 -H2A.vvvv(b,a,b,a)-H2B.vvvv(c,a,c,a)-H2B.vvvv(c,b,c,b);
            D2 = -D2;

            D3 = -D3A_O(a,i,j)-D3B_O(a,i,k)-D3B_O(a,j,k)...
                 -D3A_O(b,i,j)-D3B_O(b,i,k)-D3B_O(b,j,k)...
                 -D3C_O(c,i,k)-D3C_O(c,j,k)...
                 +D3A_V(a,i,b)+D3B_V(a,i,c)+D3B_V(b,i,c)...
                 +D3A_V(a,j,b)+D3B_V(a,j,c)+D3B_V(b,j,c)...
                 +D3C_V(a,k,c)+D3C_V(b,k,c);
    end

    if strcmp(SPIN,'C')
            DMP = sys.fa_vv(a,a) + sys.fb_vv(b,b) + sys.fb_vv(c,c) ...
                  -sys.fa_oo(i,i) - sys.fb_oo(j,j) - sys.fb_oo(k,k);

            D1 = H1A.vv(a,a) + H1B.vv(b,b) + H1B.vv(c,c) ...
                -H1A.oo(i,i) - H1B.oo(j,j) - H1B.oo(k,k);

            D2 = -H2A.voov(a,i,i,a)+H2B.vovo(b,i,b,i)+H2B.vovo(c,i,c,i)...
                 +H2B.ovov(j,a,j,a)-H2C.voov(b,j,j,b)-H2C.voov(c,j,j,c)...
                 +H2B.ovov(k,a,k,a)-H2C.voov(b,k,k,b)-H2C.voov(c,k,k,c)...
                 -H2B.oooo(j,i,j,i)-H2B.oooo(k,i,k,i)-H2C.oooo(k,j,k,j)...
                 -H2B.vvvv(b,a,b,a)-H2B.vvvv(c,a,c,a)-H2C.vvvv(c,b,c,b);
            D2 = -D2;

            D3 = -D3B_O(a,i,j)-D3B_O(a,i,k)...
                 -D3C_O(b,i,j)-D3C_O(b,i,k)-D3C_O(b,j,k)...
                 -D3C_O(c,i,j)-D3C_O(c,i,k)-D3D_O(c,j,k)...
                 +D3B_V(a,i,b)+D3B_V(a,i,c)...
                 +D3C_V(a,j,b)+D3C_V(a,j,c)+D3D_V(b,j,c)...
                 +D3C_V(a,k,b)+D3C_V(a,k,c)+D3D_V(b,k,c);
    end

    if strcmp(SPIN,'D')
        DMP = sys.fb_vv(a,a) + sys.fb_vv(b,b) + sys.fb_vv(c,c) ...
              -sys.fb_oo(i,i) - sys.fb_oo(j,j) - sys.fb_oo(k,k);

        D1 = H1B.vv(a,a) + H1B.vv(b,b) + H1B.vv(c,c) ...
            -H1B.oo(i,i) - H1B.oo(j,j) - H1B.oo(k,k);

        D2 = -H2C.voov(a,i,i,a)-H2C.voov(b,i,i,b)-H2C.voov(c,i,i,c)...
             -H2C.voov(a,j,j,a)-H2C.voov(b,j,j,b)-H2C.voov(c,j,j,c)...
             -H2C.voov(a,k,k,a)-H2C.voov(b,k,k,b)-H2C.voov(c,k,k,c)...
             -H2C.oooo(j,i,j,i)-H2C.oooo(k,i,k,i)-H2C.oooo(k,j,k,j)...
             -H2C.vvvv(b,a,b,a)-H2C.vvvv(c,a,c,a)-H2C.vvvv(c,b,c,b);
        D2 = -D2;

        D3 = -D3D_O(a,i,j)-D3D_O(a,i,k)-D3D_O(a,j,k)...
             -D3D_O(b,i,j)-D3D_O(b,i,k)-D3D_O(b,j,k)...
             -D3D_O(c,i,j)-D3D_O(c,i,k)-D3D_O(c,j,k)...
             +D3D_V(a,i,b)+D3D_V(a,i,c)+D3D_V(b,i,c)...
             +D3D_V(a,j,b)+D3D_V(a,j,c)+D3D_V(b,j,c)...
             +D3D_V(a,k,b)+D3D_V(a,k,c)+D3D_V(b,k,c);
    end


    D_A = omega-DMP;
    D_B = omega-D1;
    D_C = omega-(D1+D2);
    D_D = omega-(D1+D2+D3);

    %occ_sort = sort([i,j,k],'ascend');
    %unocc_sort = sort([A,B,C],'ascend');
    occ_sort = [i,j,k];
    unocc_sort = [A,B,C];
    if strcmp(SPIN,'A')
        fprintf('%d  %d  %d->  %d %d %d\n',occ_sort(1),occ_sort(2),occ_sort(3),unocc_sort(1),unocc_sort(2),unocc_sort(3))
        fprintf('M3A*L3A = %4.12f\n',M3*L3)
        fprintf('M3A = %4.12f\n',M3)
        fprintf('L3A = %4.12f\n',L3)
        fprintf('D_A = %4.12f\n',D_A)
        fprintf('D_B = %4.12f\n',D_B)
        fprintf('D_C = %4.12f\n',D_C)
        fprintf('D_D = %4.12f\n\n',D_D)
    end
    if strcmp(SPIN,'B')
        fprintf('%d  %d  %d->  %d %d %d\n',occ_sort(1),occ_sort(2),occ_sort(3),unocc_sort(1),unocc_sort(2),unocc_sort(3))
        fprintf('M3B*L3B = %4.12f\n',M3*L3)
        fprintf('M3B = %4.12f\n',M3)
        fprintf('L3B = %4.12f\n',L3)
        fprintf('D_A = %4.12f\n',D_A)
        fprintf('D_B = %4.12f\n',D_B)
        fprintf('D_C = %4.12f\n',D_C)
        fprintf('D_D = %4.12f\n\n',D_D)
    end
    if strcmp(SPIN,'C')
        fprintf('%d  %d  %d->  %d %d %d\n',occ_sort(1),occ_sort(2),occ_sort(3),unocc_sort(1),unocc_sort(2),unocc_sort(3))
        fprintf('M3C*L3C = %4.12f\n',M3*L3)
        fprintf('M3C = %4.12f\n',M3)
        fprintf('L3C = %4.12f\n',L3)
        fprintf('D_A = %4.12f\n',D_A)
        fprintf('D_B = %4.12f\n',D_B)
        fprintf('D_C = %4.12f\n',D_C)
        fprintf('D_D = %4.12f\n\n',D_D)
    end
    if strcmp(SPIN,'D')
        fprintf('%d  %d  %d->  %d %d %d\n',occ_sort(1),occ_sort(2),occ_sort(3),unocc_sort(1),unocc_sort(2),unocc_sort(3))
        fprintf('M3D*L3D = %4.12f\n',M3*L3)
        fprintf('M3D = %4.12f\n',M3)
        fprintf('L3D = %4.12f\n',L3)
        fprintf('D_A = %4.12f\n',D_A)
        fprintf('D_B = %4.12f\n',D_B)
        fprintf('D_C = %4.12f\n',D_C)
        fprintf('D_D = %4.12f\n\n',D_D)
    end
end

