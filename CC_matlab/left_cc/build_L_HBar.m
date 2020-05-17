function [lambda1, lambda2] = build_L_HBar(L1, L2, HBar, omega, FM, VM, flag_ground, flag_jacobi)

    addpath(genpath('/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/utils'));
       
    [Nocc,Nunocc] = size(HBar{1}{1,2});
    Zocc = ones(Nocc)-eye(Nocc);
    Zunocc = ones(Nunocc)-eye(Nunocc);
    occ = 1:Nocc; unocc = Nocc+1:(Nocc+Nunocc);
    
    %eps_occ = diag(FM(occ,occ)); eps_unocc = diag(FM(unocc,unocc));

    % Lambda 1 equations
    % delta(mu,0)<0|HBar_op|{i,a}> + <0|(L1_op*HBar_op)_C|{i,a}> +
    % <0|(L2_op*HBar_op)_C|{i,a}> = omega_mu*<0|L1_op|{i,a}>
    
    if flag_ground == 1
        D1_L1 = HBar{1}{1,2};
    else
        D1_L1 = 0.0;
    end
        
    if flag_jacobi == 1
        D2_L1 = einsum(HBar{1}{2,2}.*Zunocc,L1,'ea,ie->ia');
        D3_L1 = -einsum(HBar{1}{1,1}.*Zocc,L1,'im,ma->ia');
        
        % MP partioning
%         D2L1 = einsum(HBar{1}{2,2} - diag(eps_unocc),L1,'ea,ie->ia');
%         D3L1 = -einsum(HBar{1}{1,1} - diag(eps_occ),L1,'im,ma->ia');
    else
        D2_L1 = einsum(HBar{1}{2,2},L1,'ea,ie->ia');
        D3_L1 = -einsum(HBar{1}{1,1},L1,'im,ma->ia');
    end 

    D4_L1 = 0.5*einsum(HBar{2}{2,2,2,1},L2,'efam,imef->ia');
    D5_L1 = -0.5*einsum(HBar{2}{1,2,1,1},L2,'iemn,mnae->ia');
    D6_L1 = einsum(HBar{2}{1,2,2,1},L1,'ieam,me->ia');
    D7_L1 =  0.25*einsum(HBar{3}{2,2,1,1,1,2},L2,'efimna,mnef->ia');
    
    
    lambda1 = D1_L1 + D2_L1 + D3_L1 + D4_L1 + D5_L1 + D6_L1 + D7_L1;
    %lambda1 = lambda1 - omega*lambda1;
    
        
    % LAMBDA 2 equations
    % delta(mu,0)<0|HBar_op|{i,j,a,b}> + <0|(L1_op*HBar_op)_C|{i,j,a,b}> +
    % <0|(L1_op*HBar_op)_DC|{i,j,a,b}> + <0|(L2_op*HBar_op)_C|{i,j,a,b}> =
    % omega_mu*<0|L2_op|{i,j,a,b}>
    
    if flag_ground == 1
        D1_L2 = HBar{2}{1,1,2,2};
    else
        D1_L2 = 0.0;
    end
        
    D2_L2 = einsum(HBar{2}{2,1,2,2},L1,'ejab,ie->ijab') - einsum(HBar{2}{2,1,2,2},L1,'eiab,je->ijab');
    D3_L2 = -einsum(HBar{2}{1,1,1,2},L1,'ijmb,ma->ijab') + einsum(HBar{2}{1,1,1,2},L1,'ijma,mb->ijab');
    
    D4_L2 = 0.5*einsum(HBar{2}{2,2,2,2},L2,'efab,ijef->ijab');
    D5_L2 = 0.5*einsum(HBar{2}{1,1,1,1},L2,'ijmn,mnab->ijab');
    
    D6_L2 = einsum(HBar{2}{1,2,2,1},L2,'jebm,imae->ijab') - ...
         einsum(HBar{2}{1,2,2,1},L2,'jeam,imbe->ijab') - ...
         einsum(HBar{2}{1,2,2,1},L2,'iebm,jmae->ijab') + ...
         einsum(HBar{2}{1,2,2,1},L2,'ieam,jmbe->ijab');
     
    D7_L2 = 0.5*einsum(HBar{3}{2,2,1,1,2,2},L2,'fejmab,imef->ijab') - ...
         0.5*einsum(HBar{3}{2,2,1,1,2,2},L2,'feimab,jmef->ijab');
         
    D8_L2 = -0.5*einsum(HBar{3}{2,1,1,1,1,2},L2,'eijnmb,mnae->ijab') + ...
          0.5*einsum(HBar{3}{2,1,1,1,1,2},L2,'eijnma,mnbe->ijab');
    
    if flag_jacobi == 1
        D9_L2 = einsum(HBar{1}{2,2}.*Zunocc,L2,'ea,ijeb->ijab') - ...
             einsum(HBar{1}{2,2}.*Zunocc,L2,'eb,ijea->ijab');
        D10_L2 = -einsum(HBar{1}{1,1}.*Zocc,L2,'im,mjab->ijab') + ...
               einsum(HBar{1}{1,1}.*Zocc,L2,'jm,miab->ijab');
           
        % MP paritioning
%         D9L2 = einsum(HBar{1}{2,2}-diag(eps_unocc),L2,'ea,ijeb->ijab') - ...
%              einsum(HBar{1}{2,2}-diag(eps_unocc),L2,'eb,ijea->ijab');
%         D10L2 = -einsum(HBar{1}{1,1}-diag(eps_occ),L2,'im,mjab->ijab') + ...
%                einsum(HBar{1}{1,1}-diag(eps_occ),L2,'jm,miab->ijab');
    else
        D9_L2 = einsum(HBar{1}{2,2},L2,'ea,ijeb->ijab') - ...
             einsum(HBar{1}{2,2},L2,'eb,ijea->ijab');
        D10_L2 = -einsum(HBar{1}{1,1},L2,'im,mjab->ijab') + ...
               einsum(HBar{1}{1,1},L2,'jm,miab->ijab');
    end
    
    % disconnected piece
    D11_L2 = einsum(HBar{1}{1,2},L1,'ia,jb->ijab') - ...
          einsum(HBar{1}{1,2},L1,'ja,ib->ijab') - ...
          einsum(HBar{1}{1,2},L1,'ib,ja->ijab') + ...
          einsum(HBar{1}{1,2},L1,'jb,ia->ijab');
    
    lambda2 = D1_L2 + D2_L2 + D3_L2 + D4_L2 + D5_L2 +...
              D6_L2 + D7_L2 + D8_L2 + D9_L2 + D10_L2 +...
              D11_L2;
          
    %lambda2 = lambda2 - omega*lambda2;
          
end

