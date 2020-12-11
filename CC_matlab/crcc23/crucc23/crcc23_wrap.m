function [Ecrcc23A, Ecrcc23B, Ecrcc23C, Ecrcc23D] = crcc23_wrap(cc_t,HBar_t,sys,omega)
    fprintf('\n==================================++Entering CR-CC(2,3) Routine++=============================\n')

    Ecorr_ccsd = ucc_energy(cc_t.t1a,cc_t.t1b,cc_t.t2a,cc_t.t2b,cc_t.t2c,sys);

    if nargin < 4
        fprintf('\n')
        fprintf('-------------------------------------------------------------------------------------\n')
        fprintf('Performing correction for ground state -\n')
        [deltaA,deltaB,deltaC,deltaD] = crucc23(cc_t,HBar_t,sys);
        Ecorr_crcc23A = Ecorr_ccsd + deltaA;
        Ecorr_crcc23B = Ecorr_ccsd + deltaB;
        Ecorr_crcc23C = Ecorr_ccsd + deltaC;
        Ecorr_crcc23D = Ecorr_ccsd + deltaD;
        Ecrcc23A(1) = sys.Escf + Ecorr_crcc23A;
        Ecrcc23B(1) = sys.Escf + Ecorr_crcc23B;
        Ecrcc23C(1) = sys.Escf + Ecorr_crcc23C;
        Ecrcc23D(1) = sys.Escf + Ecorr_crcc23D;
        fprintf('\n')
        fprintf('CR-UCC(2,3)_A = %4.12f Eh     Ecorr_A = %4.12f Eh     Delta_A = %4.12f Eh\n',Ecrcc23A(1),Ecorr_crcc23A,deltaA)
        fprintf('CR-UCC(2,3)_B = %4.12f Eh     Ecorr_B = %4.12f Eh     Delta_B = %4.12f Eh\n',Ecrcc23B(1),Ecorr_crcc23B,deltaB)
        fprintf('CR-UCC(2,3)_C = %4.12f Eh     Ecorr_C = %4.12f Eh     Delta_C = %4.12f Eh\n',Ecrcc23C(1),Ecorr_crcc23C,deltaC)
        fprintf('CR-UCC(2,3)_D = %4.12f Eh     Ecorr_D = %4.12f Eh     Delta_D = %4.12f Eh\n',Ecrcc23D(1),Ecorr_crcc23D,deltaD)
        fprintf('-------------------------------------------------------------------------------------\n')
    else
        Ecrcc23A = zeros(1,length(omega)+1);
        Ecrcc23B = zeros(1,length(omega)+1);
        Ecrcc23C = zeros(1,length(omega)+1);
        Ecrcc23D = zeros(1,length(omega)+1);
    

        for iroot = 0:length(omega)

            if iroot == 0
                fprintf('\n')
                fprintf('-------------------------------------------------------------------------------------\n')
                fprintf('Performing correction for ground state -\n')
                [deltaA,deltaB,deltaC,deltaD] = crucc23(cc_t,HBar_t,sys);
                Ecorr_crcc23A = Ecorr_ccsd + deltaA;
                Ecorr_crcc23B = Ecorr_ccsd + deltaB;
                Ecorr_crcc23C = Ecorr_ccsd + deltaC;
                Ecorr_crcc23D = Ecorr_ccsd + deltaD;
                Ecrcc23A(1) = sys.Escf + Ecorr_crcc23A;
                Ecrcc23B(1) = sys.Escf + Ecorr_crcc23B;
                Ecrcc23C(1) = sys.Escf + Ecorr_crcc23C;
                Ecrcc23D(1) = sys.Escf + Ecorr_crcc23D;
                fprintf('\n')
                fprintf('CR-UCC(2,3)_A = %4.12f Eh     Ecorr_A = %4.12f Eh     Delta_A = %4.12f Eh\n',Ecrcc23A(1),Ecorr_crcc23A,deltaA)
                fprintf('CR-UCC(2,3)_B = %4.12f Eh     Ecorr_B = %4.12f Eh     Delta_B = %4.12f Eh\n',Ecrcc23B(1),Ecorr_crcc23B,deltaB)
                fprintf('CR-UCC(2,3)_C = %4.12f Eh     Ecorr_C = %4.12f Eh     Delta_C = %4.12f Eh\n',Ecrcc23C(1),Ecorr_crcc23C,deltaC)
                fprintf('CR-UCC(2,3)_D = %4.12f Eh     Ecorr_D = %4.12f Eh     Delta_D = %4.12f Eh\n',Ecrcc23D(1),Ecorr_crcc23D,deltaD)
                fprintf('-------------------------------------------------------------------------------------\n')
            else
                fprintf('\n')
                fprintf('-------------------------------------------------------------------------------------\n')
                fprintf('Performing correction for root %d -\n',iroot)

                %[deltaA,deltaB,deltaC,deltaD] = creomucc23(cc_t,omega,HBar_t,sys,iroot);
                [deltaA,deltaB,deltaC,deltaD] = creomucc23_rhf(cc_t,omega,HBar_t,sys,iroot);

                omega_crcc23A = omega(iroot) + deltaA;
                omega_crcc23B = omega(iroot) + deltaB;
                omega_crcc23C = omega(iroot) + deltaC;
                omega_crcc23D = omega(iroot) + deltaD;

                Etot_eomccsd = sys.Escf + Ecorr_ccsd + omega(iroot);
                Ecrcc23A(iroot+1) = sys.Escf + Ecorr_ccsd + omega_crcc23A;
                Ecrcc23B(iroot+1) = sys.Escf + Ecorr_ccsd + omega_crcc23B;
                Ecrcc23C(iroot+1) = sys.Escf + Ecorr_ccsd + omega_crcc23C;
                Ecrcc23D(iroot+1) = sys.Escf + Ecorr_ccsd + omega_crcc23D;
                
                VEE_A = Ecrcc23A(iroot+1)-Ecrcc23A(1);
                VEE_B = Ecrcc23B(iroot+1)-Ecrcc23B(1);
                VEE_C = Ecrcc23C(iroot+1)-Ecrcc23C(1);
                VEE_D = Ecrcc23D(iroot+1)-Ecrcc23D(1);

                HatoeV = 27.2113957;

                fprintf('\n')
                fprintf('EOMCCSD = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n\n',Etot_eomccsd,omega(iroot),omega(iroot)*HatoeV)
                fprintf('CR-EOMCC(2,3)_A = %4.12f Eh     Delta_A = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Ecrcc23A(iroot+1),deltaA,VEE_A,VEE_A*HatoeV)
                fprintf('CR-EOMCC(2,3)_B = %4.12f Eh     Delta_B = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Ecrcc23B(iroot+1),deltaB,VEE_B,VEE_B*HatoeV)
                fprintf('CR-EOMCC(2,3)_C = %4.12f Eh     Delta_C = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Ecrcc23C(iroot+1),deltaC,VEE_C,VEE_C*HatoeV)
                fprintf('CR-EOMCC(2,3)_D = %4.12f Eh     Delta_D = %4.12f Eh     VEE = %4.12f Eh   (%4.8f eV)\n',Ecrcc23D(iroot+1),deltaD,VEE_D,VEE_D*HatoeV)
                fprintf('-------------------------------------------------------------------------------------\n')
            end

        end

    end
end
