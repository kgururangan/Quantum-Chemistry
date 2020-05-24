module gaussian_integration
    implicit none

    real(16), parameter :: pi = 4 * atan (1.0_16)

    contains

        function overlap_integral(alpha, rA, beta, rB) result(sAB)

            real, intent(in) :: alpha, beta
            real, dimension(1:3), intent(in) :: rA, rB
            real :: norm_factor
            real :: sAB

            norm_factor = (2.0*alpha/pi)**(3.0/4.0)*(2.0*beta/pi)**(3.0/4.0)

            sAB = norm_factor*(pi/(alpha+beta))**(3.0/2.0)*exp(-alpha*beta/(alpha+beta)*sum((rA-rB)**2.0))

        end function overlap_integral

        function KE_integral(alpha, rA, beta, rB) result(tAB)

            real, intent(in) :: alpha, beta
            real, dimension(1:3), intent(in) :: rA, rB
            real :: x1, x2, norm_factor
            real :: tAB


            norm_factor = (2.0*alpha/pi)**(3.0/4.0)*(2.0*beta/pi)**(3.0/4.0)
            x1 = alpha*beta/(alpha+beta)
            x2 = sum((rA-rB)**2.0)
            tAB = norm_factor*x1*(3.0-2.0*x1*x2)*(pi/(alpha+beta))**(3.0/2.0)*exp(-x1*x2)
        end function KE_integral

        function PE_integral(alpha,rA,beta,rB,rC,zC) result(vAB)
            use utils_fcn, only: boys

            real, intent(in) :: alpha, beta
            real, intent(in) :: zC
            real, dimension(1:3), intent(in) :: rA, rB, rC

            real, dimension(1:3) :: rP
            real :: p, x1, rAB, rPC, norm_factor
            real :: vAB

            norm_factor = (2.0*alpha/pi)**(3.0/4.0)*(2.0*beta/pi)**(3.0/4.0)
            p = alpha + beta
            rP = (alpha*rA + beta*rB)/p 
            x1 = alpha*beta/p 
            rAB = sum((rA-rB)**2.0)
            rPC = sum((rP-rC)**2.0)

            vAB = -norm_factor*2.0*pi/p*zC*exp(-x1*rAB)*boys(p*rPC)

        end function PE_integral

        function two_body_integral(xi1,rA,xi2,rB,xi3,rC,xi4,rD) result(vABCD)
            use utils_fcn, only : boys

            real, intent(in) :: xi1, xi2, xi3, xi4
            real, dimension(1:3), intent(in) :: rA, rB, rC, rD

            real, dimension(1:3) :: rP, rQ
            real :: p, q, x1, x2, x3, rAB, rCD, rPQ, norm_factor
            real :: vABCD

            norm_factor = (2.0*xi1/pi)**(3.0/4.0)*(2.0*xi2/pi)**(3.0/4.0)*(2.0*xi3/pi)**(3.0/4.0)*(2.0*xi4/pi)**(3.0/4.0)

            p = xi1 + xi2
            q = xi3 + xi4
            rP = (xi1*rA + xi2*rB)/p 
            rQ = (xi3*rC + xi4*rD)/q 

            x1 = xi1*xi2/p 
            x2 = xi3*xi4/q
            x3 = p*q/(p+q) 
            rAB = sum((rA-rB)**2.0)
            rCD = sum((rC-rD)**2.0)
            rPQ = sum((rP-rQ)**2.0) 

            vABCD = norm_factor*(2.0*pi**(5.0/2.0))/(p*q*sqrt(p+q))*exp(-x1*rAB-x2*rCD)*boys(x3*rPQ)

        end function two_body_integral
end module gaussian_integration
