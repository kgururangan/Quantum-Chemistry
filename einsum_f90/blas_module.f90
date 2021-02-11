module blas_module

      implicit none

      contains

              function kgemv(A,x) result(c)
                        
                      real, intent(in) :: A(:,:), x(:)
                      real :: c(size(x,1))
                      integer :: M, N

                      M = size(A,1)
                      N = size(A,2)

                      call dgemv('n',size(A,1),size(A,2),1.0,A,size(A,1),x,1,0.0,c,1)

              end function kgemv

              function kgemm(A,B) result(C)

                     ! integer, intent(in) :: dimA_1, dimA_2, dimB_1, dimB_2
                      real, intent(in) :: A(:,:), B(:,:)
                      real :: C(size(A,1),size(B,2))
                      
                     ! if (size(A,1) /= size(B,1)) then
                     !   print*, 'DGEMM error: Inconsistent contraction dimensions'
                     ! else
                        call dgemm('n','n',size(A,1),size(B,2),size(A,2),1.0,&
                                    A,size(A,1),B,size(B,1),0.0,C,size(A,1))

                     ! endif

               end function kgemm

                REAL FUNCTION DDOT(N,DX,INCX,DY,INCY)
                        
                       INTEGER :: INCX, INCY, N
                       REAL :: DX(*), DY(*)
                       REAL :: DTEMP
                       INTEGER :: I, IX, IY, M, MP1

                       INTRINSIC MOD

                       DDOT = 0.0d0
                       DTEMP = 0.0d0
                       IF (N.LE.0) RETURN
                       IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
                               M = MOD(N,5)
                               IF (M.NE.0) THEN
                                  DO I = 1,M
                                       DTEMP = DTEMP + DX(I)*DY(I)
                                  ENDDO
                                  IF (N.LT.5) THEN
                                       DDOT=DTEMP
                                       RETURN
                                  ENDIF
                                ENDIF
                                MP1 = M + 1
                                DO I = MP1,N,5
                                   DTEMP = DTEMP + DX(I)*DY(I)+DX(I+1)*DY(I+1) +&
                                           DX(I+2)*DY(I+2)+DX(I+3)*DY(I+3)+&
                                           DX(I+4)*DY(I+4)
                                ENDDO
                        ELSE
                                IX = 1
                                IY = 1
                                IF (INCX.LT.0) THEN 
                                        IX = (-N+1)*INCX+1
                                ENDIF
                                IF (INCY.LT.0) THEN
                                        IY = (-N+1)*INCY+1
                                ENDIF
                                DO I = 1,N
                                        DTEMP = DTEMP+DX(IX)*DY(IY)
                                        IX = IX + INCX
                                        IY = IY + INCY
                                ENDDO
                         ENDIF
                         DDOT = DTEMP
                         RETURN
               END FUNCTION DDOT

                      
               function kgdot(a,b) result (c)
                       !integer, intent(in) :: n
                       real, intent(in) :: a(:), b(:)
                       real :: c
                       c = ddot(size(a),a,1,b,1)
               end function kgdot


               subroutine kgeig(A,VR,VL,WR,WI)

                       real, intent(in) :: A(:,:)
                       real, intent(out) :: VR(size(A,1),size(A,1)), VL(size(A,1),size(A,1))
                       real, intent(out) :: WR(size(A,1)), WI(size(A,1))
                       integer :: info
                       real :: work(5*size(A,1))

                       call dgeev('V','V',size(A,1),A,size(A,1),WR,WI,VL,size(A,1),VR,size(A,1),work,5*size(A,1),info)

!                       if (info > 0) then
!                              print*,'QR algorithm failed at',info,'root.'&
!                                      'eigenvectors not calculated'
!                       elseif (info < 0) then
!                              print*, 'The',info,'argument has an illegal value'
!                       endif
               end subroutine kgeig 

               subroutine kgeig_sym(A,W)
                ! On entry, A contains the symmetric matrix
                ! on exit, A contains the eigenvectors

                      !integer, intent(in) :: n
                       real, intent(inout) :: A(:,:)
                       real, intent(out) :: W(size(A,1))
                       real :: work(5*size(A,1))
                       integer :: info

                       call dsyev('V','U',size(A,1),A,size(A,1),W,work,5*size(A,1),info)

              end subroutine kgeig_sym

              function matinv(A) result(Ainv)
                          real, dimension(:,:), intent(in) :: A
                          real, dimension(size(A,1),size(A,2)) :: Ainv

                          real, dimension(size(A,1)) :: work  ! work array for LAPACK
                          integer, dimension(size(A,1)) :: ipiv   ! pivot indices
                          integer :: n, info

                          ! External procedures defined in LAPACK
                          external DGETRF
                          external DGETRI

                          ! Store A in Ainv to prevent it from being overwritten by LAPACK
                          Ainv = A
                          n = size(A,1)

                          ! DGETRF computes an LU factorization of a general M-by-N matrix A
                          ! using partial pivoting with row interchanges.
                          call DGETRF(n, n, Ainv, n, ipiv, info)

                          if (info /= 0) then
                             stop 'Matrix is numerically singular!'
                          end if

                          ! DGETRI computes the inverse of a matrix using the LU factorization
                          ! computed by DGETRF.
                          call DGETRI(n, Ainv, n, ipiv, work, n, info)

                          if (info /= 0) then
                             stop 'Matrix inversion failed!'
                          end if

              end function matinv

 end module blas_module
