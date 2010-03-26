!!$   Notes:
!!$
!!$   C++ will pass transposed matrices to Fortran. This is good so that the inner index (the spin
!!$   of the MPS) will become the last index in Fortran, and code will run faster because I can multiply
!!$   inner indices (bond indices) and sum over external indices.
!!$   IMPORTANT: The matrix multiplications will need to have the opposite order to compensate for the
!!$   fact that matrices are transposed:
!!$                    (L A R)^T = R^T A^T L^T and these are the matrices I have
!!$    
!!$   Mathematica has a way of passing complex numbers to C++ as 64 bit reals. Since C++ does not have
!!$   a way to handle complex, maybe I am better off passing first real part and then complex part and
!!$   adding them up inside Fortran, which can be compiled to do it fast (maybe even with threads).
!!$   This should be a minor issue anyway since the matrices are not likely to be more than 100x100
!!$   (but about 20 or 30 of them).
!!$
!!$   I HAVE TO ACCOMODATE LEFT AND RIGHT MPS TENSORS FROM THE START.
!!$   Maybe for this it will be good to have sizeL and sizeR differentiated in my tensors, as if I could
!!$   Have an inhomogeneous (in bond dimension) MPS.

module ArpackDriver

  implicit none

!!$c Common variables for the whole module

!!$c Parameters
  complex,parameter :: II=(0.0d0,1.0d0)
  real(8),parameter :: sqrt12=0.7071067811865476D0  !1/sqrt(2)
  complex(8),parameter :: one=1.0d0,zero=0.0d0
  character*1,parameter :: normal='N'  !This specifies that the matrix multiplication is normal and not transposed (using a trick actually)

!!$c     %---------------------------------------------------%
!!$c     | Define maximum dimensions for all arrays.         |
!!$c     | MAXN:   Maximum dimension of the matrix allowed.  |
!!$c     | MAXNEV: Maximum eigenvectors allowed              |
!!$c     | MAXNCV: Maximum Arnoldi vectors allowed           |
!!$c     %---------------------------------------------------%

  integer,parameter :: maxnev=4, maxncv=16

!!$c These are used to shape and reshape the tensor MPS. It is assumed the MPS has one physical index and two bond indices
  integer VectorShape(1),MatrixShape(3)    
  complex(8) sigma(2,2,3)  !Pauli matrices, sigma(:,:,0) is identity and the rest are X,Y, and Z

  contains

!!$ The subroutine expects 14+6r matrices, where r is the range of the interaction in the chain
!!$ The matrices are (where C denotes Chi the bond dimension and a is replaced by L or R, left and right sizes)
!!$    Left , Right, left and right tensors with no operators: dimensions (Ca,Ca)
!!$
!!$    DLeft and DRight are dipolar tensors that are added: dimensions (Ca,Ca,3)
!!$
!!$    hLeft and hRight are the field terms: dimensions (Ca,Ca,3)
!!$
!!$    vLeft and vRight are the list of matrices that are used for interaction of site with rest of chain
!!$                      dimensions (Ca,Ca,ra,3)
!!$ 
!!$    Ham is part of the Hamiltonian matrix: dimensions (rL+1,rL+rR+1,3)
!!$
!!$ 
!!$    The previous tensor can be passed to accelerate convergence of Arnoldi
!!$
!!$    A dimensions (CR,CL,spin) Notice the transposition of dimensions (this program accomodates only spin=2 matrices)
!!$
!!$    EVERYTHING IS PASSED IN REAL AND IMAGINARY PARTS, SO FIRST THING IS TO SUM EVERYTHING INTO COMPLEX
!!$    FORTRAN OBJECTS FIRST.
!!$
!!$    The result is placed in the Ar and Ai tensors.
!!$    Communication back to C++ and Mathematica is done by the MyInfo variable.
!!$    Here it will be said wether ARPACK was successful or not, and what caused the error.
!!$    In the future maybe implement returning the error of the vectors.
  !!  C++ call:
  !!  arpackdriver_mp_computegroundstate_(Ar,Ai,Leftr,Lefti,Rightr,Righti,DLeftr,DLefti,DRightr,DRighti,
  !!	  hLeftr, hLefti, hRightr, hRighti,vLeftr, vLefti, vRightr, vRighti,
  !!      Ham, &spin, &sizeL, &sizeR, &rangeL, &rangeR, &info);
	  
    subroutine ComputeGroundState(Ar,Ai,DLeftr,Dlefti,DRightr,DRighti, &
         & hLeftr,hLefti,hRightr,hRighti,vLeftr,vLefti,vRightr,vRighti,Ham,spin,sizeL,sizeR,rLeft,rRight,MyInfo,energy)
      integer spin,sizeL,sizeR,rLeft,rRight,MyInfo
      real(8) Ar(sizeR,sizeL,spin),Ai(sizeR,sizeL,spin)
      real(8) DLeftr(sizeL,sizeL,3),Dlefti(sizeL,sizeL,3),DRightr(sizeR,sizeR,3),DRighti(sizeR,sizeR,3)
      real(8) hLeftr(sizeL,sizeL,3),hLefti(sizeL,sizeL,3),hRightr(sizeR,sizeR,3),hRighti(sizeR,sizeR,3)
      real(8) vLeftr(sizeL,sizeL,rLeft,3),vLefti(sizeL,sizeL,rLeft,3),vRightr(sizeR,sizeR,rRight,3),vRighti(sizeR,sizeR,rRight,3)
      real(8) Ham(rLeft+rRight+1,rLeft+1,3),energy

      !Now declare the same arrays but complex
      complex(8) A(sizeR,sizeL,spin)
      complex(8) DLeft(sizeL,sizeL,3),DRight(sizeR,sizeR,3)
      complex(8) hLeft(sizeL,sizeL,3),hRight(sizeR,sizeR,3)
      complex(8) vLeft(sizeL,sizeL,rLeft,3),vRight(sizeR,sizeR,rRight,3)

      !This will be the vector form of A
      !complex(8) Avector(spin*sizeL*sizeR)

      !These are variables and vectors used by ARPACK
      integer iparam(11), ipntr(14)
      logical select(maxncv)
      complex(8)  d(maxncv),v(spin*sizeL*sizeR,maxncv) 
      complex(8) workd(3*spin*sizeL*sizeR), workev(3*maxncv), resid(spin*sizeL*sizeR),workl(3*maxncv*maxncv+5*maxncv)
      real(8) rwork(maxncv), rd(maxncv,3)

      !Local scalars, some used by ARPACK
      character*2 which
      character*1 bmat
      integer ido, n, nx, nev, ncv, lworkl, info, j,ierr, nconv, maxitr, ishfts, mode
      real(8) tol
      logical rvec

      integer Lv  !This will be the size of the input tensor

      Lv=spin*sizeL*sizeR
      VectorShape(1) = Lv
      MatrixShape = (/sizeR,sizeL,spin/)

      !First things first: convert all input to complex numbers. Use Fortran 90 parallel computation
      A=Ar+II*Ai
      DLeft=DLeftr+II*Dlefti
      DRight=DRightr+II*DRighti
      hLeft=hLeftr+II*hLefti
      hRight=hRightr+II*hRighti
      vLeft=vLeftr+II*vLefti
      vRight=vRightr+II*vRighti

	      
!!$   ARPACK requires certain parameters
!!$   BMAT='I' means a normal eigenvalue problem (not generalized)      
!!$   NEV is the number of eigenvalues to be approximated.
!!$   NCV is the number of Arnoldi vectors used.
!!$   WHICH is just set to compute the smallest eigenvalue by real part
!!$ 
!!$   ARPACK efficiency changes with NEV and NCV, so one must tune this.
!!$   However, The following conditions must be satisfied:                    |
!!$             N <= MAXN                      
!!$             NEV <= MAXNEV                    
!!$             NEV + 2 <= NCV <= MAXNCV              
!!$   A rule of thumb is that NCV=4*NEV
      nev   = 1
      ncv   = 4 
      bmat  = 'I'
      which = 'SR' 
!!$   From ARPACK example documentation:
!!$c     %---------------------------------------------------%
!!$c     | The work array WORKL is used in ZNAUPD  as         | 
!!$c     | workspace.  Its dimension LWORKL is set as        |
!!$c     | illustrated below.  The parameter TOL determines  |
!!$c     | the stopping criterion. If TOL<=0, machine        |
!!$c     | precision is used.  The variable IDO is used for  |
!!$c     | reverse communication, and is initially set to 0. |
!!$c     | Setting INFO=0 indicates that a random vector is  |
!!$c     | generated to start the ARNOLDI iteration.         | 
!!$c     %---------------------------------------------------%
!!$c
!!$   We will set info to 1 so that we use the entered value of the MPS tensor as the starting point
!!$   This should accelerate convergence after many sweeps of the algorithm
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 1

      resid=Reshape(A,VectorShape)

!!$   From ARPACK example documentation:
!!$c     %---------------------------------------------------%
!!$c     | This program uses exact shift with respect to     |
!!$c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!!$c     | IPARAM(3) specifies the maximum number of Arnoldi |
!!$c     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!!$c     | (IPARAM(7) = 1). All these options can be changed |
!!$c     | by the user. For details see the documentation in |
!!$c     | ZNAUPD .                                           |
!!$c     %---------------------------------------------------%
!!$c
      ishfts = 1
      maxitr = 300
      mode   = 1
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 

!!$      Done with initialization. Now start the main loop

!%!%!-----------------------------------------------------------------
!%!%!------------- MAIN LOOP
!%!%!-----------------------------------------------------------------
      do !This is a do block without variable, must exit explicitely
         
         !! Repeatedly call the routine ZNAUPD  and take actions indicated by parameter IDO until
         !! either convergence is indicated or maxitr has been exceeded.                      
         call znaupd  ( ido, bmat, Lv, which, nev, tol, resid, ncv,v, Lv, iparam, ipntr, workd, workl, lworkl,rwork,info )
         !!  If ido is equal to 1 or -1 then we need to provide the matrix multiplication, otherwise exit
         if (ido .ne. -1 .and. ido .ne. 1) then
            exit
         endif
         !!  When arpack needs it, 
         !!  Perform vector multiplication routine that takes workd(ipntr(1)) as the input vector, 
         !! and returns the matrix vector product to workd(ipntr(2)).          
         call  MatrixTimesVector(workd(ipntr(1)),workd(ipntr(2)),DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham,spin,Lv,sizeL,sizeR,rLeft,rRight)

      enddo
!%!%!-----------------------------------------------------------------
!%!%!------------- END MAIN LOOP
!%!%!-----------------------------------------------------------------

!! At this point there is convergence or error

      if ( info .lt. 0 ) then  
         !! ERROR
         MyInfo=info  !Return the ARPACK error code
	 Ar=0.0d0
         Ai=0.0d0
	 energy=0.0d0
         return
      else  
         !! NO ERRORS
         !! Post-Process using ZNEUPD
         !! Extract eigenvalues and compute Eigenvectors (indicated by rvec = .true.)
         rvec = .true.
         call zneupd (rvec,'A',select,d,v,Lv,sigma,workev,bmat,Lv,which,nev,tol,resid,ncv,v,Lv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)

         !! Eigenvalue will be discarded, ground-state is in first column of array v.
         !! We have to reshape it as an MPS tensor and return it
         A=Reshape(v(:,1),MatrixShape)

         Ar=Real(A)
         Ai=aImag(A)
         energy=Real(d(1))
	      
         !! Return convoluted ARPACK error code
         MyInfo=info+100*ierr  

      end if

    end subroutine ComputeGroundState






!%!%!-----------------------------------------------------------------
!%!%!-------------  MATRIX VECTOR MULTIPLICATION ROUTINE
!%!%!-----------------------------------------------------------------
!!    This routine takes a MPS site tensor and multiplies it by the
!!    effective Hamiltonian at that site given by all the D, v and h matrices.
!!    It is assumed that the state is canonized to the left and right of the
!!    site, which makes it easier because the L and R matrices are identity.
!!
!!    Notice that multiplications L*A*R are actually computed as R'*A'*L'
!!    Because the calling from C++ transposes the matrices and because in 
!!    Fortran it makes more sense to have the spin index in the last position
!!    because of the Fortran memory storage (column major, first index changes faster)
!!
    subroutine MatrixTimesVector(input,output,DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham,spin,Lv,sizeL,sizeR,rLeft,rRight)
      integer spin,Lv,sizeL,sizeR,rLeft,rRight
      complex(8) input(Lv),output(Lv)
      complex(8) DLeft(sizeL,sizeL,3),DRight(sizeR,sizeR,3)
      complex(8) hLeft(sizeL,sizeL,3),hRight(sizeR,sizeR,3)
      complex(8) vLeft(sizeL,sizeL,rLeft,3),vRight(sizeR,sizeR,rRight,3)
      real(8) Ham(rLeft+rRight+1,rLeft+1,3)
      ! Internal matrices
      complex(8) temp(sizeR,sizeL),A(sizeR,sizeL,spin),Aop(sizeR,sizeL,spin),out(sizeR,sizeL,spin),identity(sizeR,sizeL)
      integer axis,x,y,j,s
      ! BLAS Level 3 function to compute matrix multiplication
      external ZGEMM

      ! Initialize internal matrices
      A=Reshape(input,(/sizeR,sizeL,spin/))
      out=0.0d0      

      !First add up terms that do not have an operator at the site
       do s=1,spin
          do axis=1,3
            !Multiply L * (I*A) * hR[axis]
            call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,hRight(:,:,axis),sizeR,A(:,:,s),sizeR,one,out(:,:,s),sizeR)
 
            !Multiply L * (I*A) * DR[axis]
            call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,DRight(:,:,axis),sizeR,A(:,:,s),sizeR,one,out(:,:,s),sizeR)
 
            !Multiply hL[axis] * (I*A) * Right
            call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,A(:,:,s),sizeR,hLeft(:,:,axis),sizeL,one,out(:,:,s),sizeR)

            !Multiply DL[axis] * (I*A) * Right
            call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,A(:,:,s),sizeR,DLeft(:,:,axis),sizeL,one,out(:,:,s),sizeR)

            do x=1,rLeft
               do y=1,rRight
                  !Check for distance, has to be less than maximum 
                  if ((y+x-1).lt.max(rLeft,rRight)) then
                    !Multiply vLeft[axis,0-x] * (I*A) * vRight[axis,0+y] * Ham[rL+1-x,rL+1+y]
                    !First do vRight' *  A' and store in temp
                    call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,vRight(:,:,y,axis),sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
                    !Now add to out computing out=temp*vL'+out
					temp=temp*Ham(rLeft+1+y,rLeft+1-x,axis)
                    call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,vLeft(:,:,x,axis),sizeL,one,out(:,:,s),sizeR)
                  endif
               enddo
            enddo

         enddo
      enddo

      !Now add up terms that DO have an operator at the site
	!Define Pauli matrices
	sigma=0.0d0
	sigma(1,2,1)=0.5d0
	sigma(2,1,1)=0.5d0
	sigma(1,2,2)=-0.5d0*II
	sigma(2,1,2)=0.5d0*II
	sigma(1,1,3)=0.5d0
	sigma(2,2,3)=-0.5d0
	do axis=1,3
         ! Multiply matrix by operator
		Aop=0.0d0
         do s=1,spin
            temp=0.0d0
            do j=1,spin
               temp=temp+sigma(s,j,axis)*A(:,:,j)
            enddo
            Aop(:,:,s)=temp
         enddo
         !Multiply L * (Sigma[axis]*A) * R * Ham(site,site)
         out=out+Ham(rLeft+1,rLeft+1,axis)*Aop

         do s=1,spin
            
            do x=1,rLeft
               !Multiply vLeft[axis,rL+1-x] * (Sigma[axis]*A) * Right * Ham(rL+1,rL+1-x)
!				This is how this function should be, but it does not work and I don't know why
!               call ZGEMM(normal,normal,sizeR,sizeL,sizeL,Ham(rLeft+1,rLeft+1-x,axis),Aop(:,:,s),sizeR,vLeft(:,:,x,axis),sizeL,one,out(:,:,s),sizeR)
				!But this is how it works
				call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,Aop(:,:,s),sizeR,vLeft(:,:,x,axis),sizeL,zero,temp,sizeR)
				out(:,:,s)=out(:,:,s)+temp*Ham(rLeft+1,rLeft+1-x,axis)
            enddo

            do x=1,rRight   
               !Multiply vLeft[axis,0-x] * (Sigma[axis]*A) * Right 
         		call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,vRight(:,:,x,axis),sizeR,Aop(:,:,s),sizeR,zero,temp,sizeR)
                out(:,:,s)=out(:,:,s)+temp*Ham(rLeft+1+x,rLeft+1,axis)
!				This is how this function should be, but it does not work and I don't know why
!         		call ZGEMM(normal,normal,sizeR,sizeL,sizeR,Ham(rLeft+1+x,rLeft+1,axis),vRight(:,:,x,axis),sizeR,Aop(:,:,s),sizeR,one,out(:,:,s),sizeR)
            enddo

          enddo

       enddo

        !Reshape into vector form				  
      output=Reshape(out,(/Lv/))
	  
			  
    end subroutine MatrixTimesVector




!!   This is mostly a helper routine. It is used from C so that Mathematica can call the matrix vector
!!   multiplication routine directly without going through the ground state calculation.
!!   It takes separate real and imaginary parts, combines them, and calls the multiplication routine above.
!!   Returns separate real and imag parts.
    subroutine MatrixTimesVectorFromC(inputr,inputi,outputr,outputi,DLeftr,DLefti,DRightr,DRighti, & 
	    & hLeftr,hLefti,hRightr, hRighti,vLeftr, vLefti,vRightr, vRighti,Ham,spin,Lv,sizeL,sizeR,rLeft,rRight)
	integer spin,Lv,sizeL,sizeR,rLeft,rRight
	real(8) inputr(sizeR,sizeL,spin),inputi(sizeR,sizeL,spin)
	real(8) outputr(sizeR,sizeL,spin),outputi(sizeR,sizeL,spin)
	real(8) DLeftr(sizeL,sizeL,3),Dlefti(sizeL,sizeL,3),DRightr(sizeR,sizeR,3),DRighti(sizeR,sizeR,3)
	real(8) hLeftr(sizeL,sizeL,3),hLefti(sizeL,sizeL,3),hRightr(sizeR,sizeR,3),hRighti(sizeR,sizeR,3)
	real(8) vLeftr(sizeL,sizeL,rLeft,3),vLefti(sizeL,sizeL,rLeft,3),vRightr(sizeR,sizeR,rRight,3),vRighti(sizeR,sizeR,rRight,3)
	real(8) Ham(rLeft+rRight+1,rLeft+1,3)

	complex(8) input(Lv),output(Lv)
	complex(8) DLeft(sizeL,sizeL,3),DRight(sizeR,sizeR,3)
	complex(8) hLeft(sizeL,sizeL,3),hRight(sizeR,sizeR,3)
	complex(8) vLeft(sizeL,sizeL,rLeft,3),vRight(sizeR,sizeR,rRight,3)

	input=Reshape(inputr+II*inputi,(/Lv/))
	output=Reshape(outputr+II*outputi,(/Lv/))
        DLeft=DLeftr+II*Dlefti
	DRight=DRightr+II*DRighti
	hLeft=hLeftr+II*hLefti
	hRight=hRightr+II*hRighti
	vLeft=vLeftr+II*vLefti
	vRight=vRightr+II*vRighti
	
	call MatrixTimesVector(input,output,DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham,spin,Lv,sizeL,sizeR,rLeft,rRight)
	
	outputr=real(Reshape(output,(/sizeR,sizeL,spin/)))
	outputi=aimag(Reshape(output,(/sizeR,sizeL,spin/)))
		
    end subroutine MatrixTimesVectorFromC

end module ArpackDriver


!! The following is just some legacy code, left here in case some checking is needed in the future.
!! The only diference with the above code is that it includes the L and R matrices that are forced to be unity
!! So that this code is slower by a large factor.

!     subroutine MatrixTimesVectorGOOD(input,output,DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham,spin,Lv,sizeL,sizeR,rLeft,rRight)
!       integer spin,Lv,sizeL,sizeR,rLeft,rRight
!       complex(8) input(Lv),output(Lv)
!       complex(8) Left(sizeL,sizeL),Right(sizeR,sizeR)
!       complex(8) DLeft(sizeL,sizeL,3),DRight(sizeR,sizeR,3)
!       complex(8) hLeft(sizeL,sizeL,3),hRight(sizeR,sizeR,3)
!       complex(8) vLeft(sizeL,sizeL,rLeft,3),vRight(sizeR,sizeR,rRight,3)
!       real(8) Ham(rLeft+rRight+1,rLeft+1,3)
!       ! Internal matrices
!       complex(8) temp(sizeR,sizeL),A(sizeR,sizeL,spin),Aop(sizeR,sizeL,spin),out(sizeR,sizeL,spin),identity(sizeR,sizeL)
!       integer axis,x,y,j,s
!       ! BLAS Level 3 function to compute matrix multiplication
!       external ZGEMM
! 
!       ! Initialize internal matrices
!       A=Reshape(input,(/sizeR,sizeL,spin/))
!       out=0.0d0
!       Left=0.0d0
!       Right=0.0d0
!       do x=1,sizeR
!         Right(x,x)=1.0d0
!       enddo
!       do x=1,sizeL
!         Left(x,x)=1.0d0
!       enddo
! 
!       !First add up terms that do not have an operator at the site
!        do s=1,spin
!           do axis=1,3
! 
!             !Multiply L * (I*A) * hR[axis]
!             !First do hR' *  A' and store in temp
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,hRight(:,:,axis),sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
!             !Now add to out computing out=temp*L'+out
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,Left,sizeL,one,out(:,:,s),sizeR)
! 
!             !Multiply L * (I*A) * DR[axis]
!             !First do DR' *  A' and store in temp
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,DRight(:,:,axis),sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
!             !Now add to out computing out=temp*L'+out
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,Left,sizeL,one,out(:,:,s),sizeR)
! 
!             !Multiply hL[axis] * (I*A) * Right
!             !First do R' *  A' and store in temp
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,Right,sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
!             !Now add to out computing out=temp*hR'+out
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,hLeft(:,:,axis),sizeL,one,out(:,:,s),sizeR)
! 
!             !Multiply DL[axis] * (I*A) * Right
!             !First do R' *  A' and store in temp
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,Right,sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
!             !Now add to out computing out=temp*DL'+out
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,DLeft(:,:,axis),sizeL,one,out(:,:,s),sizeR)
! 
!             do x=1,rLeft
!                do y=1,rRight
!                   !Check for distance, has to be less than maximum 
!                   if ((y+x).lt.max(rLeft,rRight)) then
!                     !Multiply vLeft[axis,0-x] * (I*A) * vRight[axis,0+y] * Ham[rL+1-x,rL+1+y]
!                     !First do vRight' *  A' and store in temp
!                     call ZGEMM(normal,normal,sizeR,sizeL,sizeR,Ham(rLeft+1+y,rLeft+1-x,axis),vRight(:,:,axis,y),sizeR,A(:,:,s),sizeR,zero,temp,sizeR)
!                     !Now add to out computing out=temp*vL'+out
!                     call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,vLeft(:,:,axis,x),sizeL,one,out(:,:,s),sizeR)
!                   endif
!                enddo
!             enddo
! 
!          enddo
!       enddo
! 
!       !Now add up terms that DO have an operator at the site
! 	!Define Pauli matrices
! 	sigma=0.0d0
! 	sigma(1,2,1)=0.5d0
! 	sigma(2,1,1)=0.5d0
! 	sigma(1,2,2)=-0.5d0*II
! 	sigma(2,1,2)=0.5d0*II
! 	sigma(1,1,3)=0.5d0----------
! 	sigma(2,2,3)=-0.5d0d0
! 	do axis=1,3
!          ! Multiply matrix by operator
! 	 Aop=0.0d0
!          do s=1,spin
!             temp=0.0d0
!             do j=1,spin
!                temp=temp+sigma(s,j,axis)*A(:,:,j)
!             enddo
!             Aop(:,:,s)=temp
!          enddo
!          do s=1,spin
!             !Multiply L * (Sigma[axis]*A) * R * Ham(site,site)
!             !First do R' *  SA' and store in temp
!             call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,Right,sizeR,Aop(:,:,s),sizeR,zero,temp,sizeR)
! 	    !Now add to out computing out=temp*L'+out
! 	    temp=temp*Ham(rLeft+1,rLeft+1,axis)
! 	    call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,Left,sizeL,one,out(:,:,s),sizeR)
!             
!             do x=1,rLeft
!                !Multiply vLeft[axis,rL+1-x] * (Sigma[axis]*A) * Right * Ham(rL+1,rL+1-x)
!                !First do R' *  A' and store in temp
! 		call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,Right,sizeR,Aop(:,:,s),sizeR,zero,temp,sizeR)
!                !Now add to out computing out=temp*vL'+out
! 		temp=temp*Ham(rLeft+1,rLeft+1-x,axis)
!                call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,vLeft(:,:,axis,x),sizeL,one,out(:,:,s),sizeR)
!             enddo
! 
!             do x=1,rRight   
!                !Multiply vLeft[axis,0-x] * (Sigma[axis]*A) * Right 
!                !First do vRight' *  A' and store in temp
! 		call ZGEMM(normal,normal,sizeR,sizeL,sizeR,one,vRight(:,:,axis,x),sizeR,Aop(:,:,s),sizeR,zero,temp,sizeR)
!                !Now add to out computing out=temp*vL'+out
! 		temp=temp*Ham(rLeft+1+x,rLeft+1,axis)
!                call ZGEMM(normal,normal,sizeR,sizeL,sizeL,one,temp,sizeR,Left,sizeL,one,out(:,:,s),sizeR)
!             enddo
! 
!           enddo
!        enddo
! 				  
!       output=Reshape(out,(/Lv/))
!       
!     end subroutine MatrixTimesVectorGOOD
! 
