!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of ForMPS
!
!    ForMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ForMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ForMPS.  If not, see <http://www.gnu.org/licenses/>.

test_suite PEPSAlgorithms_Class

setup
  !Set testing mode
  MaxErrorAllowed=Warning
end setup

teardown

end teardown


test crossed_OVerlaps
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapCN,overlapCC,overlapNC,overlapNN,overlapNormal
    real(8) :: someNorm1,someNorm2
    integer :: x,y
    integer,parameter :: XcanonPos=2,YcanonPos=2, Xlength=4,Ylength=4
    integer,parameter :: MaxLongitudinalBond=6,MaxTransverseBond=6

    onePEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(onePEPS)
    smallPEPS1=onePEPS
    twoPEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(twoPEPS)
    smallPEPS2=twoPEPS

   call smallPEPS1%CanonizeAt(XcanonPos,YcanonPos,VERTICAL,6,6,someNorm1)
   call smallPEPS2%CanonizeAt(XcanonPos,YcanonPos,VERTICAL,6,6,someNorm2)

   overlapNN=Overlap(onePEPS,twoPEPS)
   overlapCN=Overlap(smallPEPS1,twoPEPS)
   overlapNC=Overlap(onePEPS,smallPEPS2)
   overlapCC=Overlap(smallPEPS1,smallPEPS2)

   print *,'Overlaps:'
   print *,' - NN = ',overlapNN
   print *,' - CN = ',overlapCN
   print *,' - NC = ',overlapNC
   print *,' - CC = ',overlapCC
end test



test Horiz_Vertical_Canon

   type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
   type(PEPSTensor) :: aTensor
   type(Tensor2) :: aMatrix
   integer :: l,r,u,d,s,x,y
   complex(8) :: overlapCanon,overlapNormal
   complex(8),allocatable :: anArray(:,:,:,:,:)

   onePEPS=new_PEPS(4,4,2,2)
   twoPEPS=new_PEPS(4,4,2,2)
   allocate (anArray(2,2,2,2,2))
   do s=1,2
   do l=1,2
   do r=1,2
   do u=1,2
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   do x=2,3
    do y=2,3
      call onePEPS%SetTensorAt(x,y,aTensor)
    enddo
   enddo
   aTensor=0.5*ONE*aTensor
   do x=2,3
    do y=2,3
      call twoPEPS%SetTensorAt(x,y,aTensor)
    enddo
   enddo

   !LEFT COL
   deallocate(anArray)
   allocate (anArray(1,2,2,2,2))
   do s=1,2
   do l=1,1
   do r=1,2
   do u=1,2
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   do y=2,3
     call onePEPS%SetTensorAt(1,y,aTensor)
   enddo
   aTensor=0.25*ONE*aTensor
   do y=2,3
     call twoPEPS%SetTensorAt(1,y,aTensor)
   enddo

!RIGHT COL
   deallocate(anArray)
   allocate (anArray(2,1,2,2,2))
   do s=1,2
   do l=1,2
   do r=1,1
   do u=1,2
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   do y=2,3
     call onePEPS%SetTensorAt(4,y,aTensor)
   enddo
   aTensor=0.25*ONE*aTensor
   do y=2,3
     call twoPEPS%SetTensorAt(4,y,aTensor)
   enddo

!LOWER ROW
   deallocate(anArray)
   allocate (anArray(2,2,2,1,2))
   do s=1,2
   do l=1,2
   do r=1,2
   do u=1,2
   do d=1,1
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   do x=2,3
     call onePEPS%SetTensorAt(x,1,aTensor)
   enddo
   aTensor=0.25*ONE*aTensor
   do x=2,3
     call twoPEPS%SetTensorAt(x,1,aTensor)
   enddo

!UPPER ROW
   deallocate(anArray)
   allocate (anArray(2,2,1,2,2))
   do s=1,2
   do l=1,2
   do r=1,2
   do u=1,1
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   do x=2,3
     call onePEPS%SetTensorAt(x,4,aTensor)
   enddo
   aTensor=0.25*ONE*aTensor
   do x=2,3
     call twoPEPS%SetTensorAt(x,4,aTensor)
   enddo

   !1,1
   deallocate(anArray)
   allocate (anArray(1,2,2,1,2))
   do s=1,2
   do l=1,1
   do r=1,2
   do u=1,2
   do d=1,1
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   call onePEPS%SetTensorAt(1,1,aTensor)
   aTensor=3.0*ONE*aTensor
   call twoPEPS%SetTensorAt(1,1,aTensor)

   !4,1
   deallocate(anArray)
   allocate (anArray(2,1,2,1,2))
   do s=1,2
   do l=1,2
   do r=1,1
   do u=1,2
   do d=1,1
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   call onePEPS%SetTensorAt(4,1,aTensor)
   aTensor=3.0*ONE*aTensor
   call twoPEPS%SetTensorAt(4,1,aTensor)

   !1,4
   deallocate(anArray)
   allocate (anArray(1,2,1,2,2))
   do s=1,2
   do l=1,1
   do r=1,2
   do u=1,1
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   call onePEPS%SetTensorAt(1,4,aTensor)
   aTensor=3.0*ONE*aTensor
   call twoPEPS%SetTensorAt(1,4,aTensor)

   !4,4
   deallocate(anArray)
   allocate (anArray(2,1,1,2,2))
   do s=1,2
   do l=1,2
   do r=1,1
   do u=1,1
   do d=1,2
      anArray(l,r,u,d,s)=10*(s-1)+l+r+u+d-4
   enddo
   enddo
   enddo
   enddo
   enddo
   aTensor=new_PEPSTensor(new_Tensor(AnArray))
   call onePEPS%SetTensorAt(4,4,aTensor)
   aTensor=3.0*ONE*aTensor
   call twoPEPS%SetTensorAt(4,4,aTensor)


   smallPEPS1=onePEPS
   smallPEPS2=onePEPS


   do y=1,4
      aTensor=onePEPS%GetTensorAt(1,y)
      aMatrix=aTensor%CollapseAllIndicesBut(RIGHT)
      aTensor=onePEPS%GetTensorAt(y,4)
      assert_equal_within(aMatrix.absdiff.aTensor%CollapseAllIndicesBut(DOWN),0.0d0,1.0d-10)
   enddo

!   overlapNormal=Overlap(onePEPS,twoPEPS)
!   overlapCanon=Overlap(smallPEPS1,smallPEPS2)
!   print *,'Initial overlaps:',Abs(overlapNormal),Abs(overlapCanon)

   call DecouplePEPSVerticallyAtCol(smallPEPS1,1,RIGHT,6,6)
   call DecouplePEPSHorizontallyAtRow(smallPEPS2,4,DOWN,6,6)

!   do y=1,4
!      print *,'Tensor number:',y
!      aTensor=smallPEPS1%GetTensorAt(1,y)
!      call aTensor%Print('Horizontal PEPS:')
!      aTensor=smallPEPS2%GetTensorAt(y,4)
!      call aTensor%Print('Vertical PEPS:')
!   enddo


end test

test gradual_Canonification_Horiz
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapCanon,overlapNormal
    real(8) :: someNorm1,someNorm2
    integer :: x,y
    integer,parameter :: XcanonPos=2,YcanonPos=2, Xlength=4,Ylength=4
    integer,parameter :: MaxLongitudinalBond=6,MaxTransverseBond=6

    onePEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(onePEPS)
    smallPEPS1=onePEPS
    twoPEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(twoPEPS)
    smallPEPS2=twoPEPS

    overlapNormal=Overlap(onePEPS,twoPEPS)
    overlapCanon=Overlap(smallPEPS1,smallPEPS2)
    print *,'Initial overlaps:',Abs(overlapNormal),Abs(overlapCanon)

    do x=XLength,XcanonPos+1,-1
      if (debug) print *,'Canonizing PEPS to down, col:',x
      call DecouplePEPSVerticallyAtCol(smallPEPS1,x,LEFT,MaxLongitudinalBond,MaxTransverseBond)
      call DecouplePEPSVerticallyAtCol(smallPEPS2,x,LEFT,MaxLongitudinalBond,MaxTransverseBond)
      overlapCanon=Overlap(smallPEPS1,smallPEPS2)
      print *,'---- overlap at x:',Abs(overlapNormal),Abs(overlapCanon)
      assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
    enddo

    do x=1,XcanonPos-1
      if (debug) print *,'Canonizing PEPS to up, col:',x
      call DecouplePEPSVerticallyAtCol(smallPEPS1,x,RIGHT,MaxLongitudinalBond,MaxTransverseBond)
      call DecouplePEPSVerticallyAtCol(smallPEPS2,x,RIGHT,MaxLongitudinalBond,MaxTransverseBond)
      overlapCanon=Overlap(smallPEPS1,smallPEPS2)
      print *,'---- overlap at x:',Abs(overlapNormal),Abs(overlapCanon)
      assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
    enddo

         do y=1,YcanonPos-1
            if (debug) print *,'Canonizing PEPS up, row:',y
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,XcanonPos,y,UP,MaxTransverseBond**2)
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,XcanonPos,y,UP,MaxTransverseBond**2)
            overlapCanon=Overlap(smallPEPS1,smallPEPS2)
            print *,'---- overlap at y:',Abs(overlapNormal),Abs(overlapCanon)
            assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
         enddo
         do y=YLength,YcanonPos+1,-1
            if (debug) print *,'Canonizing PEPS to left, col:',x
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,XcanonPos,y,DOWN,MaxTransverseBond**2)
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,XcanonPos,y,DOWN,MaxTransverseBond**2)
            overlapCanon=Overlap(smallPEPS1,smallPEPS2)
            print *,'---- overlap at y:',Abs(overlapNormal),Abs(overlapCanon)
            assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
         enddo

end test



test gradual_Canonification_Vertical
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapCanon,overlapNormal
    real(8) :: someNorm1,someNorm2
    integer :: x,y
    integer,parameter :: XcanonPos=2,YcanonPos=2, Xlength=4,Ylength=4
    integer,parameter :: MaxLongitudinalBond=6,MaxTransverseBond=6

    onePEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(onePEPS)
    smallPEPS1=onePEPS
    twoPEPS=new_PEPS(Xlength,Ylength,2,2)
    call Normalize(twoPEPS)
    smallPEPS2=twoPEPS

    overlapNormal=Overlap(onePEPS,twoPEPS)
    overlapCanon=Overlap(smallPEPS1,smallPEPS2)

    do y=1,YcanonPos-1
      call DecouplePEPSHorizontallyAtRow(smallPEPS1,y,UP,MaxLongitudinalBond,MaxTransverseBond)
      call DecouplePEPSHorizontallyAtRow(smallPEPS2,y,UP,MaxLongitudinalBond,MaxTransverseBond)
      overlapCanon=Overlap(smallPEPS1,smallPEPS2)
      assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
    enddo
    do y=YLength,YcanonPos+1,-1
      call DecouplePEPSHorizontallyAtRow(smallPEPS1,y,DOWN,MaxLongitudinalBond,MaxTransverseBond)
      call DecouplePEPSHorizontallyAtRow(smallPEPS2,y,DOWN,MaxLongitudinalBond,MaxTransverseBond)
      overlapCanon=Overlap(smallPEPS1,smallPEPS2)
      assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
    enddo

         do x=1,XcanonPos-1
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,x,YcanonPos,RIGHT,MaxTransverseBond**2)
            call CanonizeCoreSiteAndPushNorm(smallPEPS2,x,YcanonPos,RIGHT,MaxTransverseBond**2)
            overlapCanon=Overlap(smallPEPS1,smallPEPS2)
            assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
         enddo
         do x=XLength,XcanonPos+1,-1
            call CanonizeCoreSiteAndPushNorm(smallPEPS1,x,YcanonPos,LEFT,MaxTransverseBond**2)
            call CanonizeCoreSiteAndPushNorm(smallPEPS2,x,YcanonPos,LEFT,MaxTransverseBond**2)
            overlapCanon=Overlap(smallPEPS1,smallPEPS2)
            assert_equal_within(Abs(overlapNormal)-Abs(overlapCanon),0.0d0,1.0d-10)
         enddo

end test





test Canonical_Overlap_VerticalCanon
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapBS,overlapNormal
    real(8) :: someNorm1,someNorm2
    integer,parameter :: XcanonPos=2,YcanonPos=3

    onePEPS=new_PEPS(4,4,2,2)
    call Normalize(onePEPS)

    smallPEPS1=onePEPS
    call smallPEPS1%CanonizeAt(XcanonPos,YcanonPos,VERTICAL,6,6,someNorm1)
    overlapBS=Overlap(smallPEPS1,onePEPS)
    print *, 'overlap BTW Canonized PEPS and full PEPS A: ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    overlapBS=Overlap(smallPEPS1,smallPEPS1,CorePosition=YcanonPos, CoreDirection=HORIZONTAL)
    print *, 'overlap BTW Canonized PEPSs (CANON METHOD): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)
    overlapBS=Overlap(smallPEPS1,smallPEPS1)
    print *, 'overlap BTW Canonized PEPSs (EXACT): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    twoPEPS=new_PEPS(4,4,2,2)
    call Normalize(twoPEPS)
    smallPEPS2=twoPEPS
    call smallPEPS2%CanonizeAt(XcanonPos,YcanonPos,VERTICAL,6,6,someNorm2)
    overlapBS=Overlap(smallPEPS2,twoPEPS)
    print *, 'overlap BTW Canonized PEPS and full PEPS B: ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    overlapNormal=Overlap(onePEPS,twoPEPS)
    overlapBS=Overlap(smallPEPS1,smallPEPS2,CorePosition=YcanonPos, CoreDirection=HORIZONTAL)
    print *,abs(overlapNormal),abs(overlapBS)
    print *, 'difference in overlap calculated normally and canonically: ',abs(overlapBS)-abs(overlapNormal)
    assert_equal_within(abs(overlapBS)-abs(overlapNormal),0.0d0,1.0d-10)

end test

test Canonical_Overlap_HorizontalCanon
    type(PEPS) :: onePEPS,twoPEPS,smallPEPS1,smallPEPS2
    complex(8) :: overlapBS,overlapNormal
    real(8) :: someNorm1,someNorm2
    integer,parameter :: XcanonPos=2,YcanonPos=3

    onePEPS=new_PEPS(4,4,2,2)
    call Normalize(onePEPS)

    smallPEPS1=onePEPS
    call smallPEPS1%CanonizeAt(XcanonPos,YcanonPos,HORIZONTAL,6,6,someNorm1)
    overlapBS=Overlap(smallPEPS1,onePEPS)
    print *, 'overlap BTW Canonized PEPS and full PEPS A: ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    overlapBS=Overlap(smallPEPS1,CorePosition=XcanonPos, CoreDirection=VERTICAL)
    print *, 'overlap BTW Canonized PEPSs (CANON METHOD): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)
    overlapBS=Overlap(smallPEPS1,smallPEPS1)
    print *, 'overlap BTW Canonized PEPSs (EXACT): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    twoPEPS=new_PEPS(4,4,2,2)
    call Normalize(twoPEPS)
    smallPEPS2=twoPEPS
    call smallPEPS2%CanonizeAt(XcanonPos,YcanonPos,HORIZONTAL,6,6,someNorm2)
    overlapBS=Overlap(smallPEPS2,twoPEPS)
    print *, 'overlap BTW Canonized PEPS and full PEPS B: ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)

    overlapNormal=Overlap(onePEPS,twoPEPS)
    overlapBS=Overlap(smallPEPS1,smallPEPS2,CorePosition=XcanonPos, CoreDirection=VERTICAL)
    print *,abs(overlapNormal),abs(overlapBS)
    print *, 'difference in overlap calculated normally and canonically: ',abs(overlapBS)-abs(overlapNormal)
    assert_equal_within(abs(overlapBS)-abs(overlapNormal),0.0d0,1.0d-10)


end test


test PEPS_Canonization_Routines
    type(PEPS) :: aPEPS,smallPEPS
    type(PEPSTensor) :: aTensor
    integer :: dims(4)
    complex(8) :: overlapBS
    real(8) :: theNorm

    aPEPS=new_PEPS(4,4,2,2)
    call Normalize(aPEPS)
    print *,'Norm of big overlap is',Overlap_PEPS(aPEPS,aPEPS)
    smallPEPS=aPEPS
    call smallPEPS%CanonizeAt(2,2,VERTICAL,8,8,theNorm)
    call smallPEPS%PrintbondDimensions()
    assert_true(smallPEPS%IsPEPSWellFormed())
    overlapBS=Overlap_PEPS(smallPEPS,smallPEPS)
    print *,'Norm of canonical state is',overlapBS,theNorm
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-10)
    overlapBS=Overlap_PEPS(smallPEPS,aPEPS)
    print *, 'overlap BTW peps and Canonized PEPS (EXACT): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-9)

    aPEPS=new_PEPS(4,4,2,3)
    call Normalize(aPEPS)
    smallPEPS=aPEPS
    print *,'About to canonize...'
    call smallPEPS%CanonizeAt(3,2,VERTICAL,6,6)
    overlapBS=Overlap_PEPS(aPEPS,smallPEPS)
    print *, 'overlap BTW peps and Canonized PEPS (H small bond): ',abs(overlapBS)**2
    assert_equal_within(abs(overlapBS)**2,1.0d0,1.0d-5)

end test


test OverlapAlgorithm
  type(PEPS) :: aPEPS
  integer :: length=4,width=4,spin=2,bond=2, error
  complex(8) :: overlap12

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  print *,'ABOUT TO COMPUTE OVERLAP'
  overlap12 = Overlap_PEPS(aPEPS)
  print *,'PREVIOUS OVERLAP',overlap12
  assert_false(abs(overlap12)**2.eq.1.0d0)

  print *,'About to NORMALIZE ----------'
  call Normalize(aPEPS)
  overlap12 = Overlap_PEPS(aPEPS,aPEPS)
  print *,'OVERLAP AFTER NORMALIZATION',overlap12
  assert_equal_within(abs(overlap12)**2,1.0d0,1.0d-8)

  assert_false(WasThereError())

  error= aPEPS%Delete()

end test

test Progressive_truncation
  type(PEPS) :: aPEPS,smallPEPS
  integer :: length=4,width=4,spin=2,bond=3, error
  complex(8) :: overlap12
  integer :: smallbond

  aPEPS=new_PEPS(length,width,spin,bond)
  call aPEPS%ScaleBy(ONE/(4.0d0)**(1.0d0/2.0d0))
  call Normalize(aPEPS)
  do smallbond=bond-1,1,-1
      print *,'Reducing Bond ',bond,' to ',smallbond
      smallPEPS=ReduceMAXPEPSBond(aPEPS,smallbond)
      call Normalize(smallPEPS)
      overlap12 = Overlap_PEPS(aPEPS,smallPEPS)
      print *,'Bond: ',smallbond,', overlap: ',abs(overlap12)**2
      assert_true(abs(overlap12)**2.gt.0.9)
  enddo

  assert_false(WasThereError())

end test




test ExptValueProdHamilt
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.d0
    real(8),parameter :: theta=0.5d0
    integer,parameter :: Xsize=4,Ysize=4,row=2
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=1,OperatorDim=2,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,integerONE,BondDim),AR(OperatorDim,BondDim,integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart(OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    !Prepare a product PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=pauliZ
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,cos(theta)**16,1.0d-3)

end test


test ExptValueLocalSumHamtn
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=1.0d0
    real(8),parameter :: theta=4.0d0
    integer,parameter :: Xsize=4,Ysize=4,row=2
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=2,OperatorDim=2,SpinDim=2
    integer :: n,m,k,i,j,p,dims(8)
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim,BondDim,BondDim),AL(OperatorDim,integerONE,BondDim),AR(OperatorDim,BondDim,integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart(OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ

    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,2,2)=ONE
        A(2,1,2)=ONE;
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;
        AR(1,2,1)=ONE;      AR(2,1,1)=ONE;

    !Prepare a product PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=2,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    deallocate(localMatrix)
    allocate(localMatrix(2,2,integerONE,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,1,m,1,1)=localMatrix(i,j,1,m,1,1)+AL(k,1,m)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(1,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,integerONE,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,m,1,1,1)=localMatrix(i,j,m,1,1,1)+AR(k,m,1)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(4,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(2,2,BondDim,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do m=1,BondDim
    do p=1,BondDim
     do i=1,SpinDim
     do j=1,SpinDim
      do k=1,OperatorDim
        localMatrix(i,j,m,p,1,1)=localMatrix(i,j,m,p,1,1)+A(k,m,p)*spinPart(k,i,j)
      enddo
     enddo
     enddo
    enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(2,1,localTensor)
    call theH%SetTensorAt(3,1,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,4*cos(theta),1.0d-3)

end test




test ExptValueHamiltonian
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.17d0
    real(8),parameter :: theta=0.71d0
    integer,parameter :: Xsize=4,Ysize=4,row=1
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim, BondDim, BondDim),AL( OperatorDim, integerONE, BondDim),AR( OperatorDim, BondDim, integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart( OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=ONE;  AR(3,1,1)=field;
    !Prepare an identity PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now set the second row to have the Ising like Hamiltonian
    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE, BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,1
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+AL(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(1,row,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,BondDim,integerONE,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,1
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+AR(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(Xsize,row,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,BondDim,BondDim,integerONE,integerONE))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,n,m,1,1)=localMatrix(:,:,n,m,1,1)+A(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(2,row,localTensor)
    call theH%SetTensorAt(3,row,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)

    assert_equal_within(energy,3*cos(theta)**2+4*field*sin(theta),1.0d-3)

end test


test ExptValueHamiltonianCol
    type(PEPO) :: theH
    type(PEPS) :: theState
    real(8),parameter :: field=0.17d0
    real(8),parameter :: theta=0.71d0
    integer,parameter :: Xsize=4,Ysize=4,col=1
    type(PEPOTensor) :: localTensor
    type(PEPSTensor) :: localPEPS
    complex(8),allocatable :: localMatrix(:,:,:,:,:,:)
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=3,OperatorDim=3,SpinDim=2
    integer :: n,m,k
    complex(8) :: identity(SpinDim,SpinDim),pauliZ(SpinDim,SpinDim),pauliX(SpinDim,SpinDim)
    complex(8) :: A(OperatorDim, BondDim, BondDim),AL( OperatorDim, integerONE, BondDim),AR( OperatorDim, BondDim, integerONE)
    complex(8) :: AIdentity(1,1,1),spinPart( OperatorDim,2,2)
    real(8) :: energy,overlap12

    identity=ZERO;  identity(1,1)=ONE; identity(2,2)=ONE
    pauliZ=ZERO;   pauliZ(1,1)=ONE;  pauliZ(2,2)=-ONE
    pauliX=ZERO;   pauliX(1,2)=ONE;  pauliX(2,1)=ONE

    spinPart(1,:,:)=identity; spinPart(2,:,:)=PauliZ; spinPart(3,:,:)=PauliX
    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,3,3)=ONE
        A(2,1,2)=ONE;       A(2,2,3)=ONE
        A(3,1,3)=field
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;  AL(3,1,3)=field;
        AR(1,3,1)=ONE;      AR(2,2,1)=ONE;  AR(3,1,1)=field;
    !Prepare an identity PEPO
    theH=new_PEPO(Xsize,Ysize,SpinDim,BondDim)
    allocate(localMatrix(2,2,integerONE,integerONE,integerONE,integerONE))
    localMatrix(:,:,1,1,1,1)=identity
    localTensor=new_PEPOTensor(localMatrix)
    do m=1,Ysize
      do n=1,Xsize
        call theH%SetTensorAt(n,m,localTensor)
      enddo
    enddo

    !Now set the first col to have the Ising like Hamiltonian
    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,BondDim,integerONE))
    localMatrix=ZERO
    do n=1,1
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+AL(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(col,1,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,integerONE,BondDim))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,1
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+AR(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(col,Ysize,localTensor)

    deallocate(localMatrix)
    allocate(localMatrix(SpinDim,SpinDim,integerONE,integerONE,BondDim,BondDim))
    localMatrix=ZERO
    do n=1,BondDim
      do m=1,BondDim
        do k=1,OperatorDim
          localMatrix(:,:,1,1,m,n)=localMatrix(:,:,1,1,m,n)+A(k,n,m)*spinPart(k,:,:)
        enddo
      enddo
    enddo
    localTensor=new_PEPOTensor(localMatrix)
    call theH%SetTensorAt(col,2,localTensor)
    call theH%SetTensorAt(col,3,localTensor)

    !call theH%PrintBondDimensions()

    !Now generate PEPS cos(theta/2) |0> + sin(theta/2) |1> in each site (product state)
    allocate(localState(SpinDim,integerONE,integerONE,integerONE,integerONE))
    localState=ZERO
    localState(1,1,1,1,1)=Cos(theta/2.0d0)
    localState(2,1,1,1,1)=Sin(theta/2.0d0)
    theState=New_PEPS(Xsize,Ysize,SpinDim,integerONE)
    localPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=1,Xsize
      do m=1,Xsize
        call theState%SetTensorAt(n,m,localPEPS)
      enddo
    enddo

    overlap12 = Overlap_PEPS(theState,theState)
    assert_equal_within(overlap12,1.0d0,1.0d-8)

    !Now compute energy
    energy=ExpectationValue(theState,theH)
    assert_equal_within(energy,3*cos(theta)**2+4*field*sin(theta),1.0d-3)

end test


!!!!!!!!! THIS IS A VERY UGLY TEST, it just needs to be done by hand or spend some time
!!!!!  refactoring...not for much

test ControlledOverlap

    type(PEPS) :: onePEPS,twoPEPS
    complex(8) :: overlapSingle,overlapDouble
    complex(8),allocatable :: localState(:,:,:,:,:)
    integer,parameter :: BondDim=2,SpinDim=2,XSIZE=4,YSIZE=4
    integer :: l,r,u,d,s,n,m
    complex(8) :: A(SpinDim,BondDim,BondDim),AL(SpinDim,integerONE,BondDim),AR(SpinDim,BondDim,integerONE)
    type(PEPSTensor) :: tempPEPS
    complex(8) :: AIdentity(1,1,1)
    real(8) :: overlap12

    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,2,2)=0.5d0*II
        A(2,1,2)=ONE;
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;
        AR(1,2,1)=ONE;      AR(2,1,1)=ONE;

    onePEPS=new_PEPS(Xsize,Ysize,SpinDim,BondDim)
    allocate(localState(2,BondDim,BondDim,BondDim,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do m=2,Ysize-1
      do n=2,Xsize-1
        call onePEPS%SetTensorAt(n,m,tempPEPS)
      enddo
    enddo

    deallocate(localState)
    allocate(localState(2,integerONE,BondDim,BondDim,BondDim))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Ysize-1
        call onePEPS%SetTensorAt(1,n,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,integerONE,BondDim,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Ysize-1
        call onePEPS%SetTensorAt(Xsize,n,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,BondDim,integerONE,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Xsize-1
        call onePEPS%SetTensorAt(n,Ysize,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,BondDim,BondDim,integerONE))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Xsize-1
        call onePEPS%SetTensorAt(n,1,tempPEPS)
    enddo


   deallocate(localState)
    allocate(localState(2,integerONE,BondDim,BondDim,integerONE))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call onePEPS%SetTensorAt(1,1,tempPEPS)

   deallocate(localState)
    allocate(localState(2,integerONE,BondDim,integerONE,BondDim))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call onePEPS%SetTensorAt(1,Ysize,tempPEPS)

   deallocate(localState)
    allocate(localState(2,BondDim,integerONE,integerONE,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call onePEPS%SetTensorAt(Xsize,Ysize,tempPEPS)


   deallocate(localState)
    allocate(localState(2,BondDim,integerONE,BondDim,integerONE))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call onePEPS%SetTensorAt(Xsize,1,tempPEPS)

    overlap12 = Overlap_PEPS(onePEPS,onePEPS)
    print *,overlap12
    assert_equal_within(overlap12,0,005859375.0d0,1.0d-12)



    A=ZERO    !BondDim = 2
        A(1,1,1)=ONE;       A(1,2,2)=0.3d0*II
        A(2,1,2)=ONE;
    AL=ZERO;   AR=ZERO;
        AL(1,1,1)=ONE;      AL(2,1,2)=ONE;
        AR(1,2,1)=ONE;      AR(2,1,1)=ONE;

   deallocate (localState)
    twoPEPS=new_PEPS(Xsize,Ysize,SpinDim,BondDim)
    allocate(localState(2,BondDim,BondDim,BondDim,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do m=2,Ysize-1
      do n=2,Xsize-1
        call twoPEPS%SetTensorAt(n,m,tempPEPS)
      enddo
    enddo

    deallocate(localState)
    allocate(localState(2,integerONE,BondDim,BondDim,BondDim))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Ysize-1
        call twoPEPS%SetTensorAt(1,n,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,integerONE,BondDim,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,BondDim
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*A(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Ysize-1
        call twoPEPS%SetTensorAt(Xsize,n,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,BondDim,integerONE,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Xsize-1
        call twoPEPS%SetTensorAt(n,Ysize,tempPEPS)
    enddo

    deallocate(localState)
    allocate(localState(2,BondDim,BondDim,BondDim,integerONE))
    localState=ZERO
    do l=1,BondDim
     do r=1,BondDim
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=A(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    do n=2,Xsize-1
        call twoPEPS%SetTensorAt(n,1,tempPEPS)
    enddo


   deallocate(localState)
    allocate(localState(2,integerONE,BondDim,BondDim,integerONE))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call twoPEPS%SetTensorAt(1,1,tempPEPS)

   deallocate(localState)
    allocate(localState(2,integerONE,BondDim,integerONE,BondDim))
    localState=ZERO
    do l=1,1
     do r=1,BondDim
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AL(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call twoPEPS%SetTensorAt(1,Ysize,tempPEPS)

   deallocate(localState)
    allocate(localState(2,BondDim,integerONE,integerONE,BondDim))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,1
      do d=1,BondDim
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*AL(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call twoPEPS%SetTensorAt(Xsize,Ysize,tempPEPS)


   deallocate(localState)
    allocate(localState(2,BondDim,integerONE,BondDim,integerONE))
    localState=ZERO
    do l=1,BondDim
     do r=1,1
     do u=1,BondDim
      do d=1,1
      do s=1,SpinDim
        localState(s,l,r,u,d)=AR(s,l,r)*AR(s,u,d)
      enddo
      enddo
     enddo
     enddo
    enddo
    tempPEPS=new_PEPSTensor(localState(1,:,:,:,:),localState(2,:,:,:,:))
    call twoPEPS%SetTensorAt(Xsize,1,tempPEPS)

    overlap12 = Overlap_PEPS(onePEPS,twoPEPS)
    print *,overlap12
    assert_equal_within(overlap12,0,000273375d0,1.0d-12)

end test




end test_suite
