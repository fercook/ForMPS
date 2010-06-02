test_suite Operator_Class

setup
  !Set testing mode
  MaxErrorAllowed=CriticalError
end setup

teardown

end teardown

test ApplyOperator
   type(MPSTensor) :: A,B,C
   integer,parameter :: DleftT=4, DrightT=3,spinT=2
   complex(8) :: data(DleftT,DrightT,spinT)
   complex(8) :: Correct(DleftT,DrightT,spinT)
   integer error

   !Initialization
   data(:,:,1)=one*Reshape( [-0.120802, -0.736476, 0.639716, -0.0536866, -0.304505, -0.320898, &
   			   & -0.325939, -0.200223, -0.488209, 0.0946794, 0.17894, 0.76573], [DleftT,DrightT])
   data(:,:,2)=one*Reshape( [-0.258579, -0.424792, -0.512473, -0.130085, -0.442283, -0.00921493, &
   			   & -0.298086, 0.183901, -0.625986, 0.406363, 0.317842, -0.565636], [DleftT,DrightT])
   Correct(:,:,1)=data(:,:,1)
   Correct(:,:,2)=(-1.0d0)*data(:,:,2)

     A=new_MPSTensor(SpinT,DleftT,DrightT,data)
     B=new_MPSTensor(SpinT,DleftT,DrightT,Correct)
     C=PauliSigma(3).ApplyTo.(A)
    assert_equal_within(B.absdiff.C, 0.0d0, 1.0e-5)
     assert_false(WasThereError())
end test

test Site_and_Bond_terms
     type(SiteTerm) :: S
     type(SiteTermList) :: SList
     integer :: site=5
     complex(8) :: amp=1.0d0
     S=new_SiteTerm(PauliSigma(3),site,amp)
     SList=new_SiteTermList(PauliSigma(3),site,amp)
     !call ReduceAmplitude(SList)
     assert_false(WasThereError())
end test

end test_suite

