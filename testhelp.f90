
  module TestHelper
       integer,save :: testsOK=0,testsNOT=0
       logical,parameter :: PASSED=.true.,FAILED=.false.
    contains
      subroutine PrintTestMsg(message, status)
        character*(*) message
        logical status
        if (status) then
           print *,message,' = PASSED '
           testsOK=testsOK+1
        else
           print *,message,' = FAILED '
           testsNOT=testsNOT+1
        endif
      end subroutine PrintTestMsg
    end module TestHelper
