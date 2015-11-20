module papi_module
!    implicit none
   #define PAPI 
    real :: rtime, ptime, mflops
    integer*8 :: flpops
    
    integer*8 :: rpn2calls = 0, &
                 rpt2calls = 0, &
                 step2calls= 0

    integer :: check
    logical :: flopscounted = .false.
    double precision :: avg_flpops = 0.d0
    save
!contains
end module papi_module
