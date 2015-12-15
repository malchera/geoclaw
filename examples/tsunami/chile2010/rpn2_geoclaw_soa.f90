!======================================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx, &
        ql_aos,qr_aos,auxl_aos,auxr_aos,fwave,s,amdq_aos,&
        apdq_aos)
!======================================================================
!
! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!     with topography:
!     #        h_t + (hu)_x + (hv)_y = 0                           #
!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.
!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
!  From the basic clawpack routines, this routine is called with
!     ql = qr
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!       André Malcher: Changed to FORTRAN 90                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use geoclaw_module, only: g => grav, drytol => dry_tolerance
    use geoclaw_module, only: earth_radius, deg2rad
    use amr_module, only: mcapa

    implicit none

    !temporary arrays to do AoS->SoA tranform
    double precision, intent(in) :: ql_aos(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in) :: qr_aos(meqn, 1-mbc:maxm+mbc)
    double precision, intent(in) :: auxl_aos(maux,1-mbc:maxm+mbc)
    double precision, intent(in) :: auxr_aos(maux,1-mbc:maxm+mbc)
    
    double precision, intent(out) :: apdq_aos(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: amdq_aos(meqn,1-mbc:maxm+mbc)

    !input
    integer :: maxm,meqn,maux,mwaves,mbc,mx,ixy
    
    !INPUT
    ! TODO: Question: Is meqn,mwaves,maux ever != 3? Inconsistent in code
    ! If not, some optimizations possible
    double precision :: hL(1-mbc:maxm+mbc-1), hR(1-mbc:maxm+mbc-1)
    double precision :: huL(1-mbc:maxm+mbc-1), huR(1-mbc:maxm+mbc-1)
    double precision :: hvL(1-mbc:maxm+mbc-1), hvR(1-mbc:maxm+mbc-1)
    double precision :: uL(1-mbc:maxm+mbc-1), uR(1-mbc:maxm+mbc-1)
    double precision :: vL(1-mbc:maxm+mbc-1), vR(1-mbc:maxm+mbc-1)
    double precision :: phiL(1-mbc:maxm+mbc-1), phiR(1-mbc:maxm+mbc-1)
    double precision :: bL(1-mbc:maxm+mbc-1), bR(1-mbc:maxm+mbc-1)
    double precision :: sE1(1-mbc:maxm+mbc-1), sE2(1-mbc:maxm+mbc-1)

    !OUTPUT
    ! TODO: Change order of s and fwave arrays -> Both are output
    ! arrays, i.e. we only need to copy them in the end 
    double precision, intent(out) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    double precision, intent(out) :: s(mwaves, 1-mbc:maxm+mbc)
    double precision :: apdq(1-mbc:maxm+mbc,meqn)
    double precision :: amdq(1-mbc:maxm+mbc,meqn)

    !local only
    integer :: m,i,mw,maxiter,mu,nv
    double precision :: wall(3)
    double precision :: fw(3,3)
    double precision :: sw(3)

    ! TODO: Check whether operations on these variables can be vectorized at all
    ! if they're not declared as arrays!!! => Maybe work on original q? Shouldn't
    ! have any consequences as long as q is copied like here...
    double precision :: sqrt_ghL, sqrt_ghR ! Temporary variables for sqrt(g*h)
    double precision :: sL,sR,sRoe1,sRoe2,uhat,chat
    double precision :: s1m,s2m
    double precision :: hstar,hstartest,hstarHLL,sLtest,sRtest
    double precision :: tw,dxdc

    logical :: rare1,rare2
    
    !dir$ assume_aligned fwave:64
    !dir$ assume_aligned s:64
    !dir$ assume_aligned apdq:64
    !dir$ assume_aligned amdq:64
    
    ! General:
    ! TODO: Check all SIMD compiler directives for correctness and reasonability
    ! !$OMP SIMD -> This one seems to be extremely unefficient due
    ! to strided accesses and misalignment

    !set normal direction
    if (ixy.eq.1) then
        mu=2
        nv=3
    else
        mu=3
        nv=2
    endif
!!! AoS to SoA SECTION !!!    
    do m = 1,meqn
        do i = 1-mbc,maxm+mbc-1
            hL(i) = ql_aos(1,i)
            hR(i) = ql_aos(1, i+1) 
            huL(i) = ql_aos(mu,i)
            huR(i) = ql_aos(mu,i+1)
            hvL(i) = ql_aos(nv,i)
            hvR(i) = ql_aos(nv,i+1)
        enddo
    enddo

    bL(:) = auxl_aos(1,1-mbc:mx+mbc-1)
    bR(:) = auxl_aos(1,2-mbc:mx+mbc)
!!! AoS to SoA SECTION !!!
    
    fwave = 0.d0 ! TODO: Does this and the next line work as expected? (i.e. sets whole array to 0?)
    s = 0.d0
   
    !TODO: inform of a bad riemann problem from the start
    !zero (small) negative values if they exist => TODO: use built-in "where (q < 0.d0) q = 0.d0"?
    do i=1-mbc,mx+mbc
        if (hL(i) .lt. drytol) then
            ! This case uses <drytol rather than <0.d0, eases next loop
            hL(i) = 0.d0
            huL(i) = 0.d0
            hvL(i) = 0.d0
            uL(i) = 0.d0
            vL(i) = 0.d0
            phiL(i) = 0.d0
        else ! This could also be else part of if() in the previous loop
            ! Need to check whether this would be efficient though (no blending operation!)
            ! TODO: Next six lines really required for masking?
            hL(i) = hL(i)
            huL(i) = huL(i)
            hvL(i) = hvL(i)
            uL(i)   = huL(i) / hL(i)
            vL(i)   = hvL(i) / hL(i)
            phiL(i) = 0.5d0*g*hL(i)*hL(i) + huL(i)*huL(i)/hL(i)
        endif
        if (hR(i) .lt. drytol) then
            ! This case uses <drytol rather than <0.d0, eases next loop
            hR(i) = 0.d0
            huR(i) = 0.d0
            hvR(i) = 0.d0
            uR(i) = 0.d0
            vR(i) = 0.d0
            phiR(i) = 0.d0
        else ! This could also be else part of if() in the previous loop
            ! Need to check whether this would be efficient though (no blending operation!)
            ! TODO: Next six lines really required for masking?
            hR(i) = hR(i)
            huR(i) = huR(i)
            hvR(i) = hvR(i)
            uR(i)   = huR(i) / hR(i)
            vR(i)   = hvR(i) / hR(i)
            phiR(i) = 0.5d0*g*hR(i)*hR(i) + huR(i)*huR(i)/hR(i)
        endif
    enddo

    !loop through Riemann problems at each grid cell
    !$omp simd private(wall, hstar,hstartest, s1m, s2m, rare1, rare2, &
    !$omp&      sqrt_ghL, sqrt_ghR, sL, sR, chat, uhat, sRoe1, SRoe2, maxiter, &
    !$omp&      drytol, g, sw, fw)
    do i=1-mbc,mx+mbc-1 ! mx+2mbc-1 Riemann problems, where i corresponds to i-th Riemann problem
        !-----------------------Initializing-----------------------------------

        !skip problem if in a completely dry area
        ! Check: cycle/continue problem for vectorization?
        ! => YES! See https://software.intel.com/en-us/node/524555
        !if (hL(i) <= drytol .and. hR(i) <= drytol) then
        !    cycle
        !endif

        !Riemann problem variables
        wall(1:3) = 1.d0
        
        if (hR(i).le.drytol) then
            !dir$ forceinline
            call riemanntype(hL(i),hL(i),uL(i),-uL(i),hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hL(i),hstar)

            if (hstartest+bL(i).lt.bR(i)) then !right state should become ghost values that mirror left for wall problem
                wall(2)=0.d0
                wall(3)=0.d0
                hR(i)=hL(i)
                huR(i)=-huL(i)
                bR(i)=bL(i)
                phiR(i)=phiL(i)
                uR(i)=-uL(i)
                vR(i)=vL(i)
            elseif (hL(i)+bL(i).lt.bR(i)) then
                bR(i)=hL(i)+bL(i)
            endif
        elseif (hL(i).le.drytol) then ! right surface is lower than left topo
            !dir$ forceinline
            call riemanntype(hR(i),hR(i),-uR(i),uR(i),hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hR(i),hstar)

            if (hstartest+bR(i).lt.bL(i)) then  !left state should become ghost values that mirror right
                wall(1)=0.d0
                wall(2)=0.d0
                hL(i)=hR(i)
                huL(i)=-huR(i)
                bL(i)=bR(i)
                phiL(i)=phiR(i)
                uL(i)=-uR(i)
                vL(i)=vR(i)
            elseif (hR(i)+bR(i).lt.bL(i)) then
                bL(i)=hR(i)+bR(i)
            endif
        endif

        !determine wave speeds
        sqrt_ghL = sqrt(g*hL(i))
        sqrt_ghR = sqrt(g*hR(i))
        sL=uL(i)-sqrt_ghL ! 1 wave speed of left state
        sR=uR(i)+sqrt_ghR ! 2 wave speed of right state

        uhat=(sqrt_ghL*uL(i) + sqrt_ghR*uR(i))/(sqrt_ghR+sqrt_ghL) ! Roe average
        chat=sqrt(g*0.5d0*(hR(i)+hL(i))) ! Roe average
        sRoe1=uhat-chat ! Roe wave speed 1 wave
        sRoe2=uhat+chat ! Roe wave speed 2 wave

        sE1(i) = min(sL,sRoe1) ! Eindfeldt speed 1 wave
        sE2(i) = max(sR,sRoe2) ! Eindfeldt speed 2 wave

        !--------------------end initializing...finally----------
        !solve Riemann problem.

        maxiter = 1

        call riemann_fwave(meqn,mwaves,hL(i),hR(i),huL(i),huR(i),hvL(i),hvR(i), &
            bL(i),bR(i),uL(i),uR(i),vL(i),vR(i),phiL(i),phiR(i),sE1(i),sE2(i),drytol,g,sw,fw)

        ! eliminate ghost fluxes for wall
        do mw=1,3
            ! TODO: What if we use sw(mw) as fw(4,mw), allowing contiguous access
            sw(mw)   = sw(mw)*wall(mw)
            fw(1,mw) = fw(1,mw)*wall(mw) 
            fw(2,mw) = fw(2,mw)*wall(mw)
            fw(3,mw) = fw(3,mw)*wall(mw)
        enddo

        ! TODO: Perhaps split loop with if/else ixy to make it suitable for
        ! vectorization via blending?
        do mw=1,mwaves
            s(mw,i)        = sw(mw)
            fwave(1,mw,i)  = fw(1,mw)
            fwave(mu,mw,i) = fw(2,mw)
            fwave(nv,mw,i) = fw(3,mw)
        enddo
    enddo
    !$omp end simd

    !==========Capacity for mapping from latitude longitude to physical space====
    if (mcapa.gt.0) then
        do i=1-mbc,mx+mbc-1
            if (ixy.eq.1) then
                dxdc=(earth_radius*deg2rad)
            else
                dxdc=earth_radius*cos(auxr_aos(i,3))*deg2rad
            endif

            do mw=1,mwaves
                s(mw,i)=dxdc*s(mw,i)
                fwave(1,mw,i)=dxdc*fwave(1,mw,i)
                fwave(2,mw,i)=dxdc*fwave(2,mw,i)
                fwave(3,mw,i)=dxdc*fwave(3,mw,i)
            enddo
        enddo
    endif
    !===============================================================================

    !============= compute fluctuations=============================================
    amdq(:,1:3) = 0.d0
    apdq(:,1:3) = 0.d0
    do i=2-mbc,mx+mbc
        do  mw=1,mwaves
            if (s(mw,i) < 0.d0) then
                amdq(i,1:3) = amdq(i,1:3) + fwave(1:3,mw,i)
            else if (s(mw,i) > 0.d0) then
                apdq(i,1:3)  = apdq(i,1:3) + fwave(1:3,mw,i)
            else
                amdq(i,1:3) = amdq(i,1:3) + 0.5d0 * fwave(1:3,mw,i)
                apdq(i,1:3) = apdq(i,1:3) + 0.5d0 * fwave(1:3,mw,i)
            endif
        enddo
    enddo

    !!! SoA to AoS SECTION !!!
    !copy SoA arrays for ql/qr back into AoS array for consistency
    !(very unefficient!)
    do m = 1,meqn
        apdq_aos(m,:) = apdq(:,m)
        amdq_aos(m,:) = amdq(:,m)
    enddo
    !!! SoA to AoS SECTION !!!
    return
end subroutine

