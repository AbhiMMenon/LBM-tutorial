!----------------------------------------------------------------------------------------
! Fortran implementation of Philip Mocz's Lattice Boltzam tutorial: flow around a
! cylinder (in 2D) in a periodic domain. Reflective BCs are used for the cylinder,
! non-reflective BCs to inlet and outlet in this implementation.
!
! Notes:
! - DQ29 LBM
! - Equilibrium function for isothermal flow, i.e., constant speed of sound and specific
! gas constant
! - It replicates a dynamic viscosity as 
!     mu = rho*( tau - 1/2)*Delta_T
!----------------------------------------------------------------------------------------
program LBMSolver
    use, intrinsic :: iso_fortran_env, only:  dp=>real64
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    implicit none


    ! NL is the no. of DOFs in velocity phase-space (9), F is the density function, tau
    ! is the relaxating time. The arrays idxs, cxs, cys and refl are, respectively:
    ! numbering indices for each node in a voxel, x-velocity for each node, y-velocity
    ! for each node, and a 'reverse' numbering that implements a reflective BC.
    

    integer, parameter  ::  Nx = 800, Ny= 200, NL = 9, Nt= 20000, Nsave=20
    real(dp), parameter ::  rho_0 = 1.00_dp, tau = 0.7_dp, PI=4.D0*DATAN(1.D0)
    
    integer             ::  idxs(9), cxs(9), cys(9), refl(9)
    integer             ::  kx,tx,filex ! ix , shift,  iterators
    real(dp)            ::  weights(9), F(Nx,Ny,NL),F_eq(Nx,Ny,NL), &
                            x(Nx,Ny), y(Nx,Ny), &
                            ux(Nx,Ny), uy(Nx,Ny), rho(Nx,Ny), vDotU(Nx,Ny), &
                            bound(Nx,Ny,Nl), omega(Nx,Ny), &
                            nanVal
    logical             ::  mask(Nx,Ny)

    nanVal = IEEE_Value(nanVal, IEEE_QUIET_NAN)
    ! The numbering system here for D2Q9 here is as follows, it differs from Peter's
    ! python version
    !                            7   3   6
    !                              \ | /       
    !                            4 - 1 - 2
    !                              / | \
    ! y                          8   5   9
    ! ↑ 
    ! . → x

    idxs = [ 1, 2, 3, 4, 5, 6, 7, 8, 9]
    refl = [ 1, 4, 5, 2, 3, 8, 9, 6, 7]
    cxs  = [ 0, 1, 0,-1, 0, 1,-1,-1, 1]
    cys  = [ 0, 0, 1, 0,-1, 1, 1,-1,-1]
    weights = [16., 4., 4., 4., 4., 1., 1.,1.,1.]/36._dp

    ! Initial conditions
    F = 1.0_dp
    call random_number(F_eq)
    F = F + 0.1*F_eq

    x = spread([(kx, kx=0,Nx-1)], 2, Ny) 
    y = spread([(kx, kx=0,Ny-1)], 1, Nx) 

    ! Initial distribution fuction for along nodal points 2, i.e., having lattice velocity
    ! (+1,0)
    F(:,:, 2) = F(:,:,2)+ 2*(1. + 0.5*cos(2*PI*x/Nx*4))
    mask = ((x-Nx/4.)**2 + (y-Ny/2.)**2).lt.(Ny/8.)**2 !< logical mask for where the cylinder is 

    ! Normalize with density
    rho = sum(F, dim=3)
    do kx = 1,NL
        F(:,:,kx) = F(:,:,kx)*rho_0/rho
    end do
    
    call calcFlow(F,rho, ux,uy)


    filex=0 
    do tx = 1,Nt
        print '("Iteration no.",i5)',tx

        ! -- step 1: drift/adection in velocity phase-space. 

        ! The use of cshift for small arrays here as lattice velocity and lattice side
        ! lenghts are unity, i.e, there is a unit shift of the distribution function
        ! determined by cxs and cys per time advancement

        do kx = 1,9
            F(:,:,kx) = cshift(F(:,:,kx), shift=-cxs(kx), dim=1)
            F(:,:,kx) = cshift(F(:,:,kx), shift=-cys(kx), dim=2) 

        ! ** Manual version for larger arrays prevents creation of temp arrays, check
        ! ** with -Warray-temporaries, re-uses the F_eq buffer
            
        !   shift = cxs(kx)
        !   do ix=2,Nx-1
        !       F_eq(ix,:,kx) = F(ix-shift,:,kx)
        !   end do 
        !   F_eq(Nx,:,kx) = F(1,:,kx)  ! < periodic BC
        !   F_eq(1,:,kx)  = F(Nx,:,kx) ! < "
        !   F(:,:,kx)     = F_eq(:,:,kx)
        !   
        !   shift = cys(kx)
        !   do ix=2,Ny-1
        !       F_eq(:,ix,kx) = F(:,ix-shift,kx)
        !   end do 
        !   F_eq(:,Ny,kx) = F(:,1,kx) ! < periodic BC
        !   F_eq(:,1,kx)  = F(:,Ny,kx) ! < "
        !  
        !   F(:,:,kx) = F_eq(:,:,kx)
            

        end do

        call calcFlow(F,rho, ux,uy)


        ! -- compute equilibrium fuction, includes cs=1/sqrt(3)
        F_eq = 0._dp
        do kx = 1,9
            vDotU = ux*cxs(kx) + uy*cys(kx)
            F_eq(:,:,kx) = rho*weights(kx)*(1. + 3.*vDotU + 4.5*vDotU**2- 1.5*(ux**2 + uy**2))
        end do


        ! -- Store reflections
        do kx = 1,Nl
            where(mask)
                bound(:,:,kx) = F(:, :,refl(kx))
            end where
        end do

        ! -- Collision step, relax to equilibrium. The tau here represents a viscosity that is 
        ! -- related to c_s (speed of sound) as nu = c_s^2*(tau-0.5)

        F = F+ (F_eq -F)/tau 
        
        ! -- apply reflections to cylinder
        do kx = 1,Nl
            where(mask)
                F(:, :,kx) = bound(:,:,kx)
            end where
        end do

        ! -- testing: apply non-reflective BC on inlet and outlet
        F(Nx,:,[7,4,8]) = F(Nx-1,:,[7,4,8])
        F(1,:,[2,6,9]) = F(2,:,[2,6,9])


        if (mod(tx,Nsave)==0) then 
            ! -- mask out the circle for plotting
            where(mask)
                rho = nanVal
                ux  = nanVal
                uy  = nanVal
                omega = nanVal
            end where
            
        !   call mat_write(rho, 'rho', filex)
        !   call mat_write( ux,  'ux', filex)
        !   call mat_write( uy,  'uy', filex)
            call mat_write( omega,  'w', filex)
            filex = filex + 1
        end if

    end do
    


    contains 


        !------------------------------------------------------------------------------!
        ! Write matrix M to a file 'prefix_id.dat', 
        ! M - array to write
        ! prefix - name of scalar
        ! id - time_step number
        !------------------------------------------------------------------------------!
        subroutine mat_write(M,prefix,id)
            implicit none
            real(dp), intent(in)          :: M(:,:)
            integer, intent(in)           :: id
            character(len=*), intent(in)  :: prefix
            integer                       :: io
            character(len=25)             :: fName
            
            write (fName, '(I4)') id

            fName = 'data/'//prefix//'_'//trim(adjustl(fName))//'.dat'
            fName = trim(adjustl(fName))

            open(newunit=io, file=fName,  action='write',form='unformatted',access='stream',status='replace')
                write(io) real(M)
            close(io)
        end subroutine mat_write

        !------------------------------------------------------------------------------!
        ! Compute density, ux and uy as weighted sums of particle numbers in  velocity !
        ! phase space                                                                  !
        !------------------------------------------------------------------------------!
        subroutine calcFlow(F, rho, ux, uy)
            implicit none
            real(dp), intent(in)    ::  F(:,:,:)
            real(dp), intent(out)   ::  ux(:,:), uy(:,:), rho(:,:)
            integer                 ::  ix, jx
        
            ! -- Calculate flow variables
            ux = 0.
            uy = 0.
            omega = 0.
            rho = sum(F,dim=3)
            do kx = 1,Nl
               !rho = rho+ F(:,:, kx)
                ux  = ux + F(:,:,kx)*cxs(kx)
                uy  = uy + F(:,:,kx)*cys(kx)
            end do
            
            where(mask)
                ux  = 0.
                uy  = 0.
            end where

            ux = ux/rho
            uy = uy/rho

            do jx = 2,Ny-2
                do ix = 2,Nx-2
                    omega(ix,jx)  = (ux(ix  ,jx+1) - ux(ix  ,jx-1)) - &
                                    (uy(ix+1,jx  ) - uy(ix-1,jx  ))
                end do
            end do
            omega = omega * 0.5
           
           !omega = uy(3:Nx,:)- uy(1:Nx-2,:) + uy(:,3:Ny)- uy(:,1:Ny-2)
           !omega = (cshift(ux,1,dim=2)- cshift(ux,-1,dim=2)) - &
           !        (cshift(uy,1,dim=1)- cshift(uy,-1,dim=1))

        end subroutine calcFlow

end program LBMSolver

