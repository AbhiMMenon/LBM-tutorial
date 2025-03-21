! Fortran implementation of Philip Mocz's Lattice Boltzam tutorial
! Notes:
! - DQ29 LBM scheme
! - Equilibrium function for isothermal flow, i.e., constant speed of sound and specific gas constant
! - It replicates a dynamic viscosity as 
!     mu = rho *( tau - 1/2)*Delta_T
! - Delta_T (lattice time step) and Delta_x (lattice distance) are 1 (unity) for simplicity
program LBMSolver
    use, intrinsic :: iso_fortran_env, only:  dp=>real64
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    implicit none


    ! NL is the no. of DOFs in velocity phase-space, F is the number density in
    ! velocity phase-space, tau is the relaxating time. The arrays idxs, cxs, cys and
    ! refl are, respectively: numbering indices for each node in a voxel, x-velocity
    ! for each node, y-velocity for each node, and a 'reverse' numbering that implements
    ! a reflective BC.
    

    integer, parameter  ::  Nx = 400, Ny= 100, NL = 9, Nt= 2000
    real(dp), parameter ::  rho_0 = 1.100_dp, tau = 0.62_dp, PI=4.D0*DATAN(1.D0)
    
    integer             ::  idxs(9), cxs(9), cys(9), refl(9)
    integer             ::  ix,jx,kx,tx !iterators
    real(dp)            ::  weights(9), F(Nx,Ny,NL),F_eq(Nx,Ny,NL), x(Nx,Ny), y(Nx,Ny)
    real(dp)            ::  ux(Nx,Ny), uy(Nx,Ny), rho(Nx,Ny), vDotU(Nx,Ny), bound(Nx,Ny,Nl)

    logical             ::  mask(Nx,Ny)

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
    F = F + 0.01*F_eq

    x = spread([(ix, ix=0,Nx-1)], 2, Ny) 
    y = spread([(jx, jx=0,Ny-1)], 1, Nx) 

    ! Add an initial number distribution along nodal points 2, i.e., having lattice velocity
    ! (+1,0)
    F(:,:, 2) = F(:,:,2) + 2*(1. + 0.05*cos(2*PI*x/Nx*4))
    mask = ((x-Nx/4.)**2 + (y-Ny/2.)**2).lt.(Ny/4.)**2 !< logical mask for where the cylinder is 

    ! Normalize with density
    rho = sum(F(:,:,:), dim=3)
    do kx = 1,NL
        F(:,:,kx) = F(:,:,kx)*rho_0/rho
    end do
    
    call calcFlow(F,rho, ux,uy)

    where(mask)
        rho = IEEE_Value(rho, IEEE_QUIET_NAN)
        ux  = IEEE_Value(ux, IEEE_QUIET_NAN)
        uy  = IEEE_Value(uy, IEEE_QUIET_NAN)
    end where
    
    call mat_write(rho, x,y,'rho',0)
    call mat_write(ux, x,y,'ux',0)
    call mat_write(uy, x,y,'uy',0)
    

    do tx = 1,Nt
        print '("Iteration no.",i5)',tx

        ! -- step 1: drift/adection in velocity phase-space
        do kx = 1,9
            F(:,:,kx) = cshift(F(:,:,kx), shift=-cxs(kx), dim=1)
            F(:,:,kx) = cshift(F(:,:,kx), shift=-cys(kx), dim=2) 
        end do

        call calcFlow(F,rho, ux,uy)


        ! -- compute equilibrium number densities
        F_eq = 0._dp
        do kx = 1,9
            vDotU = ux*cxs(kx) + uy*cys(kx)
            F_eq(:,:,kx) = rho*weights(kx)*(1. + 3.*vDotU + 4.5*vDotU**2- 1.5*(ux**2 + uy**2))
        end do


        ! -- Store reflective BC
        do kx = 1,Nl
            where(mask)
                bound(:,:,kx) = F(:, :,refl(kx))
            end where
        end do
        
        ! -- Collision step, relax to equilibrium
        F = F+ (F_eq -F)/tau
        
        ! -- enforce reflective BC
        do kx = 1,Nl
            where(mask)
                F(:, :,kx) = bound(:,:,kx)
            end where
        end do

    end do
    
    ! -- Print end fields
    where(mask)
        rho = IEEE_Value(rho, IEEE_QUIET_NAN)
        ux  = IEEE_Value(ux, IEEE_QUIET_NAN)
        uy  = IEEE_Value(uy, IEEE_QUIET_NAN)
    end where
    
    call mat_write(rho, x,y,'rho',1)
    call mat_write(ux, x,y,'ux',1)
    call mat_write(uy, x,y,'uy',1)
    


    contains 


        !----------------------------------------------------------!
        ! Write matrix M to a file 'id.dat', ascii format for now  !
        ! M - array to write
        ! x - x meshgrid
        ! y - y meshgrid
        ! prefix - name of scalar
        ! id - time_step number
        !----------------------------------------------------------!
        subroutine mat_write(M,X,Y,prefix,id)
            implicit none
            real(dp), intent(in)          :: M(:,:), X(:,:), Y(:,:)
            integer, intent(in)           :: id
            character(len=*), intent(in)  :: prefix
            integer                       :: io, ix, jx
            character(len=25)             :: fName
            
            write (fName, '(I4)') id

            fName = 'data/'//prefix//'_'//trim(adjustl(fName))//'.dat'
            fName = trim(adjustl(fName))

            open(newunit=io, file=fName,  action='write')
            do ix = 1,size(M,dim=1)
                do jx = 1,size(M,dim=2)
                    write (io,*) X(ix,jx), Y(ix,jx),M(ix,jx)
                end do
            end do
            close(io)
        end subroutine mat_write


        subroutine calcFlow(F, rho, ux, uy)
            implicit none
            real(dp), intent(in)          :: F(:,:,:)
            real(dp), intent(out)         :: ux(:,:), uy(:,:), rho(:,:)
        
            ! -- Calculate flow variables
            rho = sum(F(:,:,:), dim=3)
            ux = 0.
            uy = 0.
            do kx = 1,Nl
                ux = ux + F(:,:,kx)*cxs(kx)
                uy = uy + F(:,:,kx)*cys(kx)
            end do
            ux = ux/rho
            uy = uy/rho
        end subroutine calcFlow

end program LBMSolver

