! Conservation check for specified level.
! This is mostly a debugging tool and assumes grids don't overlap
!
! ******************************************************************
! conck - conservation check  for specified level
!         mostly a debugging tool
!         this assumes grids don't overlap
! 
! ******************************************************************
subroutine conck(level, nvar, naux, time, rest)

    use, intrinsic :: iso_fortran_env, only: real64, real128

    use amr_module, only: levelptr, store1, storeaux, ndilo, ndihi, ndjlo, ndjhi
    use amr_module, only: alloc, hxposs, hyposs, possk, lstart, node, nghost, outunit
    use amr_module, only: t0, init_mass, init_momentum

    ! Input
    integer, intent(in) :: level, nvar, naux
    real(kind=8), intent(in) :: time
    logical, intent(in) :: rest

    ! Locals
    integer :: mptr, nx, ny, mitot, mjtot
    real(kind=8) :: hx, hy, dt

    real(kind=real128) :: mass, momentum
    character(len=100) :: out_form_1 = "('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',e22.15)"
    character(len=100) :: out_form_2 = "('time t = ',e12.5,',  total mom  = ',e22.15, '  diff = ',e22.15)"

    ! Indexing
    iadd(ivar,i,j) = loc + ivar - 1 + nvar * ((j - 1) * mitot + i - 1)
    iaddaux(m,i,j) = locaux + m - 1 + naux * (i - 1) + naux * mitot * (j - 1)

    hx      = hxposs(level)
    hy      = hyposs(level)
    dt      = possk(level)
    mass = 0.0_real128
    momentum = 0.0_real128

    mptr = lstart(level)
    do while (mptr /= 0)
        loc    = node(store1,mptr)
        locaux = node(storeaux,mptr)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost

        if (mcapa == 0) then
            do j  = nghost+1, mjtot-nghost
                do i  = nghost+1, mitot-nghost
                    mass = mass + alloc(iadd(1,i,j)) 
                    momentum = momentum + alloc(iadd(2,i,j)) + alloc(iadd(3,i,j))
                end do
            end do
        else
            stop "Capacity array not implemented."
        end if

        mptr = node(levelptr,mptr)
    end do

    mass = mass * hx * hy
    momentum = momentum * hx * hy
    if (mass < 1d-90) mass = 0.0_real128
    if (momentum < 1d-90) momentum = 0.0_real128
    if (abs(time - t0) < 1e-8 .and. (level == 1) .and. .not. rest) then
        init_mass = mass
        init_momentum = momentum
        print *, 'Total mass at initial time: ', init_mass
        print *, 'Total momentum at initial time: ', init_momentum
    endif

    write(outunit, out_form_1) time, mass, mass - init_mass
    write(outunit, out_form_2) time, momentum, momentum - init_momentum

    ! if (abs(mass - init_mass) > 1d-20) then
    !     print *, "*** mass difference, ", mass - init_mass,"at t = ", time
    ! end if

end subroutine conck