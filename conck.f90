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

    use amr_module, only: levelptr, store1, storeaux, ndilo, ndihi, ndjlo, ndjhi
    use amr_module, only: alloc, hxposs, hyposs, possk, lstart, node, nghost, outunit

    ! Input
    integer, intent(in) :: level, nvar, naux
    real(kind=8), intent(in) :: time
    logical, intent(in) :: rest

    ! Locals
    integer :: mptr, nx, ny, mitot, mjtot
    real(kind=8) :: hx, hy, dt, totmass
    character(len=100) :: out_form = "('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',e11.4)"

    iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
    iaddaux(i,j) = locaux + mcapa - 1 + naux*(i-1) + naux*mitot*(j-1)

    hx      = hxposs(level)
    hy      = hyposs(level)
    dt      = possk(level)
    totmass = 0.d0

    mptr = lstart(level)
    if (mptr /= 0) then
        loc    = node(store1,mptr)
        locaux = node(storeaux,mptr)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost

        if (mcapa == 0) then
            do j  = nghost+1, mjtot-nghost
                do i  = nghost+1, mitot-nghost
                    totmass = totmass + alloc(iadd(1,i,j)) 
                end do
            end do
        else
            ! with capa array:
            do j  = nghost+1, mjtot-nghost
                do i  = nghost+1, mitot-nghost
                    totmass = totmass + alloc(iadd(1,i,j))*alloc(iaddaux(i,j)) 
                end do
            end do
        end if

        mptr = node(levelptr,mptr)
    end if

    totmass = totmass * hx * hy
    if (time == t0 .and. (level == 1) .and. .not. rest) then
        tmass0 = totmass
        print *, 'Total mass at initial time: ', tmass0
    endif
    write(outunit, out_form) time, totmass, totmass-tmass0

end subroutine conck