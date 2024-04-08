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
    use amr_module, only: t0, QR_K, tmass0, tmom0, tu20

    ! Input
    integer, intent(in) :: level, nvar, naux
    real(kind=8), intent(in) :: time
    logical, intent(in) :: rest

    ! Locals
    integer :: mptr, nx, ny, mitot, mjtot
    real(kind=8) :: hx, hy, dt
    real(kind=QR_K) :: totmass, totmom, totu2
    character(len=100) :: out_form_1 = "('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',e11.4)"
    character(len=100) :: out_form_2 = "('time t = ',e12.5,',  total mom  = ',e22.15, '  diff = ',e11.4)"
    character(len=100) :: out_form_3 = "('time t = ',e12.5,',  total u**2 = ',e22.15, '  diff = ',e11.4)"

    iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
    iaddaux(i,j) = locaux + mcapa - 1 + naux*(i-1) + naux*mitot*(j-1)

    hx      = hxposs(level)
    hy      = hyposs(level)
    dt      = possk(level)
    totmass = 0.0
    totmom = 0.0
    totu2 = 0.0

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
                    totmom = totmom + abs(alloc(iadd(2,i,j))) + abs(alloc(iadd(3,i,j)))
                    if (alloc(iadd(1,i,j)) > 0.d0) then
                        totu2 = totu2 + (alloc(iadd(2,i,j))**2 + alloc(iadd(3,i,j))**2) / alloc(iadd(1,i,j))**2
                    end if
                end do
            end do
        else
            ! with capa array:
            do j  = nghost+1, mjtot-nghost
                do i  = nghost+1, mitot-nghost
                    totmass = totmass + alloc(iadd(1,i,j))*alloc(iaddaux(i,j)) 
                    totmom = totmom + (alloc(iadd(2,i,j))**2 + alloc(iadd(3,i,j))**2) * alloc(iaddaux(i,j))
                    if (alloc(iadd(1,i,j)) > 0.d0) then
                        totu2 = totu2 + (alloc(iadd(2,i,j))**2 + alloc(iadd(3,i,j))**2) / alloc(iadd(1,i,j))**2   &
                                        * alloc(iaddaux(i,j))
                    end if
                end do
            end do
        end if

        mptr = node(levelptr,mptr)
    end if

    totmass = totmass * hx * hy
    if (abs(time - t0) < 1e-8 .and. (level == 1) .and. .not. rest) then
        tmass0 = totmass
        tmom0 = totmom
        tu20 = totu2
        print *, 'Total mass at initial time: ', tmass0
        print *, 'Total momentum at initial time: ', tmom0
        print *, 'Total u**2 at initial time: ', tu20
    endif
    write(outunit, out_form_1) time, totmass, totmass-tmass0
    write(outunit, out_form_2) time, totmom, totmom-tmom0
    write(outunit, out_form_3) time, totu2, totu2-tu20

end subroutine conck