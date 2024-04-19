subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

     use geoclaw_module, only: g => grav, dry_tolerance, rho
     use, intrinsic :: iso_fortran_env, only: real64

     implicit none

     ! Input
     integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx

     real(kind=real64), intent(inout) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: s(mwaves, 1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: ql(meqn, 1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: qr(meqn, 1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: apdq(meqn,1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: amdq(meqn,1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: auxl(maux,1-mbc:maxm+mbc)
     real(kind=real64), intent(inout) :: auxr(maux,1-mbc:maxm+mbc)

     ! Local
     integer :: m, i, mw, mu, nv
     real(kind=real64) :: h_R, h_L, hu_R, hu_L, hv_R, hv_L, u_R, u_L, v_R, v_L
     real(kind=real64) :: b_R, b_L, s_L, s_R, u_hat, c_hat

     ! Primary loop through Riemann problems
     fwave = 0.0_real64
     s = 0.0_real64
     do i = 2 - mbc, mx + mbc

          ! Direction
          if (ixy == 1) then
               mu = 2
               nv = 3
          else
               mu = 3
               nv = 2
          end if

          ! Extract Riemann problem
          h_L = qr(1, i-1)
          h_R = ql(1, i)
          hu_L = qr(mu, i-1)
          hu_R = ql(mu, i)
          b_L = auxr(1, i-1)
          b_R = auxl(1, i)
          hv_L = qr(nv, i-1)
          hv_R = ql(nv, i)

          ! Compute u, v
          if (h_R > dry_tolerance) then
               u_R = hu_R / h_R
               v_R = hv_R / h_R
          else
               u_R = 0.0_real64
               v_R = 0.0_real64
          end if
          if (h_L > dry_tolerance) then
               u_L = hu_L / h_L
               v_L = hv_L / h_L
          else
               u_L = 0.0_real64
               v_L = 0.0_real64
          end if

          ! Assume no dry states or other complexities
          s_L = u_L - sqrt(g * h_L)
          s_R = u_R + sqrt(g * h_R)
          u_hat = (sqrt(g * h_L) * u_L + sqrt(g * h_R) * u_R)         &
                         / (sqrt(g * h_R) + sqrt(g * h_L))
          c_hat = sqrt(0.5_real64 * g * (h_R + h_L))

          ! HLLC
          s(1, i) = s_L
          s(2, i) = s_R
          s(3, i) = 0.5_real64 * (s(1, i) + s(2, i))
          fwave(:, 1, i) = 0.0_real64
          fwave(:, 2, i) = 0.0_real64
          fwave(:, 3, i) = 0.0_real64

     end do

     ! Fluctuations
     amdq = 0.0_real64
     apdq = 0.0_real64
     do i = 2 - mbc, mx + mbc
          do mw = 1, mwaves
               if (s(mw, i) < 0.0_real64) then
                    amdq(:, i) = amdq(:, i) + fwave(:, mw, i)
               else if (s(mw, i) > 0.0_real64) then
                    apdq(:, i) = apdq(:, i) + fwave(:, mw, i)
               else
                    ! amdq(:, i) = amdq(:, i) + 0.5_real64 * fwave(:, mw, i)
                    ! apdq(:, i) = apdq(:, i) + 0.5_real64 * fwave(:, mw, i)
                    continue
               end if
          end do
     end do

end subroutine rpn2