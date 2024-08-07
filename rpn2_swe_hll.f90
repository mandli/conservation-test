subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

     use geoclaw_module, only: g => grav, dry_tolerance
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
     real(kind=real64) :: b_R, b_L, s_L, s_R, u_hat, c_hat, delta(3), beta(3)

     real(kind=real64) :: fwave_sum(3)

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

          ! Jump in flux
          delta(1) = hu_R - hu_L
          delta(2) = -u_R**2 + 0.5_real64 * g * h_R**2                &
                         - (-u_L**2 + 0.5_real64 * g * h_L**2) - 0.0_real64
          delta(3) = h_R * u_R * v_R - h_L * u_L * v_R

          ! HLLC
          ! s(1, i) = min(s_L, u_hat - c_hat)
          s(1, i) = s_L
          s(2, i) = 0.5_real64 * (s(1, i) + s(2, i))
          ! s(3, i) = max(s_R, u_hat + c_hat)
          s(3, i) = s_R
          beta(1) = (delta(2) - s(3, i) * delta(1)) / (s(1, i) - s(3, i))
          beta(3) = delta(1) - beta(1)
          beta(2) = delta(3) - v_L * beta(1) - v_R * beta(3)

          ! Wave 1
          fwave(1, 1, i) = beta(1)
          fwave(mu, 1, i) = beta(1) * s(1, i)
          fwave(nv, 1, i) = beta(1) * v_L

          ! Wave 2
          fwave(nv, 2, i) = beta(2)

          ! Wave 3
          fwave(1, 3, i) = beta(3)
          fwave(mu, 3, i) = beta(3) * s(3, i)
          fwave(nv, 3, i) = beta(3) * v_R

          ! ! Output waves
          fwave_sum = 0.0_real64
          do mw=1,3
               ! if (beta(mw) > 0.d0) then
               !      print "('fwave = (', i1,', ',i2,')')", mw, i
               !      print *, "beta(", mw, ") = ", beta(mw)
               !      print *, fwave(:, mw, i)
               ! end if
               fwave_sum = fwave_sum + fwave(:, mw, i)
          end do

          do m=1,3
               if (abs(delta(m) - fwave_sum(m)) > 1d-14) then !.and.         &
                   ! abs(delta(m) - fwave_sum(m)) > 0.d0) then
                    print *, "delta - (amdq + apdq) = ", delta - fwave_sum
                    print "('(', i1,', ',i2,')')", mw, i
                    ! stop
               end if
          end do

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
                    amdq(:, i) = amdq(:, i) + 0.5_real64 * fwave(:, mw, i)
                    apdq(:, i) = apdq(:, i) + 0.5_real64 * fwave(:, mw, i)
                    stop
               end if
          end do
     end do

     ! do i=2-mbc, mx+mbc
     !      if (amdq(1, i) > 0.0_real64) then
     !           print *, "Apdq, Amdq = ", amdq(1, i), apdq(1, i)
     !      end if
     ! end do

end subroutine rpn2