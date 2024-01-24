module common_hh
  implicit none
  integer :: nx, ny, i, j, k, im, jm, ni, nj
  integer ::j_surf, j_bound, j_periodic, j_snu, &
             j_uadvec, j_qgive, j_cadvec
  integer ::j_dens, j_dens0, j_temp, j_temp0
  real(8)::g, snu, time, rho, skt, beta_t, skc, diam, stime_surf, alpha_surf
  real(8)::dt, dx, dx2
  real(8)::tmp0, con0
  real(8)::t_vol, q_ave, hs_ave_c, u_ave_c
  real(8)::t_length, t_vol_ini, st_dens, st_temp

  real(8)::qp, slope, hs_ave, u_ave, us_ave, kappa, disam, ks, h_down, h_up_boundary
  real(8)::qp0, q_stt, q_trn, qp1, q_relax, h_plus, q_ref
  real(8)::dhdt_down, h_down_old
  integer::j_hdw, j_hup, j_qadj
  real(8)::xi1, z0, cd0, cd1
  real(8)::al_ep

end module common_hh

module common_geom
  implicit none
  real(8), dimension(:), allocatable::xi, xxi, dy, ddy
  real(8), dimension(:), allocatable::eta, eta_up, x_up
  integer, dimension(:, :), allocatable::obst
end module common_geom

module common_hyd
  implicit none
  real(8), dimension(:), allocatable::h, hn, hs, uba, usta, qc, ypsurf, cd
  real(8), dimension(:, :), allocatable::yu, yun, yv, yvn, v1, v2, omega
  real(8), dimension(:, :), allocatable::yt, ytn, yp, ypn, yc, ycn
  real(8), dimension(:, :), allocatable::yuo, yvo
  real(8), dimension(:, :), allocatable::snu_t, snu_tx
  real(8), dimension(:), allocatable::u00
end module common_hyd

module common_grad
  implicit none
  real(8), dimension(:, :), allocatable::gux, guy, gvx, gvy
  real(8), dimension(:, :), allocatable::gcx, gcy, gtx, gty
  real(8), dimension(:, :), allocatable::gux_n, guy_n, gvx_n, gvy_n
  real(8), dimension(:, :), allocatable::gcx_n, gcy_n, gtx_n, gty_n
end module common_grad

!================================================

module alloc_var_m
  use common_geom
  use common_hyd
  use common_grad
  implicit none
contains
!------------------------------------------------
  subroutine alloc_var(im, jm)
    integer::im, jm
    integer::i, j
    i = im; j = jm
    allocate (xi(0:j), xxi(0:j), dy(0:j), ddy(0:j), eta(0:i), eta_up(0:i), x_up(0:i))
    allocate (h(0:i), hn(0:i), hs(0:i), usta(0:i), qc(0:i), ypsurf(0:i) &
              , cd(0:i), uba(0:i))
    allocate (yu(0:i, 0:j), yun(0:i, 0:j), yv(0:i, 0:j), yvn(0:i, 0:j))
    allocate (yt(0:i, 0:j), ytn(0:i, 0:j), yp(0:i, 0:j), ypn(0:i, 0:j) &
              , yc(0:i, 0:j), ycn(0:i, 0:j), obst(0:i, 0:j))
    allocate (v1(0:i, 0:j), v2(0:i, 0:j), omega(0:i, 0:j))
    allocate (yuo(0:i, 0:j), yvo(0:i, 0:j))
    allocate (snu_t(0:i, 0:j), snu_tx(0:i, 0:j), u00(0:j))
    allocate (gux(0:i, 0:j), guy(0:i, 0:j), gvx(0:i, 0:j), gvy(0:i, 0:j))
    allocate (gcx(0:i, 0:j), gcy(0:i, 0:j), gtx(0:i, 0:j), gty(0:i, 0:j))
    allocate (gux_n(0:i, 0:j), guy_n(0:i, 0:j), gvx_n(0:i, 0:j) &
              , gvy_n(0:i, 0:j), gcx_n(0:i, 0:j), gcy_n(0:i, 0:j) &
              , gtx_n(0:i, 0:j), gty_n(0:i, 0:j))
  end subroutine alloc_var
end module alloc_var_m

!-------------------------------------------------
module ss_nu_m
  use common_hh
  use common_hyd
  implicit none
contains
! ------------------------------------------------
  function ss_nu(usta, hs, xi)
!-------------------------------------------------
    real(8)::ss_nu, usta, hs, xi
!
    if (j_snu == 1) then
      ss_nu = kappa*abs(usta)*hs*xi*(1.-xi)*al_ep
!      write(44,'(5f12.7)')ss_nu,usta,hs,xi,al_ep
    else
      ss_nu = snu*al_ep
!      write(44,'(f12.7)')ss_nu
    end if
  end function ss_nu
end module ss_nu_m

!================================================
module initial_val_m
  use common_hh
  use common_geom
  use common_hyd
  use common_grad
  use ss_nu_m
  implicit none
contains
!------------------------------------------------
  subroutine initial_0
    xi = 0.; xxi = 0.; dy = 0.; ddy = 0.; eta = 0.; eta_up = 0.; x_up = 0.
    h = 0.; hn = 0.; hs = 0.; uba = 0.; usta = 0.; qc = 0.; ypsurf(0:i) = 0.; cd = 0.
    yu = 0.; yun = 0.; yv = 0.; yvn = 0.; yt = 0.; ytn = 0.; yp = 0.; ypn = 0.
    v1 = 0.; v2 = 0.; omega = 0.
    yuo = 0.; yvo = 0.
    gux = 0.; guy = 0.; gvx = 0.; gvy = 0.
    gux_n = 0.; guy_n = 0.; gvx_n = 0.; gvy_n = 0.
    gcx = 0.; gcy = 0.; gtx = 0.; gty = 0.
    gcx_n = 0.; gcy_n = 0.; gtx_n = 0.; gty_n = 0.
    obst = 0
    snu_t = 0.; snu_tx = 0.; u00 = 0.
    t_vol = 0.; q_ave = 0.; hs_ave_c = 0.; u_ave_c = 0.; t_vol_ini = 0.; t_length = 0.

  end subroutine initial_0

  subroutine initial
    real(8)::zz, hs_up, us_hp, sb, q_ave_c

    do j = 0, ny + 1
      dy(j) = 1./float(ny)
    end do
    xi(0) = 0.
    do j = 1, ny + 1
      xi(j) = xi(j - 1) + dy(j)
    end do

    do j = 0, ny
      ddy(j) = (dy(j) + dy(j + 1))*.5
    end do
    ddy(ny + 1) = dy(ny)
    do j = 1, ny + 1
      xxi(j) = (xi(j) + xi(j - 1))*.5
    end do

    xi1 = dy(1)*.5

!      do j=0,ny+1
!       write(44,'(i4,4f10.5)') j,dy(j),ddy(j),xi(j),xxi(j)
!      end do

!
! Initial u,usta,cd,snu_tx
!
    sb = log(hs_ave/z0) - 1.
    do j = 1, ny
      zz = xxi(j)*hs_ave
!       u00(j)=u_ave*log(zz/z0)/sb
      u00(j) = 0.
    end do
    u00(0) = 0.
    u00(ny + 1) = u00(ny)

!      write(44,*) z0,sb,hs_ave

!      do j=1,ny
!       write(44,'(i4,3f10.4)') j,xxi(j),xxi(j)*hs_ave,u00(j)
!      end do
!
    cd1 = 1./(1./kappa*log(xi1*hs_ave/z0))
    us_ave = u00(1)*cd1

    q_ave_c = 0.
    do i = 0, nx + 1
      cd(i) = cd1**2
      usta(i) = 0.
      do j = 0, ny + 1
        yu(i, j) = 0.
        yun(i, j) = yu(i, j)
        snu_tx(i, j) = 0.; snu_t(i, j) = 0.
        gux(i, j) = 0.; guy(i, j) = 0.; gvx(i, j) = 0.; gvy(i, j) = 0.
        gux_n(i, j) = 0.; guy_n(i, j) = 0.; gvx_n(i, j) = 0.; gvy_n(i, j) = 0.
      end do
    end do

  end subroutine initial

end module initial_val_m

!================================================
module bound_p_m
  use common_hh
  use common_geom
  use common_hyd
  implicit none
contains
!------------------------------------------------
  subroutine bound_p(smg_g)
!
! 表面張力の計算 (j_surf==1の時のみ)
!
    real(8)::smg_g
!
! ----- Surface Pressure ------
!
    do i = 2, nx - 1
      ypsurf(i) = smg_g*(hn(i - 1) - 2.*hn(i) + hn(i + 1))/dx2
      ypn(i, ny + 1) = 2.*ypsurf(i) - ypn(i, ny)
    end do

    if (j_bound == 2 .and. j_periodic == 1) then
! 開水路周期境界条件
      ypn(1, ny + 1) = ypn(nx - 1, ny + 1)
      ypn(nx, ny + 1) = ypn(2, ny + 1)
    end if

! 壁条件(上流，左壁)
    if (j_bound == 1 .or. j_bound == 4) then
      ypn(1, ny + 1) = ypn(2, ny + 1)
    end if

! 壁条件(下流，右壁)
    if (j_bound == 1 .or. j_bound == 3) then
      ypn(nx, ny + 1) = ypn(nx - 1, ny + 1)
    end if
  end subroutine bound_p

  subroutine pbound
! j_bound=2 .and. j_periodic==1 周期境界の時のみ
    do j = 0, ny + 1
      ypn(1, j) = ypn(nx - 1, j)
      ypn(nx, j) = ypn(2, j)
    end do
  end subroutine pbound
end module bound_p_m

!================================================
module sor_m
  use common_hh
  use common_geom
  use common_hyd
  implicit none
contains
!------------------------------------------------
  subroutine sor(sorerr, lsor, soralpha, lmax &
                 , jc_in_l, jc_in_r, jc_in_b, jc_in_u &
                 , jc_l_s, jc_l_e, jc_r_s, jc_r_e, ic_b_s, ic_b_e, ic_u_s, ic_u_e &
                 , c_bound_l, c_bound_r, c_bound_b, c_bound_u)

    real(8)::sorerr, soralpha, div, ww, err, errx
    integer:: lsor, l, lmax

    real(8)::work(0:im, 0:jm), a_p(0:im, 0:jm), a_f(0:im, 0:jm)
    real(8)::a_w(0:im, 0:jm), a_e(0:im, 0:jm), a_n(0:im, 0:jm), a_s(0:im, 0:jm)
    real(8)::a_wx, a_ex
    real(8)::af_old, af_new
    integer::i1, i2
    integer::jc_in_l, jc_in_r, jc_in_b, jc_in_u
    integer::jc_l_s(50), jc_l_e(50), jc_r_s(50), jc_r_e(50), &
              ic_b_s(50), ic_b_e(50), ic_u_s(50), ic_u_e(50)
    real(8)::c_bound_l(50), c_bound_r(50), c_bound_b(50), c_bound_u(50)
    integer:: m

    div = 0.
    work = 0.; a_e = 0.; a_w = 0.; a_s = 0.; a_n = 0.; a_p = 0.; a_f = 0.

    i1 = 2; i2 = nx - 1

    do i = i1, i2
      do j = 1, ny
        if (obst(i, j) == 0) then
          work(i, j) = ypn(i, j)

! 壁境界(閉鎖水域または上流壁)
          if (i == i1 .and. &
              (((j_bound == 1 .or. j_bound == 4) .and. i == 2) .or. obst(i - 1, j) == 1)) then
            a_w(i, j) = 0.
          else
            a_w(i, j) = (hs(i) + hs(i - 1))/(dx2*2.)
          end if
! 壁境界(閉鎖水域または下流壁)
          if (i == i2 .and. &
              (((j_bound == 1 .or. j_bound == 3) .and. i == nx - 1) .or. obst(i + 1, j) == 1)) then
            a_e(i, j) = 0.
          else
            a_e(i, j) = (hs(i + 1) + hs(i))/(dx2*2.)
          end if

!
          if (obst(i, j - 1) == 1 .or. j == 1) then
            a_s(i, j) = 0.
          else
            a_s(i, j) = 1./(hs(i)*dy(j)**2)
          end if
          if (obst(i, j + 1) == 1 .or. j == ny) then
            a_n(i, j) = 0.
          else
            a_n(i, j) = 1./(hs(i)*dy(j)**2)
          end if

          a_p(i, j) = a_e(i, j) + a_w(i, j) + a_n(i, j) + a_s(i, j)
          if (j == ny) a_p(i, ny) = a_p(i, ny) + 2./(hs(i)*dy(ny)**2)
          div = (yu(i, j)*(hs(i) + hs(i + 1)) &
                 - yu(i - 1, j)*(hs(i) + hs(i - 1)))/(dx*2.) + (v1(i, j) - v1(i, j - 1))/dy(j)
          a_f(i, j) = -div/dt*rho &
                      - (omega(i, j)*(ypn(i + 1, j + 1) + ypn(i, j + 1) - ypn(i + 1, j - 1) - ypn(i, j - 1)) &
                         + omega(i - 1, j)*(ypn(i, j + 1) + ypn(i - 1, j + 1) - ypn(i, j - 1) - ypn(i - 1, j - 1))) &
                      /(4.*dx*dy(j))
          if (j_dens == 1) then
            if (j == 1) then
              a_f(i, j) = a_f(i, j) + rho*g/ddy(1)*(yc(i, j + 1) - yc(i, j))
            else if (j == ny) then
              a_f(i, j) = a_f(i, j) + rho*g/ddy(ny - 1)*(yc(i, j) - yc(i, j - 1))
            else
              a_f(i, j) = a_f(i, j) + rho*g/(2.*dy(j))*(yc(i, j + 1) - yc(i, j - 1))
            end if
          end if
          if (j_temp == 1) then
            if (j == 1) then
              a_f(i, j) = a_f(i, j) + rho*g/ddy(1)*(yt(i, j) - yt(i, j + 1))*beta_t
            else if (j == ny) then
              a_f(i, j) = a_f(i, j) + rho*g/ddy(ny - 1)*(yt(i, j - 1) - yt(i, j))*beta_t
            else
              a_f(i, j) = a_f(i, j) + rho*g/(2.*dy(j))*(yt(i, j - 1) - yt(i, j + 1))*beta_t
            end if
          end if
        end if
      end do
    end do
!
    do l = 1, lsor
      err = 0.0
      do i = i1, i2
        do j = 1, ny
          if (obst(i, j) == 1) then
            work(i, j) = 0.
          else
            ww = (a_e(i, j)*ypn(i + 1, j) + a_w(i, j)*ypn(i - 1, j) &
                  + a_n(i, j)*ypn(i, j + 1) + a_s(i, j)*ypn(i, j - 1) + a_f(i, j)) &
                 /a_p(i, j)
            work(i, j) = (1.0 - soralpha)*ypn(i, j) + soralpha*ww
            errx = abs(work(i, j) - ypn(i, j))
            err = err + errx
            ypn(i, j) = work(i, j)
          end if
        end do
      end do
      if (err < sorerr) goto 40
!
! 底面と水面
      do i = i1, i2
        ypn(i, 0) = ypn(i, 1)
! 水面の圧力(固定の時)
        if (j_surf == 0 .or. (j_surf == 1 .and. time < stime_surf)) ypn(i, ny + 1) = -ypn(i, ny)
      end do

!固定物の上面・下面
      do i = i1, i2
        do j = 1, ny - 1
          if (j_dens == 1 .and. obst(i, j) == 1 .and. obst(i, j + 1) == 0) &
            ypn(i, j) = ypn(i, j + 1) - rho*g*hs(i) &
                        *(beta_t*(yt(i, j + 1) - tmp0) - (yc(i, j + 1) - con0))*dy(j)
          if (j_temp == 1 .and. obst(i, j) == 1 .and. obst(i, j - 1) == 0) &
            ypn(i, j) = ypn(i, j - 1) + rho*g*hs(i) &
                        *(beta_t*(yt(i, j - 1) - tmp0) - (yc(i, j - 1) - con0))*dy(j)
        end do
      end do

! 左(上流)壁の場合
      if (j_bound == 1 .or. j_bound == 4) then
        do j = 0, ny + 1
          ypn(1, j) = ypn(2, j)
        end do
      end if

! 右(下流)壁の場合
      if (j_bound == 1 .or. j_bound == 3) then
        do j = 0, ny + 1
          ypn(nx, j) = ypn(nx - 1, j)
        end do
      end if

! 開放境界で周期境界
      if (j_bound == 2 .and. j_periodic == 1) then
        do j = 0, ny + 1
          ypn(1, j) = ypn(nx - 1, j)
          ypn(nx, j) = ypn(2, j)
        end do
      end if

!下流側開放
      if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 4) then
        do j = 0, ny + 1
          ypn(nx, j) = 0.
        end do
      end if

!上流開放
      if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 3) then
        do j = 0, ny + 1
          ypn(1, j) = 0.
        end do
      end if

! 固定物の左右境界
      do i = i1, i2
        do j = 1, ny
          if (obst(i, j) == 0) then
            if (i < nx - 1 .and. obst(i + 1, j) == 1) ypn(i + 1, j) = ypn(i, j)
            if (i > 2 .and. obst(i - 1, j) == 1) ypn(i - 1, j) = ypn(i, j)
          end if
        end do
      end do

    end do
40  continue
    lmax = l

!     write(44,*) 'sor',time
!     do i=i1,i2
!      write(44,'(i3,20e12.3)') i,(ypn(i,j),j=1,ny)
!     end do

  end subroutine sor
end module sor_m

!===============================================
module rhs_m
  use common_hh
  use common_geom
  use common_hyd
  implicit none
contains
!------------------------------------------------
  subroutine rhs

    real(8)::deds, dhsds, dhds, hs_up, dpds, dpdg, dd, tvp, cvp, gbt
    integer::i1, i2

!     write(44,*) 'rhs(1)'
!     write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)

!閉鎖水域
    if (j_bound == 1) then
      i1 = 2; i2 = nx - 2
!開放水域
    else if (j_bound == 2) then
      i1 = 1; i2 = nx - 1
!下流閉鎖
    else if (j_bound == 3) then
      i1 = 1; i2 = nx - 2
!上流閉鎖
    else
      i1 = 2; i2 = nx - 1
    end if

!     write(44,*) 'rhs'
    do i = i1, i2
!      deds=(eta(i+1)-eta(i))/dx
!      dhsds=(hs(i+1)-hs(i))/dx
      dhds = (h(i + 1) - h(i))/dx
      hs_up = (hs(i + 1) + hs(i))*.5
!      if(i<=2) then
!        write(44,*) ' i  ,deds,      dhsds,      dhds,       hs_up' &
!         ,'          h(i+1)   h(i)'
!        write(44,'(i3,6e12.3)') &
!           i,deds,dhsds,dhds,hs_up,h(i+1),h(i)
!        write(44,*) ' j   omega(i,j), dpds,       dpdg,       dd'   &
!           ,'        yu(i,j)     yun(i,j)'
!       end if
      do j = 1, ny
!       omega(i,j)=deds+xxi(j)*dhsds
        dpds = (ypn(i + 1, j) - ypn(i, j))/dx
        dpdg = (ypn(i + 1, j + 1) + ypn(i, j + 1) - ypn(i + 1, j - 1) - ypn(i, j - 1))/(4.*dy(j))
        dd = dpds - dpdg*omega(i, j)/hs_up
        yun(i, j) = yu(i, j) - (dd/rho + g*dhds)*dt
!       if(i<=2) write(44,'(i4,6e12.3)') &
!         j,omega(i,j),dpds,dpdg,dd,yu(i,j),yun(i,j)
        if (obst(i, j) == 1 .or. obst(i + 1, j) == 1) yun(i, j) = 0.
      end do
    end do

!     do i=1,2
!      write(44,'(i3,12e12.3)') i,(yu(i,j),j=1,ny,2)
!      write(44,'(i3,12e12.3)') i,(yun(i,j),j=1,ny,2)

    do i = i1, i2 + 1
      do j = 1, ny
        if (j_temp == 1 .or. j_dens == 1) then
          if (j_temp == 1) then
            tvp = (yt(i, j) + yt(i, j + 1))*.5
          else
            tvp = 0.
          end if
          if (j_dens == 1) then
            cvp = (yc(i, j) + yc(i, j + 1))*.5
          else
            cvp = 0.
          end if
          gbt = g*(beta_t*(tvp - tmp0) - (cvp - con0))
        else
          gbt = 0.
        end if
        if (j == ny .and. (j_surf == 0 .or. (j_surf == 1 .and. time < stime_surf))) then
          yvn(i, j) = 0.
        else
          yvn(i, j) = yv(i, j) + &
                      (gbt - (ypn(i, j + 1) - ypn(i, j))/(rho*ddy(j)*hs(i)))*dt
        end if
        if (obst(i, j) == 1 .or. obst(i, j + 1) == 1) yvn(i, j) = 0.
!       write(44,'(2i5,5f12.5)') i,j,omega(i,j),yun(i,j),yvn(i,j)
      end do
    end do

!     write(44,*) 'rhs(2)'
!     write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)

  end subroutine rhs

!**************************************
  subroutine omg
!**************************************
!
    real(8)::deds, dhsds

    omega = 0.; deds = 0.; dhsds = 0.

    do i = 1, nx - 1
      deds = (eta(i + 1) - eta(i))/dx
      dhsds = (hs(i + 1) - hs(i))/dx
      do j = 1, ny
        omega(i, j) = deds + xxi(j)*dhsds
        if (obst(i, j) == 1 .or. obst(i + 1, j) == 1) omega(i, j) = 0.
      end do
    end do

  end subroutine omg

end module rhs_m

!================================================
module diffusion_m
  use common_hh
  use common_geom
  use common_hyd
  implicit none
contains
!------------------------------------------------
  subroutine diffusion

    real(8)::hs_up, d1, d2, d3, d4, d2x
    real(8)::uvis(im, jm), vvis(im, jm), tvis(im, jm), cvis(im, jm)
    real(8)::snu_up, snu_vp
    real(8)::omega_vp, dodx, dhdx
    integer::i1, i2

    uvis = 0.; vvis = 0.; tvis = 0.; cvis = 0.; omega_vp = 0.; dodx = 0.; dhdx = 0.
    snu_up = 0.; snu_vp = 0.

!     write(44,*) 'diffusion(1)'
!     write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)
!
!閉鎖水域
    if (j_bound == 1) then
      i1 = 2; i2 = nx - 2
!開放水域
    else if (j_bound == 2) then
      i1 = 1; i2 = nx - 1
!下流閉鎖
    else if (j_bound == 3) then
      i1 = 1; i2 = nx - 2
!上流閉鎖
    else
      i1 = 2; i2 = nx - 1
    end if

    do i = i1, i2
!!     write(44,'(i5,3f10.4)') i,hs(i),yu(i,0),yu(i,1)
      hs_up = (hs(i) + hs(i + 1))*.5
      do j = 2, ny
        d1 = (yun(i + 1, j) - 2.*yun(i, j) + yun(i - 1, j))/dx2
        d2 = (yun(i, j + 1) - 2.*yun(i, j) + yun(i, j - 1)) &
             /(dy(j)*hs_up)**2*(1.+omega(i, j)**2)
!q      d2=(yun(i,j+1)-2.*yun(i,j)+yun(i,j-1))/(dy(j)**2*hs_up**2)
        d3 = -omega(i, j)/(hs_up*2.*dx*dy(j))*( &
             yun(i + 1, j + 1) - yun(i - 1, j + 1) - yun(i + 1, j - 1) + yun(i - 1, j - 1))
        d4 = (-(omega(i + 1, j) - omega(i - 1, j))/2.-(hs(i + 1) - hs(i)) + &
              omega(i, j)/hs_up*(hs(i + 1) - hs(i))) &
             *(yun(i, j + 1) - yun(i, j - 1))/(2.*dy(j)*dx*hs_up)
        snu_up = (snu_t(i, j) + snu_t(i + 1, j))*.5
        uvis(i, j) = snu_up*(d1 + d2 + d3 + d4)
      end do
      j = 1
      d1 = (yun(i + 1, j) - 2.*yun(i, j) + yun(i - 1, j))/dx2
      snu_up = (snu_t(i, j) + snu_t(i + 1, j))*.5
      d2x = snu_up/(dy(j)**2*hs_up**2)*(yun(i, j + 1) - yun(i, j)) &
            - cd(i)*yun(i, j)*abs(yun(i, j))/(dy(j)*hs_up)
      d3 = -omega(i, j)/(hs_up*2.*dx*dy(j))*( &
           yun(i + 1, j + 1) - yun(i - 1, j + 1) - yun(i + 1, j - 1) + yun(i - 1, j - 1))
      d4 = -1./(hs_up*4.*dx*dy(j))*(omega(i + 1, j) - omega(i - 1, j))* &
           (yun(i, j + 1) - yun(i, j - 1))
      uvis(i, j) = snu_up*(d1 + d3 + d4) + d2x
    end do

!     write(44,*) 'diffusion(2)'
!     write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)
!
    do i = i1, i2
      do j = 1, ny
        yun(i, j) = yun(i, j) + uvis(i, j)*dt
      end do
    end do
!
!     write(44,*) 'diffusion(3)'
!     write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)
!
    do j = 1, ny
      do i = i1, i2 + 1
        snu_vp = (snu_t(i, j) + snu_t(i, j + 1))*.5
        omega_vp = (omega(i, j + 1) + omega(i, j) + omega(i - 1, j + 1) + omega(i - 1, j))*.25
        dodx = (omega(i, j + 1) + omega(i, j) - omega(i - 1, j + 1) + omega(i - 1, j))/(2.*dx)
        dhdx = (hs(i + 1) - hs(i - 1))/(2.*dx)
        d1 = (yvn(i + 1, j) - 2.*yvn(i, j) + yvn(i - 1, j))/dx2
        d2 = (yvn(i + 1, j + 1) - yvn(i - 1, j + 1) - yvn(i + 1, j - 1) + yvn(i - 1, j - 1))* &
             omega_vp/(hs(i)*(dy(j + 1) + dy(j))*dx)
        d3 = (1.+omega_vp**2)/(hs(i)**2*ddy(j)**2)* &
             (yvn(i, j + 1) - 2.*yvn(i, j) + yvn(i, j - 1))
        d4 = (-(dodx - dhdx/hs(i))/hs(i) + omega_vp/hs(i)**2*dhdx) &
             *(yvn(i, j + 1) - yvn(i, j - 1))/(dy(j + 1) + dy(j))
        vvis(i, j) = snu_vp*(d1 + d2 + d3 + d4)
!q      vvis(i,j)=snu_vp*((yvn(i+1,j)-2.*yvn(i,j)+yvn(i-1,j))/dx2 &
!q       +(yvn(i,j+1)-2.*yvn(i,j)+yvn(i,j-1))/(ddy(j)**2*hs(i)**2))
      end do
    end do
!
    do i = i1 + 1, i2
      do j = 1, ny
        if (j == ny .and. (j_surf == 0 .or. (j_surf == 1 .and. time < stime_surf))) then
          yvn(i, ny) = 0.
        else
          yvn(i, j) = yvn(i, j) + vvis(i, j)*dt
        end if
      end do
    end do
!
!  --- Density -----
!
    if (j_dens == 1) then
      do i = i1, i2 + 1
        do j = 1, ny
          cvis(i, j) = skc*((yc(i + 1, j) - 2.*yc(i, j) + yc(i - 1, j))/dx2 &
                            + (yc(i, j + 1) - 2.*yc(i, j) + yc(i, j - 1))/(dy(j)**2*hs(i)**2))
          ycn(i, j) = yc(i, j) + cvis(i, j)*dt
        end do
      end do
    end if

    if (j_temp == 1) then
      do i = i1, i2 + 1
        do j = 1, ny
          tvis(i, j) = skt*((yt(i + 1, j) - 2.*yt(i, j) + yt(i - 1, j))/dx2 &
                            + (yt(i, j + 1) - 2.*yt(i, j) + yt(i, j - 1))/(dy(j)**2*hs(i)**2))
          ytn(i, j) = yt(i, j) + tvis(i, j)*dt
        end do
      end do
    end if

!
!
  end subroutine diffusion
end module diffusion_m

!==============================================
module dcip0_m
  use common_hh
  use common_hyd
  use common_geom
  implicit none
contains

  subroutine dcip0(f, gx, gy, itp)

    real(8), dimension(0:im, 0:jm)::f, gx, gy
    real(8), dimension(0:im, 0:jm)::fn, gxn, gyn, u, v
    real(8)::hs_up, dx3, a1, e1, b1, f1, d1, c1, g1, tmp, tmq, xx, yy, gxo, gyo
    integer::itp, isn, jsn, im1, jm1, i1, i2, j1, j2

    fn = 0.; gxn = 0.; gyn = 0.; u = 0.; v = 0.
    i1 = 0; i2 = 0; j1 = 0; j2 = 0
!
! itp=1...u, itp=2...v, itp=3...concentration, itp=4...temparature

! uの場合
    if (itp == 1) then
      j1 = 1; j2 = ny
      if (j_bound == 1) then
! 閉鎖水域
        i1 = 2; i2 = nx - 2  !i=1とi=nx-1はゼロ値を境界条件で与える
      else if (j_bound == 2) then
!開放水域
        i1 = 1; i2 = nx - 1 !i=0とi=nxに境界条件を与える(含む周期境界)
      else if (j_bound == 3) then
!上流開放,下流壁
        i1 = 1; i2 = nx - 2 !i=0には境界条件，i=nx-1にゼロを与える
      else
!上流壁,下流開放
        i1 = 2; i2 = nx - 1 !i=1にゼロをi=nxに境界条件を与える
      end if
      do i = i1, i2
        hs_up = (hs(i) + hs(i + 1))*.5
        do j = j1, j2
          u(i, j) = yu(i, j)
          v(i, j) = 0.25*(v2(i, j) + v2(i + 1, j) + v2(i, j - 1) + v2(i + 1, j - 1)) &
                    /hs_up
        end do
      end do
!
!      write(44,*) 'cip itp=1'
!      do i=i1,i2
!       write(44,'(a1,12f12.7)')'u',(u(i,j),j=j1,j2,2)
!       write(44,'(a1,12f12.7)')'v',(v(i,j),j=j1,j2,2)
!      end do

! vの場合
    else if (itp == 2) then
      i1 = 2; i2 = nx - 1
      j1 = 1
      if (j_surf == 1 .and. time > stime_surf) then
        j2 = ny
      else
        j2 = ny - 1
      end if
      do j = j1, j2
        do i = i1, i2
          u(i, j) = 0.25*(yu(i, j) + yu(i - 1, j) + yu(i, j + 1) + yu(i - 1, j + 1))
          v(i, j) = v2(i, j)/hs(i)
        end do
      end do
! c,tの場合
    else if (itp >= 3) then
      j1 = 1; j2 = ny
      i1 = 2; i2 = nx - 1
      do j = j1, j2
        do i = i1, i2
          u(i, j) = 0.5*(yu(i, j) + yu(i - 1, j))
          v(i, j) = 0.5*(v2(i, j) + v2(i, j - 1))/hs(i)
        end do
      end do
    else
      write (*, *) 'dcip-type error'
      stop
    end if
    dx3 = dx2*dx

    do j = j1, j2
      do i = i1, i2
        xx = -u(i, j)*dt
        yy = -v(i, j)*dt
        isn = sign(1.0, u(i, j))
        jsn = sign(1.0, v(i, j))
        im1 = i - isn
        jm1 = j - jsn
        if (itp >= 3) then
!        if(j==j1 .and. jm1==j-1) jm1=j
!        if(j==j2 .and. jm1==j+1) jm1=j
          if ((j_bound == 1 .or. j_bound == 4) .and. i == i1 .and. im1 == i - 1) im1 = i
          if ((j_bound == 1 .or. j_bound == 3) .and. i == i2 .and. im1 == i + 1) im1 = i
          if (obst(i + 1, j) == 1 .and. im1 == i + 1) im1 = i
          if (obst(i - 1, j) == 1 .and. im1 == i - 1) im1 = i
          if (obst(i, j + 1) == 1 .and. jm1 == j + 1) jm1 = j
          if (obst(i, j - 1) == 1 .and. jm1 == j - 1) jm1 = j
        end if
        a1 = ((gx(im1, j) + gx(i, j))*dx*isn - 2.0*(f(i, j) - f(im1, j)))/(dx3*isn)
        e1 = (3.0*(f(im1, j) - f(i, j)) + (gx(im1, j) + 2.0*gx(i, j))*dx*isn)/dx2
        b1 = ((gy(i, jm1) + gy(i, j))*dy(j)*jsn - 2.0*(f(i, j) - f(i, jm1)))/(dy(j)**3*jsn)
        f1 = (3.0*(f(i, jm1) - f(i, j)) + (gy(i, jm1) + 2.0*gy(i, j))*dy(j)*jsn)/dy(j)**2
        tmp = f(i, j) - f(i, jm1) - f(im1, j) + f(im1, jm1)
        tmq = gy(im1, j) - gy(i, j)
        d1 = (-tmp - tmq*dy(j)*jsn)/(dx*dy(j)**2*isn)
        c1 = (-tmp - (gx(i, jm1) - gx(i, j))*dx*isn)/(dx2*dy(j)*jsn)
        g1 = (-tmq + c1*dx2)/(dx*isn)
!--------------------------------------------
        fn(i, j) = ((a1*xx + c1*yy + e1)*xx + g1*yy + gx(i, j))*xx &
                   + ((b1*yy + d1*xx + f1)*yy + gy(i, j))*yy + f(i, j)
        gxn(i, j) = (3.0*a1*xx + 2.0*(c1*yy + e1))*xx + (d1*yy + g1)*yy + gx(i, j)
        gyn(i, j) = (3.0*b1*yy + 2.0*(d1*xx + f1))*yy + (c1*xx + g1)*xx + gy(i, j)
!       write(44,'(2i4,3f12.6)') i,j,fn(i,j),gxn(i,j),gyn(i,j)
      end do
    end do

    do j = j1, j2
      do i = i1, i2

!       if(itp>=3) then
!        if(j==1.or. obst(i,j-1)==1) then
!         if(obst(i,j+1)==1) then
!          gyo=0.
!         else
!          gyo=(f(i,j+1)-f(i,j))/dy(j)
!         end if
!        else if(j==ny.or.obst(i,j+1)==1) then
!         if(obst(i,j-1)==1) then
!          gyo=0.
!         else
!          gyo=(f(i,j)-f(i,j-1))/dy(j)
!         end if
!        else
!         gyo=(f(i,j+1)-f(i,j-1))/(2.*dy(j))
!        end if

!        if((i==2.and.(j_bound==1.or.j_bound==4)) .or. obst(i-1,j)==1) then
!         if(obst(i+1,j)==1) then
!          gxo=0.
!         else
!          gxo=(f(i+1,j)-f(i,j))/dx
!         end if
!        else if((i==nx-1.and.(j_bound==1.or.j_bound==3)).or.obst(i+1,j)==1) then
!         if(obst(i-1,j)==1) then
!          gxo=0.
!         else
!          gxo=(f(i,j)-f(i-1,j))/dx
!         end if
!        else
!         gxo=(f(i+1,j)-f(i-1,j))/(2.*dx)
!        end if
!       else
        gxo = (f(i + 1, j) - f(i - 1, j))/(2.*dx)
        gyo = (f(i, j + 1) - f(i, j - 1))/(2.*dy(j))
!       end if

        f(i, j) = fn(i, j)
        gx(i, j) = gxn(i, j) - (gxo*(u(i + 1, j) - u(i - 1, j)) &
                                + gyo*(v(i + 1, j) - v(i - 1, j)))*0.5*dt/dx
        gy(i, j) = gyn(i, j) - (gxo*(u(i, j + 1) - u(i, j - 1)) &
                                + gyo*(v(i, j + 1) - v(i, j - 1)))*0.5*dt/dy(j)

      end do
    end do

  end subroutine dcip0

!=======================================================

  subroutine upwind(f, gx, gy, itp)

    real(8), dimension(0:im, 0:jm)::f, gx, gy
    real(8), dimension(0:im, 0:jm)::fn, u, v
    real(8)::hs_up, u_dfdx, v_dfdy
    integer::itp, i1, i2, j1, j2
!
    fn = 0.; u = 0.; v = 0.
    hs_up = 0.; u_dfdx = 0.; v_dfdy = 0.
    i1 = 0; i2 = 0; j1 = 0; j2 = 0
!
    if (itp == 1) then
      j1 = 1; j2 = ny
      if (j_bound == 1) then
        i1 = 2; i2 = nx - 2
      else if (j_bound == 2) then
        i1 = 1; i2 = nx - 1
      else if (j_bound == 3) then
        i1 = 1; i2 = nx - 2
      else
        i1 = 2; i2 = nx - 1
      end if
      do i = i1, i2
        hs_up = (hs(i) + hs(i + 1))*.5
        do j = j1, j2
          u(i, j) = yu(i, j)
          v(i, j) = 0.25*(v2(i, j) + v2(i + 1, j) + v2(i, j - 1) + v2(i + 1, j - 1)) &
                    /hs_up
        end do
      end do
    else if (itp == 2) then
      i1 = 2; i2 = nx - 1
      j1 = 1
      if (j_surf == 1 .and. time > stime_surf) then
        j2 = ny
      else
        j2 = ny - 1
      end if
      do j = j1, j2
        do i = i1, i2
          u(i, j) = 0.25*(yu(i, j) + yu(i - 1, j) + yu(i, j + 1) + yu(i - 1, j + 1))
          v(i, j) = v2(i, j)/hs(i)
        end do
      end do
    else if (itp >= 3) then
      j1 = 1; j2 = ny
      i1 = 2; i2 = nx - 1
      do j = j1, j2
        do i = i1, i2
          u(i, j) = 0.5*(yu(i, j) + yu(i - 1, j))
          v(i, j) = 0.5*(v2(i, j) + v2(i, j - 1))/hs(i)
        end do
      end do
    end if
!
    do i = i1, i2
      do j = j1, j2
        u_dfdx = ((u(i, j) + abs(u(i, j)))*(f(i, j) - f(i - 1, j)) &
                  + (u(i, j) - abs(u(i, j)))*(f(i + 1, j) - f(i, j)))/(2.*dx)
        v_dfdy = ((v(i, j) + abs(v(i, j)))*(f(i, j) - f(i, j - 1)) &
                  + (v(i, j) - abs(v(i, j)))*(f(i, j + 1) - f(i, j)))/(2.*dy(j))
        fn(i, j) = f(i, j) - (u_dfdx + v_dfdy)*dt
      end do
    end do
!
    do i = 1, nx
      do j = 1, ny
        f(i, j) = fn(i, j)
        gx(i, j) = (fn(i + 1, j) - fn(i - 1, j))/(2.*dx)
        gy(i, j) = (fn(i, j + 1) - fn(i, j - 1))/(2.*dy(j))
      end do
    end do
!
  end subroutine upwind
end module dcip0_m

!=============================================
module bound_m
  use common_hh
  use common_hyd
  use common_geom
  implicit none
contains

  subroutine bound_u(u, v)

    real(8), dimension(0:im, 0:jm)::u, v
    real(8)::dvds, dhds, hs_up, dvdxi, duds, dd
    integer::i1, i2

!     write(44,*) 'Bound_u(1)'
!     write(44,'(10f10.5)') (u(0,j),j=0,ny+1,2)

    if (j_bound == 1) then
! 閉鎖水域
      i1 = 2; i2 = nx - 2
    else if (j_bound == 2) then
!開放水域
      i1 = 1; i2 = nx - 1
    else if (j_bound == 3) then
!上流開放,下流壁
      i1 = 1; i2 = nx - 2
    else
!上流壁,下流開放
      i1 = 2; i2 = nx - 1
    end if

    do i = i1, i2
      if (j_surf == 0 .or. (j_surf == 1 .and. time < stime_surf)) then
        u(i, ny + 1) = u(i, ny)
      else
        dvds = (v(i + 1, ny) - v(i, ny))/dx
        dhds = (hn(i + 1) - hn(i))/dx
        hs_up = (hs(i + 1) + hs(i))*.5
        dvdxi = (v(i + 1, ny) + v(i, ny) - v(i + 1, ny - 1) - v(i, ny - 1))/(dy(ny)*2.)
        if (i == i1) then
          duds = (u(i + 1, ny) - u(i, ny))/dx
        else if (i == i2) then
          duds = (u(i, ny) - u(i - 1, ny))/dx
        else
          duds = (u(i + 1, ny) - u(i - 1, ny))/(2.*dx)
        end if
        dd = -dvds + dhds*(duds - dvdxi/hs_up)
        u(i, ny + 1) = u(i, ny) + dd*dy(ny)*hs_up
      end if
      u(i, 0) = -u(i, 1)
    end do

!上下流のuの条件
    if (j_bound == 2) then !上下流フリーの場合
      if (j_periodic == 1) then !周期境界条件
        do j = 0, ny + 1
          u(0, j) = u(nx - 2, j)
          u(nx, j) = u(2, j)
        end do
      else  !上下流がフリーの状態
        do j = 0, ny + 1
          if (j_qgive == 1 .and. j_qadj == 0 .and. time > q_stt) then
            u(1, j) = u(1, j)*(1.-q_relax) + u00(j)*q_relax
          end if
          u(0, j) = u(1, j)
          u(nx, j) = u(nx - 1, j)
        end do
      end if
    else if (j_bound == 1) then !上下流閉鎖
      do j = 0, ny + 1
        u(1, j) = 0.
        u(nx - 1, j) = 0.
      end do
    else if (j_bound == 3) then !上流開放・下流は壁
      do j = 0, ny + 1
        u(0, j) = u(1, j)
        u(nx - 1, j) = 0.
      end do
    else if (j_bound == 4) then !下流開放・上流壁
      do j = 0, ny + 1
        u(1, j) = 0.
        u(nx, j) = u(nx - 1, j)
      end do
    end if

    do i = i1, i2
      do j = 1, ny
        if (obst(i, j) == 1 .or. obst(i + 1, j) == 1) u(i, j) = 0.
      end do
    end do

!     write(44,*) 'Bound_u(2)'
!     write(44,'(10f10.5)') (u(0,j),j=0,ny+1,2)

  end subroutine bound_u
!-------------------------------------
  subroutine bound_v(v)

    real(8), dimension(0:im, 0:jm)::v
!水面条件
    do i = 2, nx - 1
      if (j_surf == 1 .and. time > stime_surf) then
        v(i, ny + 1) = 3.*v(i, ny) - 3.*v(i, ny - 1) + v(i, ny - 2)
      else
        v(i, ny + 1) = 0.
      end if
    end do

    if (j_bound == 2) then!開放水域
      if (j_periodic == 1) then !周期境界条件
        do j = 0, ny + 1
          v(1, j) = v(nx - 1, j)
          v(nx, j) = v(2, j)
        end do
      else
        do j = 0, ny + 1
!        v(1,j)=0.
!        v(nx,j)=0.
          v(1, j) = v(2, j)
          v(nx, j) = v(nx - 1, j)
        end do
      end if
    else if (j_bound == 1) then !上下流閉鎖
      do j = 0, ny
!       v(1,j)=0.
!       v(nx,j)=0.
        v(1, j) = v(2, j)
        v(nx, j) = v(nx - 1, j)
      end do
    else if (j_bound == 3) then!上流開放下流壁
      do j = 0, ny
!       v(1,j)=0.
!       v(nx,j)=0.
        v(1, j) = v(2, j)
        v(nx, j) = v(nx - 1, j)
      end do
    else if (j_bound == 4) then!下流開放上流壁
      do j = 0, ny
!       v(1,j)=0.
!       v(nx,j)=0.
        v(1, j) = v(2, j)
        v(nx, j) = v(nx - 1, j)
      end do
    end if

!     if((j_bound==2 .or. j_bound==4).and.j_hdw==2) then
!      do j=0,ny
!       v(nx,j)=-v(nx-1,j)+dhdt_down*xi(j)
!      end do
!     end if

! 底面条件
    do i = 1, nx + 1
      v(i, 0) = 0.
    end do

!     do j=1,ny-1
!      do i=2,nx-1
!       if(obst(i,j)==1 .or. obst(i,j+1)==1) v(i,j)=0.
!      end do
!     end do
    return
  end subroutine bound_v
!--------------------------------------------
  subroutine bound_t(t, j_in_l, j_in_r, j_in_b, j_in_u &
                     , j_l_s, j_l_e, j_r_s, j_r_e, i_b_s, i_b_e, i_u_s, i_u_e &
                     , t_bound_l, t_bound_r, t_bound_b, t_bound_u)
    real(8), dimension(0:im, 0:jm)::t
    integer::m
    integer::j_in_l, j_in_r, j_in_b, j_in_u
    integer::j_l_s(50), j_l_e(50), j_r_s(50), j_r_e(50), &
              i_b_s(50), i_b_e(50), i_u_s(50), i_u_e(50)
    real(8)::t_bound_l(50), t_bound_r(50), t_bound_b(50), t_bound_u(50)

    if (j_bound == 2) then
      if (j_periodic == 1) then!周期境界条件
        do j = 1, ny
          t(1, j) = t(nx - 1, j)
          t(0, j) = t(nx - 2, j)
          t(nx, j) = t(2, j)
          t(nx + 1, j) = t(3, j)
        end do
      else
        do j = 1, ny
          t(1, j) = t(2, j)
          t(nx, j) = t(nx - 1, j)
        end do
      end if

    end if

    if (j_in_l > 0) then
! Left Boundray
      do m = 1, j_in_l
        do j = j_l_s(m), j_l_e(m)
          t(1, j) = t_bound_l(m)
          t(2, j) = t_bound_l(m)
        end do
      end do
    end if

    if (j_in_r > 0) then
! Right Boundray
      do m = 1, j_in_r
        do j = j_r_s(m), j_r_e(m)
          t(nx, j) = t_bound_r(m)
          t(nx - 1, j) = t_bound_r(m)
        end do
      end do
    end if

    if (j_in_b > 0) then
! Bottom Boundray
      do m = 1, j_in_b
        do i = i_b_s(m), i_b_e(m)
          t(i, 0) = t_bound_b(m)
          t(i, 1) = t_bound_b(m)
        end do
      end do
    end if

    if (j_in_u > 0) then
! Upper Boundray
      do m = 1, j_in_u
        do i = i_u_s(m), i_u_e(m)
          t(i, ny + 1) = t_bound_u(m)
          t(i, ny) = t_bound_u(m)
        end do
      end do
    end if

  end subroutine bound_t
!-----------------------------------------------
  subroutine bound_c(c, jc_in_l, jc_in_r, jc_in_b, jc_in_u &
                     , jc_l_s, jc_l_e, jc_r_s, jc_r_e, ic_b_s, ic_b_e, ic_u_s, ic_u_e &
                     , c_bound_l, c_bound_r, c_bound_b, c_bound_u)

    real(8), dimension(0:im, 0:jm)::c
    integer::m
    integer::jc_in_l, jc_in_r, jc_in_b, jc_in_u
    integer::jc_l_s(50), jc_l_e(50), jc_r_s(50), jc_r_e(50), &
              ic_b_s(50), ic_b_e(50), ic_u_s(50), ic_u_e(50)
    real(8)::c_bound_l(50), c_bound_r(50), c_bound_b(50), c_bound_u(50)

    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 1, ny
        c(1, j) = c(nx - 1, j)
        c(0, j) = c(nx - 2, j)
        c(nx, j) = c(2, j)
        c(nx + 1, j) = c(3, j)
      end do
    end if

    if ((j_bound == 2 .and. j_periodic == 0) .or. (j_bound == 3)) then !上流開放
      do j = 1, ny
        if (yu(1, j) > 0.) then
          c(1, j) = con0
        else
          c(1, j) = con0
        end if
      end do
    else
      do j = 1, ny
        c(1, j) = con0
      end do
    end if
    if ((j_bound == 2 .and. j_periodic == 0) .or. (j_bound == 4)) then !下流開放
      do j = 1, ny
        if (yu(nx - 1, j) < 0.) then
          c(nx, j) = con0
        else
          c(nx, j) = con0
        end if
      end do
    else
      do j = 1, ny
        c(nx, j) = con0
      end do
    end if

    if (jc_in_l > 0) then
! Left Boundray
      do m = 1, jc_in_l
        do j = jc_l_s(m), jc_l_e(m)
          c(1, j) = c_bound_l(m)
          if (j_bound == 2 .or. j_bound == 3) c(2, j) = c_bound_l(m)
        end do
      end do
    end if

    if (jc_in_r > 0) then
! Right Boundray
      do m = 1, jc_in_r
        do j = jc_r_s(m), jc_r_e(m)
          c(nx, j) = c_bound_r(m)
          if (j_bound == 2 .or. j_bound == 4) c(nx - 1, j) = c_bound_r(m)
        end do
      end do
    end if

    if (jc_in_b > 0) then
! Bottom Boundray
      do m = 1, jc_in_b
        do i = ic_b_s(m), ic_b_e(m)
          c(i, 0) = c_bound_b(m)
!        c(i,1)=c_bound_b(m)
        end do
      end do
    end if

    if (jc_in_u > 0) then
! Upper Boundray
      do m = 1, jc_in_u
        do i = ic_u_s(m), ic_u_e(m)
          c(i, ny + 1) = c_bound_u(m)
!        c(i,ny)=c_bound_u(m)
        end do
      end do
    end if

  end subroutine bound_c
end module bound_m
!=================================================
module gbound_m
  use common_hh
  use common_hyd
  implicit none
contains
!-----------------------------------------------
  subroutine gbound_u(ynx, yny)
    real(8), dimension(0:im, 0:jm)::ynx, yny

    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 1, ny
        ynx(0, j) = ynx(nx - 2, j)
        ynx(nx, j) = ynx(2, j)
        yny(0, j) = yny(nx - 2, j)
        yny(nx, j) = yny(2, j)
      end do
    else
      do j = 0, ny + 1
        ynx(1, j) = 0.
        ynx(nx - 1, j) = 0.
        yny(1, j) = 0.
        yny(nx - 1, j) = 0.
        ynx(0, j) = 0.
        ynx(nx, j) = 0.
        yny(0, j) = 0.
        yny(nx, j) = 0.
      end do
    end if

    do i = 0, nx
!p     ynx(i,0)=ynx(i,1)
!p     ynx(i,ny+1)=ynx(i,ny)
!p     yny(i,0)=yny(i,1)
!p     yny(i,ny+1)=yny(i,ny)
      ynx(i, 0) = 0.
      ynx(i, ny + 1) = 0.
      yny(i, 0) = 0.
      yny(i, ny + 1) = 0.
    end do
    return
  end subroutine gbound_u
!--------------------------------------------
  subroutine gbound_v(ynx, yny)
    real(8), dimension(0:im, 0:jm)::ynx, yny

    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 1, ny
        ynx(1, j) = ynx(nx - 1, j)
        ynx(nx, j) = ynx(2, j)
        yny(1, j) = yny(nx - 1, j)
        yny(nx, j) = yny(2, j)
      end do
    else
      do j = 1, ny
        ynx(1, j) = 0.
        ynx(nx, j) = 0.
        yny(1, j) = 0.
        yny(nx, j) = 0.
      end do
    end if

    do i = 1, nx
!p     ynx(i,0)=ynx(i,1)
!p     yny(i,0)=yny(i,1)
      ynx(i, 0) = 0.
      ynx(i, ny) = 0.
      yny(i, 0) = 0.
      yny(i, ny) = 0.
    end do
  end subroutine gbound_v
!--------------------------------------------
  subroutine gbound_t(ynx, yny)
    real(8), dimension(0:im, 0:jm)::ynx, yny

    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 1, ny
        ynx(1, j) = ynx(nx - 1, j)
        ynx(nx, j) = ynx(2, j)
        yny(1, j) = yny(nx - 1, j)
        yny(nx, j) = yny(2, j)
      end do
    else
      do j = 1, ny
        ynx(1, j) = 0.
        ynx(nx, j) = 0.
        yny(1, j) = 0.
        yny(nx, j) = 0.
      end do
    end if

!r     do i=1,nx
!r      ynx(i,0)=0.
!r      ynx(i,ny)=0.
!r      yny(i,0)=0.
!r      yny(i,ny)=0.
!r     end do

  end subroutine gbound_t

end module gbound_m

!================================================
module newgrd_m
  use common_hh
  use common_geom
  use common_hyd
  use common_grad
  use gbound_m
  implicit none
contains
!------------------------------------------------
  subroutine newgrd(yn, y, gxn, gyn, gx, gy, itp)

    real(8), dimension(0:im, 0:jm)::yn, y, gxn, gyn, gx, gy
    integer::itp, i1, i2, j1, j2

    i1 = 0; i2 = 0; j1 = 0; j2 = 0
!
    if (j_bound == 1) then
! 閉鎖水域
      i1 = 2; i2 = nx - 2  !i=1とi=nx-1はゼロ値を境界条件で与える
    else if (j_bound == 2) then
!開放水域
      i1 = 1; i2 = nx - 1 !i=0とi=nxに境界条件を与える(含む周期境界)
    else if (j_bound == 3) then
!上流開放,下流壁
      i1 = 1; i2 = nx - 2 !i=0には境界条件，i=nx-1にゼロを与える
    else
!上流壁,下流開放
      i1 = 2; i2 = nx - 1 !i=1にゼロをi=nxに境界条件を与える
    end if

!
    if (itp == 1) then
      j1 = 1
      j2 = ny - 1
      do i = i1, i2
        do j = j1, j2
          gxn(i, j) = gx(i, j) + (yn(i + 1, j) - yn(i - 1, j) - y(i + 1, j) + y(i - 1, j))*0.5/dx
          gyn(i, j) = gy(i, j) + (yn(i, j + 1) - yn(i, j - 1) - y(i, j + 1) + y(i, j - 1))*0.5/dy(j)
        end do
      end do
    else if (itp == 2) then
      j1 = 2
      j2 = ny - 1
      do i = i1, i2 + 1
        do j = j1, j2
          gxn(i, j) = gx(i, j) + (yn(i + 1, j) - yn(i - 1, j) - y(i + 1, j) + y(i - 1, j))*0.5/dx
          gyn(i, j) = gy(i, j) + (yn(i, j + 1) - yn(i, j - 1) - y(i, j + 1) + y(i, j - 1))*0.5/ddy(j)
        end do
      end do
    else if (itp >= 3) then
      j1 = 2
      j2 = ny - 1
      do i = i1, i2 + 1
        do j = j1, j2
          if (i == i1) then
            gxn(i, j) = gx(i, j) + (yn(i + 1, j) - yn(i, j) - y(i + 1, j) + y(i, j))/dx
          else if (i == i2 + 1) then
            gxn(i, j) = gx(i, j) + (yn(i, j) - yn(i - 1, j) - y(i, j) + y(i - 1, j))/dx
          else
            gxn(i, j) = gx(i, j) + (yn(i + 1, j) - yn(i - 1, j) - y(i + 1, j) + y(i - 1, j))*0.5/dx
          end if
          if (j == j1) then
            gyn(i, j) = gy(i, j) + (yn(i, j + 1) - yn(i, j) - y(i, j + 1) + y(i, j))/dy(j)
          else if (j == j2) then
            gyn(i, j) = gy(i, j) + (yn(i, j) - yn(i, j - 1) - y(i, j) + y(i, j - 1))/dy(j)
          else
            gyn(i, j) = gy(i, j) + (yn(i, j + 1) - yn(i, j - 1) - y(i, j + 1) + y(i, j - 1))*0.5/dy(j)
          end if
        end do
      end do
    end if
    if (itp == 1) then
      call gbound_u(gxn, gyn)
    else if (itp == 2) then
      call gbound_v(gxn, gyn)
    else
      call gbound_t(gxn, gyn)
    end if
  end subroutine newgrd
end module newgrd_m

!=============================================

module output_m
  use common_hh
  use common_geom
  use iric
  implicit none
contains

  subroutine output(fid, u, v, p, t, c, h, eta_up, qc, xg, yg_bed, obst)
    real(8), dimension(0:im, 0:jm)::u, v, p, t, c
    integer, dimension(0:im, 0:jm)::obst
    real(8), dimension(0:im)::h, eta_up, h_up, hs_up, qc
    real(8), dimension(0:im, 0:jm)::ux, vx, px, tx, cx
    real(8), dimension(ni, nj)::xg, yg, ug, vg, tg, cg, pg, vorg, hsg, hg, qcg
    real(8), dimension(ni)::yg_bed
    real(8), dimension(ni - 1, nj - 1)::t_cell, c_cell, p_cell
    real(8)::dudy, dvdx
    integer::fid, ier
    integer::num = 0
    real(8)::sum_p, sum_t, sum_c

    ug = 0.; vg = 0.; pg = 0.; hsg = 0.; vorg = 0.; yg = 0.; hg = 0.

    do i = 1, nx - 1
      h_up(i) = (h(i) + h(i + 1))*.5
    end do
    do i = 1, nx - 1
      hs_up(i) = h_up(i) - eta_up(i)
      do j = 0, ny
        yg(i, j + 1) = yg_bed(i) + hs_up(i)*xi(j)
        ux(i, j) = (u(i, j) + u(i, j + 1))*.5
        vx(i, j) = (v(i, j) + v(i + 1, j))*.5
        dudy = (u(i, j + 1) - u(i, j))/(dy(j)*hs_up(i))
        dvdx = (v(i + 1, j) - v(i, j))/dx
        num = 0; sum_p = 0.; sum_t = 0.; sum_c = 0.
        if ((obst(i, j) == 0) .and. (i > 1) .and. (j > 0)) then
          num = num + 1
          sum_p = sum_p + p(i, j); sum_t = sum_t + t(i, j); sum_c = sum_c + c(i, j)
        else if ((obst(i + 1, j) == 0) .and. (i < nx - 1) .and. (j > 0)) then
          num = num + 1
          sum_p = sum_p + p(i + 1, j); sum_t = sum_t + t(i + 1, j); sum_c = sum_c + c(i + 1, j)
        else if ((obst(i, j + 1) == 0) .and. (i > 1) .and. (j < ny)) then
          num = num + 1
          sum_p = sum_p + p(i, j + 1); sum_t = sum_t + t(i, j + 1); sum_c = sum_c + c(i, j + 1)
        else if ((obst(i + 1, j + 1) == 0) .and. (i < nx - 1) .and. (j < ny)) then
          num = num + 1
          sum_p = sum_p + p(i + 1, j + 1); sum_t = sum_t + t(i + 1, j + 1); sum_c = sum_c + c(i + 1, j + 1)
        end if
        if (num > 0) then
          px(i, j) = sum_p/float(num)
          tx(i, j) = sum_t/float(num)
          cx(i, j) = sum_c/float(num)
        else
          px(i, j) = 0.
          tx(i, j) = tmp0
          cx(i, j) = con0
        end if
        ug(i, j + 1) = ux(i, j); vg(i, j + 1) = vx(i, j)
        tg(i, j + 1) = tx(i, j); cg(i, j + 1) = cx(i, j); pg(i, j + 1) = px(i, j)
        vorg(i, j + 1) = dudy - dvdx
        hsg(i, j + 1) = hs_up(i)
        hg(i, j + 1) = h_up(i)
        qcg(i, j + 1) = qc(i)
      end do
    end do
    do i = 1, ni - 1
      do j = 1, nj - 1
        t_cell(i, j) = t(i + 1, j)
        c_cell(i, j) = c(i + 1, j)
        p_cell(i, j) = p(i + 1, j)
      end do
    end do
!     write(44,*)
!     write(44,'(a3,30e12.3)') 'qc ',(qc(i),i=0,nx)
!     write(44,'(a3,30e12.3)') 'qcg',(qcg(i,1),i=1,ni)

!     write(133) time
!     write(133) (hs_up(i),i=1,nx-1)
!     write(133) ((ux(i,j),i=1,nx-1),j=0,ny)
!     write(133) ((vx(i,j),i=1,nx-1),j=0,ny)
!     write(133) ((px(i,j),i=1,nx-1),j=0,ny)
!     write(133) ((tx(i,j),i=1,nx-1),j=0,ny)
!     write(133) ((cx(i,j),i=1,nx-1),j=0,ny)

    call cg_iric_write_sol_start(fid, ier)
    call cg_iric_write_sol_time(fid, time, ier)
    call cg_iric_write_sol_baseiterative_real(fid, "Input Discharge", qp0, ier)
    call cg_iric_write_sol_baseiterative_real(fid, "Averageed Discharge", q_ave, ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Upstream Water Surface Elevation ", h_up_boundary, ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Downstream Water Surface Elevation ", h_down, ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Upstream calculated discharge", qc(1), ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Downstream calculated discharge", qc(nx - 1), ier)
    call cg_iric_write_sol_baseiterative_real(fid, "Volume", t_vol, ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Average Depth", hs_ave_c, ier)
    call cg_iric_write_sol_baseiterative_real &
      (fid, "Average Velocity", u_ave_c, ier)
    call cg_iRIC_Write_Sol_Grid2d_Coords(fid, xg, yg, ier)
    call cg_iRIC_Write_Sol_Node_Real(fid, "VelocityX", ug, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "VelocityY", vg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Pressure", pg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Temperature", tg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Concentrartion", cg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Vorticity", vorg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Depth", hsg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Water surface elevation", hg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "Dischare per unit width", qcg, IER)
    call cg_iRIC_Write_Sol_Node_Real(fid, "X-Velocity", ug, IER)
    call cg_iric_write_sol_cell_real(fid, "P_cell", p_cell, ier)
    call cg_iric_write_sol_cell_real(fid, "T_cell", t_cell, ier)
    call cg_iric_write_sol_cell_real(fid, "C_cell", c_cell, ier)
    call cg_iric_write_sol_end(fid, ier)
  end subroutine output
end module output_m

!=================================================
module shift_m
  use common_hh
  implicit none
contains
!--------------------------------------------
  subroutine shift(yn, y)
    real(8), dimension(0:im, 0:jm)::yn, y

    do j = 0, ny + 1
      do i = 0, nx + 1
        y(i, j) = yn(i, j)
      end do
    end do

  end subroutine shift
!--------------------------------------------
  subroutine hshift(h, hn, hs, eta)
    real(8), dimension(0:im)::h, hn, hs, eta
    do i = 1, nx
      h(i) = hn(i)
      hs(i) = hn(i) - eta(i)
    end do
  end subroutine hshift
end module shift_m
!=================================================

module v12cal_m
  use common_hh
  use common_hyd
  use common_geom
  use ss_nu_m
  implicit none
contains
!--------------------------------------------
  subroutine v1cal(u, v)
    real(8), dimension(0:im, 0:jm)::u, v
    real(8)::hs_up, dhsds, deds, u_vp, ome
    integer::i1, i2

    if (j_bound == 1) then !閉鎖水域
      i1 = 2; i2 = nx - 1
    else if (j_bound == 2) then !開放水域
      i1 = 1; i2 = nx
    else if (j_bound == 3) then !上流開放,下流壁
      i1 = 1; i2 = nx - 1
    else ! 上流壁,下流開放
      i1 = 2; i2 = nx
    end if
    do i = i1, i2 - 1
      hs_up = (hs(i) + hs(i + 1))*.5
      u(i, ny + 1) = u(i, ny) - (v(i + 1, ny) - v(i, ny))*dy(ny)*hs_up/dx
    end do

!     write(44,*) 'v1cal'
    do i = i1, i2
      dhsds = (hs(i + 1) - hs(i - 1))/(2.*dx)
      deds = (eta(i + 1) - eta(i - 1))/(2.*dx)
      do j = 0, ny
        u_vp = (u(i, j) + u(i, j + 1) + u(i - 1, j) + u(i - 1, j + 1))*.25
        ome = xi(j)*dhsds + deds
        v1(i, j) = v(i, j) - ome*u_vp
!       if((i<5 .or. i>nx-5).and.(j==1 .or. j>ny-5)) then
!        write(44,'(2i5,4f12.5)') i,j,u_vp,dhsds,deds,v1(i,j)
!        write(44,'(4f12.5)') u(i,j),u(i,j+1),u(i-1,j),u(i-1,j+1)
!       end if
      end do
    end do
    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 0, ny
        v1(1, j) = v1(nx - 1, j)
        v1(nx, j) = v1(2, j)
      end do
    end if
  end subroutine v1cal
!--------------------------------------------
  subroutine v2cal
    integer::i1, i2

    if (j_bound == 1) then !閉鎖水域
      i1 = 2; i2 = nx - 1
    else if (j_bound == 2) then!開放水域
      i1 = 1; i2 = nx
    else if (j_bound == 3) then!上流開放,下流壁
      i1 = 1; i2 = nx - 1
    else ! 上流壁,下流開放
      i1 = 2; i2 = nx
    end if
!     write(44,*) 'v2cal',i1,i2

    do i = i1, i2
      do j = 0, ny
!       write(44,'(2i4,3f10.6)') i,j,v1(i,j),xi(j),v1(i,ny)
        v2(i, j) = v1(i, j) - xi(j)*v1(i, ny)
      end do
    end do
    if (j_bound == 2) then
      do j = 0, ny
        v2(1, j) = v2(nx - 1, j)
        v2(nx, j) = v2(2, j)
      end do
    end if
  end subroutine v2cal
!--------------------------------------------
  subroutine qcal(u)
    real(8), dimension(0:im, 0:jm)::u
    real(8)::hs_up, us_hp
    integer::i1, i2
    integer::nq

    if (j_bound == 1) then
! 閉鎖水域
      i1 = 2; i2 = nx - 2  !i=1とi=nx-1はゼロ値を境界条件で与える
    else if (j_bound == 2) then
!開放水域
      i1 = 1; i2 = nx - 1 !i=0とi=nxに境界条件を与える(含む周期境界)
    else if (j_bound == 3) then
!上流開放,下流壁
      i1 = 1; i2 = nx - 2 !i=0には境界条件，i=nx-1にゼロを与える
    else
!上流壁,下流開放
      i1 = 2; i2 = nx - 1 !i=1にゼロをi=nxに境界条件を与える
    end if

    q_ave = 0.; nq = 0; t_vol = 0.; u_ave_c = 0.; hs_ave_c = 0.; t_length = 0.
    do i = i1, i2
      qc(i) = 0.
      hs_up = (hs(i) + hs(i + 1))*.5
      do j = 1, ny
        qc(i) = qc(i) + u(i, j)*hs_up*ddy(j)
      end do
      uba(i) = qc(i)/hs_up
      cd1 = 1./(1./kappa*log(xi1*hs_up/z0))
      cd(i) = cd1**2
      usta(i) = u(i, 1)*cd1
      q_ave = q_ave + qc(i)
      u_ave_c = u_ave_c + uba(i)
      nq = nq + 1
      if (i > i1) then
        t_vol = t_vol + hs(i)*dx
        t_length = t_length + dx
      end if
      do j = 0, ny
        snu_tx(i, j) = ss_nu(usta(i), hs_up, xi(j))
      end do
    end do
    q_ave = q_ave/float(nq)
    hs_ave_c = t_vol/t_length
    u_ave_c = u_ave_c/float(nq)

    if (j_bound == 2) then
      if (j_periodic == 1) then
        do i = 0, 1
          qc(i) = qc(nx - 2 + i)
          usta(i) = usta(nx - 2 + i)
        end do
        do i = nx, nx + 1
          qc(i) = qc(i - nx + 2)
          usta(i) = usta(i - nx + 2)
        end do
      else
        qc(0) = qc(1)
        usta(0) = usta(1)
        qc(nx) = qc(nx - 1)
        usta(nx) = usta(nx - 1)
      end if
    end if
!     write(44,'(a4,30e12.3)') 'qcal',(qc(i),i=0,nx)
!
    do i = i1 + 1, i2
      us_hp = (usta(i) + usta(i - 1))*.5
      do j = 1, ny
        snu_t(i, j) = ss_nu(us_hp, hs(i), xxi(j))
      end do
    end do

    if (j_bound == 2 .and. j_periodic == 1) then
      do j = 1, ny
        snu_t(1, j) = snu_t(nx - 1, j)
        snu_t(nx, j) = snu_t(2, j)
      end do
    else
      do j = 1, ny
        snu_t(1, j) = snu_t(2, j)
        snu_t(nx, j) = snu_t(nx - 1, j)
      end do
    end if
  end subroutine qcal
!--------------------------------------------
  subroutine hcal
    real(8)::dhdt

    do i = 2, nx - 1
      dhdt = v1(i, ny)
      hn(i) = h(i) + dhdt*dt*alpha_surf
      hs(i) = hn(i) - eta(i)
    end do
!
! ----- 上下流下流水位境界条件 -----
!
    if (j_bound == 2 .and. j_periodic == 1) then !周期境界条件
      hs(1) = hs(nx - 1)
      hn(1) = eta(1) + hs(1)
      hs(nx) = hs(2)
      hn(nx) = eta(nx) + hs(nx)
    end if

    if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 4) then !下流端開放
      if (j_hdw == 3) then
        hn(nx) = hn(nx - 1)
      else
        hn(nx) = h_down
      end if
      hs(nx) = hn(nx) - eta(nx)
    else  !下流壁
      hs(nx) = hs(nx - 1)
      hn(nx) = eta(nx) + hs(nx)
    end if
    h_down = hn(nx)
    hs(nx) = hn(nx) - eta(nx)

    if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 3) then !上流開放
      if (j_qgive == 0) then
        if (j_hup == 3) then
          hn(1) = hn(2)
        else
          hn(1) = h_up_boundary
        end if
        hs(1) = hn(1) - eta(1)
      else if (j_qgive == 1) then
        if (j_qadj == 0) then! discharge adjust by velocity
          hn(1) = hn(2)
        else !discharge adjust by water surface elevation at upstream end
          hn(1) = hn(2) + h_plus
          h_up_boundary = hn(1)
        end if
      else  !上流壁
        hs(1) = hs(2)
        hn(1) = eta(1) + hs(1)
      end if
      h_up_boundary = hn(1)
      hs(1) = hn(1) - eta(1)
    end if

    if (j_bound == 1 .or. j_bound == 4) then
      hn(1) = hn(2)
      hs(1) = hn(1) - eta(1)
    end if
    if (j_bound == 1 .or. j_bound == 3) then
      hn(nx) = hn(nx - 1)
      hs(nx) = hn(nx) - eta(nx)
    end if

  end subroutine hcal
!--------------------------------------------
end module v12cal_m
!============================================

!*************************************************
!
! Main
!
!*************************************************
Program Shimizu
  use common_hh
  use common_geom
  use common_hyd
  use common_grad
  use alloc_var_m
  use initial_val_m
  use bound_p_m
  use sor_m
  use rhs_m
  use diffusion_m
  use newgrd_m
  use dcip0_m
  use bound_m
  use gbound_m
  use output_m
  use shift_m
  use v12cal_m
  use iric

  implicit none

  integer, parameter :: strMax = 1000
  character(len=strMax) :: condFile

  real(8)::chl, h0, tuk, etime, sorerr, soralpha, surf_tension
  real(8)::smg_g
  integer::lsor, icount, itout, itt
  integer::m, hloop, hloop0, lmax
  real(8)::hloop_err, hs_up0, u_ave0, zz, sb

!------ Values for iRIC ------
  integer::ircount, ier, fid, istatus, fin
!------ for total error check -----
  real(8)::diff_ratio
!------ Values for iRIC CGS Setting -----
  real(8), dimension(:, :), allocatable::xg, yg, zg
  real(8), dimension(:), allocatable::yg_bed
  real(8), dimension(:, :), allocatable::tg, cg
  integer, dimension(:, :), allocatable::obs_g
  real(8)::total_1, total_2
!
! ====== Values for Bounday Temperature(CGNS I/O) ======
  integer::j_in, indexmax
  integer, dimension(:), allocatable :: j_inlen, j_size
  integer, dimension(:, :, :), allocatable:: indices
  double precision, dimension(:), allocatable :: boundary_tmp_value
  character(250), dimension(:), allocatable :: flowname
!-----------------------------------------
  integer::j_in_l, j_in_r, j_in_b, j_in_u
  integer::j_l_s(50), j_l_e(50), j_r_s(50), j_r_e(50), &
            i_b_s(50), i_b_e(50), i_u_s(50), i_u_e(50)
  real(8)::t_bound_l(50), t_bound_r(50), t_bound_b(50), t_bound_u(50)
!
! ======= Values for Bounday Concentration(CGNS I/O) ======
  integer::jc_in, indexmax_c
  integer, dimension(:), allocatable :: jc_inlen, jc_size
  integer, dimension(:, :, :), allocatable:: indices_c
  double precision, dimension(:), allocatable :: boundary_con_value
  character(250), dimension(:), allocatable :: flowname_c
!-----------------------------------------
  integer::jc_in_l, jc_in_r, jc_in_b, jc_in_u
  integer::jc_l_s(50), jc_l_e(50), jc_r_s(50), jc_r_e(50), &
            ic_b_s(50), ic_b_e(50), ic_u_s(50), ic_u_e(50)
  real(8)::c_bound_l(50), c_bound_r(50), c_bound_b(50), &
            c_bound_u(50)

!------ Background Condition -----
  integer::j_bcv
  integer::j_bottom_shape
  real(8)::dune_height
!-------- Values for DSWSE----------------------------
  integer::j_set_hdw
  real(8)::h_down_given, h_down_ini, amp_alpha &
            , hd_amp, hd_wl, hd_st, hd_ap &
            , hu_amp, hu_wl, hu_st, hu_ap
  real(8)::pi
!-------- Values for USWSE----------------------------
  integer::j_set_hup
  real(8)::h_up_given, h_up_ini

!------------------------------------------
!     open(44,file='c:\tmp\tmp.d',status='unknown')
!     open(133,file='c:\tmp\out.d',status='unknown',form='unformatted')
!------------------------------------------

  ircount = 0; ier = 0; fid = 0; ni = 0; nj = 0
  tuk = 0.; etime = 0.; dt = 0.
  sorerr = 0.; lsor = 0; soralpha = 0.
  snu = 0.; skt = 0.; skc = 0.; beta_t = 0.; rho = 0.; surf_tension = 0.
  hloop = 0; hloop0 = 0; hloop_err = 0.; m = 0
  diff_ratio = 0.; h0 = 0.
! ----- T boundary values -----
  j_in_l = 0; j_in_r = 0; j_in_b = 0; j_in_u = 0
  j_l_s = 0; j_l_e = 0; j_r_s = 0; j_r_e = 0
  i_b_s = 0; i_b_e = 0; i_u_s = 0; i_u_e = 0
  t_bound_l = 0.; t_bound_r = 0.; t_bound_b = 0.; t_bound_u = 0.
! ------ C boundary values ------
  jc_in_l = 0; jc_in_r = 0; jc_in_b = 0; jc_in_u = 0
  jc_l_s = 0; jc_l_e = 0; jc_r_s = 0; jc_r_e = 0
  ic_b_s = 0; ic_b_e = 0; ic_u_s = 0; ic_u_e = 0
  c_bound_l = 0.; c_bound_r = 0.; c_bound_b = 0.; c_bound_u = 0.
! -------------------------------------------
  j_bcv = 0
  total_1 = 0.; total_2 = 0.
  dhdt_down = 0.; h_down_old = 0.

  lmax = 0

  pi = 4.*atan(1.0d0)

!==========================================================
! 格子生成データファイルを開く
  call cg_iric_open("./nays2dv/iRICZone/gridcreate.cgn", IRIC_MODE_MODIFY, fin, ier)
  if (ier /= 0) stop "*** Open error of CGNS file ***"

!--- Read from grid generator file -----
!  格子生成条件データファイルからの読み出し
  call cg_iric_read_real(fin, 'h0', hs_ave, ier)
  call cg_iric_read_real(fin, 'b_slope', slope, ier)
  call cg_iric_read_integer(fin, 'j_bottom_shape', j_bottom_shape, ier)
  call cg_iric_read_real(fin, 'dune_height', dune_height, ier)
!     if(j_bottom_shape==2) hs_ave=hs_ave+dune_height/2.

!     write(*,'(a8,3f10.5)') 'h,slope=',hs_ave,slope,dune_height

  call cg_iric_close(fin, ier)
!=========================================================

  ircount = nargs()
  if (ircount == 2) then
    call getarg(1, condFile, ier)
  else
    write (*, *) "Input File not specified."
    stop
  end if

! 計算データファイルを開く
  call cg_iric_open(condFile, IRIC_MODE_MODIFY, fid, ier)
  if (ier /= 0) then
    write (*, *) "*** Open error of CGNS file ***"
  end if

!
!guiにcgnsファイルを読込みであることを知らせるファイルを生成
  call iric_initoption(IRIC_OPTION_CANCEL, ier)

!
! 2次元格子データーのサイズを読み込む ni,nj
!
  call cg_iRIC_Read_Grid2d_Str_Size(fid, ni, nj, ier)
!     write(*,*) 'ni,nj=',ni,nj

  if (ier /= 0) then
    write (*, *) "*** error:cg_iRIC_Read_Grid2d_Str_Size ***"
    stop
  end if
!
  allocate (xg(ni, nj), yg(ni, nj), zg(ni, nj), yg_bed(ni))
  allocate (tg(ni - 1, nj - 1), cg(ni - 1, nj - 1), obs_g(ni - 1, nj - 1))
  xg = 0.; yg = 0.; zg = 0.; tg = 0.; cg = 0.; yg_bed = 0.

! 2次元格子データーを読み込む xg(ni,nj),yg(ni,nj)
  call cg_iRIC_Read_Grid2d_Coords(fid, xg, yg, ier)

!     write(44,*) 'Grid Data'
!     do i=1,ni
!      write(44,'(i3,30e12.3)') i,(xg(i,j),j=1,nj)
!      write(44,'(i3,30e12.3)') i,(yg(i,j),j=1,nj)
!     end do
!
  if (ier /= 0) then
    write (*, *) "*** error:cg_iRIC_Read_Grid2d_Coords_f ***"
    stop
  end if
! 仮想・河床データの読み込み
  call cg_iric_read_grid_real_node(fid, "Elevation", zg, ier)
  if (ier /= 0) then
    write (*, *) "error:cg_iric_read_grid_real_node_f-Elevation "
    stop
  end if

! 初期密度・温度データの読み込み
  call cg_iric_read_grid_real_cell(fid, 'Ini_tmp', tg, ier)
  if (ier /= 0) stop
  call cg_iric_read_grid_real_cell(fid, 'Ini_con', cg, ier)
  if (ier /= 0) stop

! 障害物セルデータの読み込み
  call cg_iric_read_grid_integer_cell(fid, 'Obstacle', obs_g, ier)
  if (ier /= 0) stop

!     do i=1,ni-1
!      do j=1,nj-1
!       write(44,'(2i5,2f10.4)') i,j,tg(i,j),cg(i,j)
!      end do
!     end do

  call cg_iric_read_integer(fid, 'j_bcv', j_bcv, ier)
!  j_ccv 0=Min, 1=Max, 2=Average, 3=Give
  if (j_bcv == 3) then
    call cg_iric_read_real(fid, 'tmp0', tmp0, ier)
    call cg_iric_read_real(fid, 'con0', con0, ier)
  end if

! 背景温度・背景濃度tmp0,con0の計算

  if (j_bcv < 3) then
    if (j_bcv == 0) then
      tmp0 = 9999.
      con0 = 9999.
    else if (j_bcv == 1) then
      tmp0 = -9999.
      con0 = -9999.
    else if (j_bcv == 2) then
      tmp0 = 0.
      con0 = 0.
    end if

    do i = 1, ni - 1
      do j = 1, nj - 1
        if (j_bcv == 0) then
          tmp0 = min(tmp0, tg(i, j))
          con0 = min(con0, cg(i, j))
        else if (j_bcv == 1) then
          tmp0 = max(tmp0, tg(i, j))
          con0 = max(con0, cg(i, j))
        else if (j_bcv == 2) then
          tmp0 = tmp0 + tg(i, j)
          con0 = con0 + cg(i, j)
        end if
      end do
    end do
    if (j_bcv == 2) then
      tmp0 = tmp0/float((ni - 1)*(nj - 1))
      con0 = con0/float((ni - 1)*(nj - 1))
    end if
  end if

!     write(*,*) 'tmp0,con0=',tmp0,con0
!     write(44,*) 'tmp0,con0=',tmp0,con0

!     do i=1,ni
!      do j=1,nj,nj-1
!       write(44,'(2i4,2e12.3)') i,j,xg(i,j),yg(i,j)
!      end do
!     end do
!     write(44,*)
!     do i=1,ni-1
!      write(44,'(21f3.0)') (tg(i,j),j=1,nj-1)
!     end do
!     write(44,*)
!     do i=1,ni-1
!      write(44,'(21f3.0)') (cg(i,j)*100.,j=1,nj-1)
!     end do

!--- Read Computational Parameters  -----

  call cg_iric_read_integer(fid, 'j_dens', j_dens0, ier)
  call cg_iric_read_real(fid, 'st_dens', st_dens, ier)
!
  call cg_iric_read_integer(fid, 'j_temp', j_temp0, ier)
  call cg_iric_read_real(fid, 'st_temp', st_temp, ier)

  call cg_iric_read_integer(fid, 'j_uadvec', j_uadvec, ier)
! j_uadvec=1; CIP Scheme
! j_uadvec=2; Upwind Scheme
  call cg_iric_read_integer(fid, 'j_cadvec', j_cadvec, ier)
! j_cadvec=1; Density CIP Scheme
! j_cadvec=2; Density Upwind Scheme
  call cg_iric_read_integer(fid, 'j_bound', j_bound, ier)
! j_bound=1; closed
! j_bound=2; open
! j_bound=3; upstream open and downstream closed
! j_bound=4; downstream open and upstream closed

  call cg_iric_read_integer(fid, 'j_periodic', j_periodic, ier)
! j_periodic=0; non periodic open channel
! j_periodic=1; periodic bundary

  qp = 0.; qp0 = 0.; q_stt = 0.; q_trn = 0.
  if (j_bound == 2) then
    call cg_iric_read_integer(fid, 'j_qgive', j_qgive, ier)
    if (j_qgive == 1) then
      call cg_iric_read_real(fid, 'qp', qp, ier)
      call cg_iric_read_integer(fid, 'j_qadj', j_qadj, ier)
      call cg_iric_read_real(fid, 'q_relax', q_relax, ier)
      call cg_iric_read_real(fid, 'q_stt', q_stt, ier)
      call cg_iric_read_real(fid, 'q_trn', q_trn, ier)
    end if
  end if

  if (j_bound == 1 .or. j_bound == 3) j_hdw = 3
  if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 4) then
!下流水位の条件
    call cg_iric_read_integer(fid, 'j_hdw', j_hdw, ier)
! 下流端水位条件 1==一定, 2==Sine波形, 3==自由(dH/dx=0)
    if (j_hdw == 1) then
      call cg_iric_read_integer(fid, 'j_set_hdw', j_set_hdw, ier)
! 一定値の条件 0==初期条件と同じ, 1==値の入力
      if (j_set_hdw == 1) then
        call cg_iric_read_real(fid, 'h_down_given', h_down_given, ier)
! 下流端水位の値
      end if
    else if (j_hdw == 2) then
! Sine波形の
      call cg_iric_read_real(fid, 'hd_amp', hd_amp, ier) !波高
      call cg_iric_read_real(fid, 'hd_wl', hd_wl, ier) !波長
      call cg_iric_read_real(fid, 'hd_st', hd_st, ier) !開始時刻
      call cg_iric_read_real(fid, 'hd_ap', hd_ap, ier) !開始時刻
    end if
  end if
  if (j_bound == 1 .or. j_bound == 4) j_hup = 3
  if (j_qgive == 1) then !上流から流量を与える場合
    j_hup = 3
    j_set_hup = 0
    hu_amp = 0.; hu_wl = 0.; hu_st = 0.; hu_ap = 0.
  else if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 3) then
!上流水位の条件
    call cg_iric_read_integer(fid, 'j_hup', j_hup, ier)
!上流水位の条件　1=一定， 2=Sine波, 3=自由(dH/dx=0)
    if (j_hup == 1) then
      call cg_iric_read_integer(fid, 'j_set_hup', j_set_hup, ier)
! 一定値の条件 0==初期条件と同じ, 1==値の入力
      if (j_set_hup == 1) then
        call cg_iric_read_real(fid, 'h_up_given', h_up_given, ier)
! 上流端水位の値
      end if
    else if (j_hup == 2) then
! Sine波形の
      call cg_iric_read_real(fid, 'hu_amp', hu_amp, ier) !波高
      call cg_iric_read_real(fid, 'hu_wl', hu_wl, ier) !波長
      call cg_iric_read_real(fid, 'hu_st', hu_st, ier) !開始時刻
      call cg_iric_read_real(fid, 'hu_ap', hu_ap, ier) !開始時刻
    end if
  end if

  call cg_iric_read_real(fid, 'diam', diam, ier)
  call cg_iric_read_integer(fid, 'j_snu', j_snu, ier)
  call cg_iric_read_real(fid, 'al_ep', al_ep, ier)

  call cg_iric_read_real(fid, 'tuk', tuk, ier)
  call cg_iric_read_real(fid, 'etime', etime, ier)
  call cg_iric_read_real(fid, 'dt', dt, ier)
!     write(*,'(3f10.5)') tuk,etime,dt
  call cg_iric_read_real(fid, 'sorerr', sorerr, ier)
  call cg_iric_read_integer(fid, 'lsor', lsor, ier)
! lsor=sorの反復回数(5回程度) ----------
  call cg_iric_read_real(fid, 'soralpha', soralpha, ier)
! soralpha = sor法の緩和回数(=1.5で加速緩和) ------
!     write(*,'(f10.5,i5,f10.5)') sorerr,lsor,soralpha

  call cg_iric_read_integer(fid, 'j_surf', j_surf, ier)
! 自由水面計算 j_surf=1 (Yes)  j_suef=0 (no)
! j_surf=0--->固定水面
  call cg_iric_read_real(fid, 'alpha_surf', alpha_surf, ier)
  call cg_iric_read_real(fid, 'stime_surf', stime_surf, ier)

!     write(*,*) tuk,etime,dt,sorerr,lsor,soralpha,j_surf

  call cg_iric_read_integer(fid, 'hloop', hloop, ier)
! 自由水面計算 繰り返し回数
  if (j_surf == 0) hloop = 1

  call cg_iric_read_real(fid, 'hloop_err', hloop_err, ier)
! 自由水面計算 繰り返し計算打ち切り誤差

  call cg_iric_read_real(fid, 'snu', snu, ier)
  call cg_iric_read_real(fid, 'skt', skt, ier)
  call cg_iric_read_real(fid, 'skc', skc, ier)
  call cg_iric_read_real(fid, 'beta_t', beta_t, ier)
  call cg_iric_read_real(fid, 'rho', rho, ier)
  call cg_iric_read_real(fid, 'surf_tension', surf_tension, ier)

!     write(*,'(7e12.3)') snu,skt,skc,beta_t,rho,surf_tension
! ---------------------------------------------------
! 境界条件の読み込み
! Read Boundray Condition of Temperature
! ---------------------------------------------------
!     write(44,*) 'ni,nj=',ni,nj
! 温度 j_in 境界条件の個数
  call cg_iric_read_bc_count(fid, 'b_tmp', j_in)
!     write(44,*) 'j_in=',j_in
! ---------------------------------------------------
! 濃度 j_in_c 境界条件の個数
  call cg_iric_read_bc_count(fid, 'b_con', jc_in)
!     write(44,*) 'jc_in=',jc_in
! ---------------------------------------------------

! -----------温度境界条件の読み込み -----

  allocate (j_inlen(j_in))  ! ede number*2
  allocate (boundary_tmp_value(j_in))  !Boundary Temperature Value
  allocate (flowname(j_in)) ! Boundary Name

  indexmax = 0 ! max of j_inlen(i)

  do i = 1, j_in
    call cg_iric_read_bc_indicessize(fid, 'b_tmp', i, j_inlen(i), ier)
!      write(44,*) 'i,j_inlen(i)=',i,j_inlen(i)    !境界要素の数
    if (indexmax < j_inlen(i)) indexmax = j_inlen(i)
    call cg_iric_read_bc_string(fid, 'b_tmp', i, '_caption', flowname(i), ier)
    !境界条件要素の名前
!      write(44,*) flowname(i)(1:30)
  end do

!     write(44,*) 'indexmax=',indexmax
  allocate (indices(j_in, 2, indexmax))

  do i = 1, j_in
    call cg_iric_read_bc_indices(fid, 'b_tmp', i, indices(i:i, :, :), ier)
    call cg_iric_read_bc_real(fid, 'b_tmp', i, 'b_tmp_val', &
                              boundary_tmp_value(i), ier)
!      write(44,*) 'tmp=',boundary_tmp_value(i)
!      do j=1,2
!       do k=1,j_inlen(i)
!        write(44,*) 'i,j,k,indices=',i,j,k,indices(i,j,k)
!       end do
!      end do
  end do
!
!-------境界条件の場所を探す-----
! Assign Boudary Condition
!
!     write(44,*) j_in,j_in_l,j_in_r,j_in_b
!     do i=1,j_in
!      do j=1,2
!       do k=1,j_inlen(i)
!        write(44,*) i,j,k,indices(i,j,k)
!       end do
!      end do
!     end do

  do i = 1, j_in
    do j = 1, j_inlen(i) - 1
! Check i=1 Boundary (Left Boundary)
      if (indices(i, 1, 1) == 1 .and. indices(i, 1, j_inlen(i)) == 1) then
        if (j == 1) then
          j_in_l = j_in_l + 1
          t_bound_l(j_in_l) = boundary_tmp_value(i)
!         write(44,*) 'Bound L',j_in_l,t_bound_l(j_in_l)
          j_l_s(j_in_l) = indices(i, 2, 1)
          j_l_e(j_in_l) = indices(i, 2, j_inlen(i) - 1)
!         write(44,*) 'Bound L start end=',j_l_s(j_in_l),j_l_e(j_in_l)
        end if
! Check i=ni Boundary (Right Boundary)
      else if (indices(i, 1, 1) == ni .and. indices(i, 1, j_inlen(i)) == ni) then
        if (j == 1) then
          j_in_r = j_in_r + 1
          t_bound_r(j_in_r) = boundary_tmp_value(i)
!         write(44,*) 'Bound R',j_in_r,t_bound_r(j_in_r)
          j_r_s(j_in_r) = indices(i, 2, 1)
          j_r_e(j_in_r) = indices(i, 2, j_inlen(i) - 1)
!         write(44,*) 'Bound R start end=',j_r_s(j_in_r),j_r_e(j_in_r)
        end if
! Check j=1 Boundary (Bottom Boundary)
      else if (indices(i, 2, j) == 1) then
        if (j == 1) then
          j_in_b = j_in_b + 1
          t_bound_b(j_in_b) = boundary_tmp_value(i)
!         write(44,*) 'Bound B',j_in_b,t_bound_b(j_in_b)
          i_b_s(j_in_b) = indices(i, 1, 1) + 1
          i_b_e(j_in_b) = indices(i, 1, j_inlen(i) - 1) + 1
!         write(44,*) 'Bound B start end=',i_b_s(j_in_b),i_b_e(j_in_b)
        end if
! Check j=nj Boundary (Upper Boundary)
      else if (indices(i, 2, j) == nj) then
        if (j == 1) then
          j_in_u = j_in_u + 1
          t_bound_u(j_in_u) = boundary_tmp_value(i)
!         write(44,*) 'Bound U',j_in_u,t_bound_u(j_in_u)
          i_u_s(j_in_u) = indices(i, 1, 1) + 1
          i_u_e(j_in_u) = indices(i, 1, j_inlen(i) - 1) + 1
!         write(44,*) 'Bound U start end=',i_u_s(j_in_u),i_u_e(j_in_u)
        end if
      end if
    end do
  end do

! -----------濃度境界条件の読み込み -----

  allocate (jc_inlen(jc_in))  ! ede number*2
  allocate (boundary_con_value(jc_in))  !Boundary Concentration Value
  allocate (flowname_c(jc_in)) ! Boundary Name

  flowname_c = ''

  indexmax_c = 0 ! max of jc_inlen(i)

  do i = 1, jc_in
    call cg_iric_read_bc_indicessize(fid, 'b_con', i, jc_inlen(i), ier)
!      write(44,*) 'i,jc_inlen(i)=',i,jc_inlen(i)
    if (indexmax_c < jc_inlen(i)) indexmax_c = jc_inlen(i)
    call cg_iric_read_bc_string(fid, 'b_con', i, '_caption', flowname_c(i), ier)
!      write(44,*) flowname_c(i)(1:50)
  end do

!     write(44,*) 'indexmax_c=',indexmax_c
  allocate (indices_c(jc_in, 2, indexmax_c))

  do i = 1, jc_in
    call cg_iric_read_bc_indices(fid, 'b_con', i, indices_c(i:i, :, :), ier)
    call cg_iric_read_bc_real(fid, 'b_con', i, 'b_con_val', &
                              boundary_con_value(i), ier)
!      write(44,*) 'concentration=',boundary_con_value(i)
!      do j=1,2
!       do k=1,jc_inlen(i)
!        write(44,*) 'i,j,k,indices_c=',i,j,k,indices_c(i,j,k)
!       end do
!      end do
  end do
!
!-------境界条件の場所を探す-----
! Assign Boudary Condition
!
!     write(44,*) 'jc_in,jc_in_l,jc_in_r,jc_in_b=', &
!                   jc_in,jc_in_l,jc_in_r,jc_in_b
!     do i=1,jc_in
!      do j=1,2
!       do k=1,jc_inlen(i)
!        write(44,*) i,j,k,indices_c(i,j,k)
!       end do
!      end do
!     end do

  do i = 1, jc_in
    do j = 1, jc_inlen(i) - 1
! Check i=1 Boundary (Left Boundary)
      if (indices_c(i, 1, 1) == 1 .and. indices_c(i, 1, jc_inlen(i)) == 1) then
        if (j == 1) then
          jc_in_l = jc_in_l + 1
          c_bound_l(jc_in_l) = boundary_con_value(i)
!         write(44,*) '[con]Bound L',jc_in_l,c_bound_l(jc_in_l)
          jc_l_s(jc_in_l) = indices_c(i, 2, 1)
          jc_l_e(jc_in_l) = indices_c(i, 2, jc_inlen(i) - 1)
!         write(44,*) '[con]Bound L start end=',jc_l_s(jc_in_l),jc_l_e(jc_in_l)
        end if
! Check i=ni Boundary (Right Boundary)
      else if (indices_c(i, 1, 1) == ni .and. indices_c(i, 1, jc_inlen(i)) == ni) then
        if (j == 1) then
          jc_in_r = jc_in_r + 1
          c_bound_r(jc_in_r) = boundary_con_value(i)
!         write(44,*) '[con]Bound R',jc_in_r,c_bound_r(jc_in_r)
          jc_r_s(jc_in_r) = indices_c(i, 2, 1)
          jc_r_e(jc_in_r) = indices_c(i, 2, jc_inlen(i) - 1)
!         write(44,*) '[con]Bound R start end=',jc_r_s(jc_in_r),jc_r_e(jc_in_r)
        end if
! Check j=1 Boundary (Bottom Boundary)
      else if (indices_c(i, 2, j) == 1) then
        if (j == 1) then
          jc_in_b = jc_in_b + 1
          c_bound_b(jc_in_b) = boundary_con_value(i)
!         write(44,*) '[con]Bound B',jc_in_b,c_bound_b(jc_in_b)
          ic_b_s(jc_in_b) = indices_c(i, 1, 1) + 1
          ic_b_e(jc_in_b) = indices_c(i, 1, jc_inlen(i) - 1) + 1
!         write(44,*) '[con]Bound B start end=',ic_b_s(jc_in_b),ic_b_e(jc_in_b)
        end if
! Check j=nj Boundary (Upper Boundary)
      else if (indices_c(i, 2, j) == nj) then
        if (j == 1) then
          jc_in_u = jc_in_u + 1
          c_bound_u(jc_in_u) = boundary_con_value(i)
!         write(44,*) '[con]Bound U',jc_in_u,c_bound_u(jc_in_u)
          ic_u_s(jc_in_u) = indices_c(i, 1, 1) + 1
          ic_u_e(jc_in_u) = indices_c(i, 1, jc_inlen(i) - 1) + 1
!         write(44,*) '[con]Bound U start end=',ic_u_s(jc_in_u),ic_u_e(jc_in_u)
        end if
      end if
    end do
  end do
!     write(44,*) jc_in_l
!     write(44,*) jc_in_r
!     write(44,*) jc_in_b
!     write(44,*) jc_in_u
  nx = ni + 1
  ny = nj - 1
  chl = xg(ni, 1) - xg(1, 1)

!------------------------------------------------

  im = nx + 5
  jm = ny + 5

  call alloc_var(im, jm)
  call initial_0

!------------------------------------------------

! yg(i,j)--->yg_bed(i)
  do i = 1, ni
    yg_bed(i) = yg(i, 1)
  end do

!  tg(i,j)=ytn(i+1,j) --->  yt(i,j)=tg(i-1,j)

  do i = 1, nx - 1
    qc(i) = 0.
  end do

!
!     write(44,*) 'yg'
!     do i=1,ni
!      write(44,'(i3,2e12.3)') i,yg(i,1),yg(i,nj)
!     end do

  do i = 2, nx - 1
    h(i) = (yg(i - 1, nj) + yg(i, nj))*.5
    eta(i) = (yg(i - 1, 1) + yg(i, 1))*.5
  end do

!     write(44,*)'h,eta'
!     do i=2,nx-1
!      write(44,'(i3,2e12.3)') i,h(i),eta(i)
!     end do

  if (j_bound == 2 .or. j_bound == 3) then
!      eta(1)=eta(2)+(eta(2)-eta(3))
!      eta(0)=eta(1)+(eta(1)-eta(2))
!      h(1)=h(2)+(h(2)-h(3))                  !上流端水位の初期値
!      h(0)=h(1)+(h(1)-h(2))                  !上流端水位の初期値
    eta(1) = eta(2)
    eta(0) = eta(1)
    h(1) = h(2)
    h(0) = h(1)
  else
    eta(1) = eta(2)
    eta(0) = eta(1)
    h(1) = h(2)
    h(0) = h(1)
  end if
  if (j_bound == 2 .or. j_bound == 4) then
!      eta(nx)=eta(nx-1)-(eta(nx-2)-eta(nx-1))
!      eta(nx+1)=eta(nx)-(eta(nx-1)-eta(nx))
!      h(nx)=h(nx-1)-(h(nx-2)-h(nx-1))        !下流端水位の初期値
!      h(nx+1)=h(nx)-(h(nx-1)-h(nx))        !下流端水位の初期値
    eta(nx) = eta(nx - 1)
    eta(nx + 1) = eta(nx)
    h(nx) = h(nx - 1)
    h(nx + 1) = h(nx)
  else
    eta(nx) = eta(nx - 1)
    eta(nx + 1) = eta(nx)
    h(nx) = h(nx - 1)
    h(nx + 1) = h(nx)
  end if
  h_up_ini = h(1)
  h_down_ini = h(nx)
!------------------------------------------------
  do i = 0, nx + 1
    hs(i) = h(i) - eta(i)    !初期水深
    hn(i) = h(i)
  end do

!     do i=0,3
!      write(44,'(i3,5e12.3)') i,eta(i),hs(i),h(i),hn(i),eta(i)+hs(i)
!     end do

  hs_ave = 0.
  do i = 1, nx
    hs_ave = hs_ave + hs(i)
  end do
  hs_ave = hs_ave/float(nx)  ! hs_ave==格子データから計算した平均水深の初期値
  ! h0==格子データから読んできた平均水深
  u_ave = 0.

  do i = 1, nx - 1
    eta_up(i) = yg(i, 1)
  end do
!------------------------------------------------

  do i = 2, nx - 1
    do j = 1, ny
!       write(*,'(2i4)') i,j
      yt(i, j) = tg(i - 1, j)
      yc(i, j) = cg(i - 1, j)
      obst(i, j) = obs_g(i - 1, j)
    end do
!      write(44,'(50i1)')(obst(i,j),j=1,ny)
    yt(i, ny + 1) = tmp0; yc(i, ny + 1) = con0
    yt(i, 0) = tmp0; yc(i, 0) = con0
  end do
  do j = 0, ny + 1
    yt(1, j) = tmp0; yt(nx, j) = tmp0
    yc(1, j) = con0; yc(nx, j) = con0
  end do

  do i = 1, nx
    obst(i, 0) = 1
  end do
  if (j_bound == 1 .or. j_bound == 3) then
    do j = 0, ny
      obst(nx, j) = 1
    end do
  end if
  if (j_bound == 1 .or. j_bound == 4) then
    do j = 0, ny
      obst(1, j) = 1
    end do
  end if

!     do i=1,nx
!      do j=0,ny+1
!       ytn(i,j)=yt(i,j); ycn(i,j)=yc(i,j)
!       write(44,'(2i4,2e12.3)') i,j,yt(i,j),yc(i,j)
!      end do
!     end do

  dx = chl/float(nx - 2)
  itout = int(tuk/dt)
  time = 0.
  itt = itout
  g = 9.8
  kappa = 0.4
  dx2 = dx*dx
!
  smg_g = -surf_tension/(rho*g)
!     write(*,*) 'smg_g=',smg_g

!
! ------ Values for Velocity and Shear Cal. ------
! 流れの計算用パラメータの設定
!     write(*,*) 'j_bound=',j_bound
! ---------------------------------------------------
  ks = 2.*diam
  z0 = ks/30.

  if (j_qgive == 1) then
    if (abs(qp) > 1e-6) then
      q_ref = qp
    else
      q_ref = hs_ave*0.1
    end if
  end if

  call initial

!     write(44,'(a23,4f12.6)') 'u_ave,us_ave,qp,hs_ave=' &
!        ,u_ave,us_ave,qp,hs_ave

  icount = 0
  call qcal(yu)
  t_vol_ini = t_vol

2000 continue
!
! Main Roop Starting Point
!--------------------------------------------------
  icount = icount + 1

  call cg_iric_check_update(fid, ier)

!
! 密度流(濃度)計算フラグのON/Off
!
  if (j_dens0 == 0 .or. time < st_dens) then
    j_dens = 0
  else
    j_dens = 1
  end if
!
! 密度流(温度)計算フラグのON/Off
!
  if (j_temp0 == 0 .or. time < st_temp) then
    j_temp = 0
  else
    j_temp = 1
  end if
!
! ---- Setup Upstream and Downstream Water Surface Elevation -----
! 上下流境界条件水位の設定
!
  if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 4) then
    if (j_hdw == 1) then
      if (j_set_hdw == 0) then
        h_down = h_down_ini
      else
        h_down = h_down_given
      end if
    else if (j_hdw == 2) then
      if (time < hd_st) then
        h_down = h_down_ini
      else
        if (time < (hd_st + hd_ap)) then
          amp_alpha = 0.5*(1.-cos(pi/hd_ap*(time - hd_st)))
        else
          amp_alpha = 1.0
        end if
        h_down = h_down_ini + hd_amp*sin(2.*pi*(time - hd_st)/hd_wl)*amp_alpha
      end if
    end if
  end if

  if (icount == 1) then
    h_down_old = h_down
    dhdt_down = 0.
  else
    dhdt_down = (h_down - h_down_old)/dt
    h_down_old = h_down
  end if

  if ((j_bound == 2 .and. j_periodic == 0) .or. j_bound == 3) then
    if (j_hup == 1) then
      if (j_set_hup == 0) then
        h_up_boundary = h_up_ini
      else
        h_up_boundary = h_up_given
      end if
    else if (j_hup == 2) then
      if (time < hu_st) then
        h_up_boundary = h_up_ini
      else
        if (time < (hu_st + hu_ap)) then
          amp_alpha = 0.5*(1.-cos(pi/hu_ap*(time - hu_st)))
        else
          amp_alpha = 1.0
        end if
        h_up_boundary = h_up_ini + hu_amp*sin(2.*pi*(time - hu_st)/hu_wl)*amp_alpha
      end if
    end if
  end if
!
! ----- 上流端流量の計算 -----
!
  if (j_bound == 2 .and. j_qgive == 1) then

! 目標流量の設定(qp0)
    if (time <= q_stt) then
      qp0 = 0.
    else if (time <= (q_stt + q_trn)) then
      qp0 = qp*0.5*(1 - cos(pi/q_trn*(time - q_stt)))
    else if (time > (q_stt + q_trn)) then
      qp0 = qp
    end if

! 目標の平均流速
!
    if (j_qadj == 0) then
!       hs_up0=(hs(1)+hs(2))*.5
      hs_up0 = hs(1)
      sb = log(hs_up0/z0) - 1.
      u_ave0 = qp0/hs_up0
      qp1 = 0.
      do j = 1, ny
        zz = xxi(j)*hs_up0
        u00(j) = u_ave0*log(zz/z0)/sb
        qp1 = qp1 + u00(j)*hs_up0*ddy(j)
      end do
      u00(0) = -u00(1); u00(ny + 1) = u00(ny)
    else
      h_plus = (qp0 - qc(1))/q_ref*hs(1)*q_relax
    end if
  end if

!
! ---- Nan Check ------
!
  do i = 1, nx
    do j = 0, ny
      if (isnan(yu(i, j)) .or. isnan(yp(i, j))) then
        write (*, *) 'Calculation is Falure (Nan found)', time
        call cg_iric_close(fid, ier)
        stop
      end if
    end do
30  if (isnan(hs(i)) .or. isnan(qc(i)) .or. isnan(usta(i))) then
      write (*, *) 'Calculation is Falure (Nan found)', time
      call cg_iric_close(fid, ier)
      stop
    end if
  end do

  call iric_check_cancel(istatus)
  if (istatus == 1) then
    write (*, *) &
      "Solver is stopped because the STOP button was clicked."
    call cg_iric_close(fid, ier)
    stop
  end if

!--------------------------------------------------
  if (itt == itout) then
    itt = 0
    call output(fid, yu, yv, yp, yt, yc, h, eta_up, qc, xg, yg_bed, obst)
    write (*, '(a5,7f10.5)') &
      'time=', time, qp, qp0, qp1, q_ave, h_down, h_up_boundary
    write (44, '(a5,7f10.5)') &
      'time=', time, qp, qp0, qp1, q_ave, h_down, h_up_boundary
!      do j=ny,1,-1
!       write(44,'(30f10.3)') (yp(i,j),i=1,nx,2)
!      end do
!      write(44,'(5f12.8)') h(nx-1),h(nx),h_down,hs(nx-1),hs(nx)
!      write(44,*) '(1)yu'
!      write(44,'(10f10.5)') (yu(0,j),j=0,ny+1,2)
!      write(44,*) '(2)yun'
!      write(44,'(10f10.5)') (yun(0,j),j=0,ny+1,2)
  end if
  itt = itt + 1
!--------------------------------------------------
!     yuo=yu
!     yvo=yv

  if (j_surf == 0 .or. (j_surf == 1 .and. time < stime_surf)) then
    hloop0 = 1
  else
    hloop0 = hloop
  end if
  do m = 1, hloop0
!      call v1cal(yu,yv)    !<---dune の計算はこっち, 密度流もこっちみたい
    call v1cal(yun, yvn)  !<---これにすると自由水面セイシュOK
!-------------------------------------------------
    if (j_surf == 1 .and. time > stime_surf) then
      call hcal                          ! h-->hn,hs
      call bound_p(smg_g)                ! hn,ypn-->ypn(i,ny+1)
    end if
!--------------------------------------------------
    call omg
    call sor(sorerr, lsor, soralpha, lmax, &
             jc_in_l, jc_in_r, jc_in_b, jc_in_u &
             , jc_l_s, jc_l_e, jc_r_s, jc_r_e, ic_b_s, ic_b_e, ic_u_s, ic_u_e &
             , c_bound_l, c_bound_r, c_bound_b, c_bound_u)

!
!ypn,yun,v1-->ypn
!      if(j_bound==2 .and. j_periodic==1) call pbound
!--------------------------------------------------

    call rhs  !yu,yv,ypn-->yun,yvn
    call bound_u(yun, yvn)  !bound_u(1)
    call bound_v(yvn)
!--------------------------------------------------
    call diffusion !yun,yvn-->yun,yvn  yt,yc-->ytn,ycn
    call bound_u(yun, yvn)  !bound_u(2)
    call bound_v(yvn)
    call newgrd(yun, yu, gux_n, guy_n, gux, guy, 1)
    if (j_temp == 1) then
      call bound_t(ytn, j_in_l, j_in_r, j_in_b, j_in_u &
                   , j_l_s, j_l_e, j_r_s, j_r_e, i_b_s, i_b_e, i_u_s, i_u_e &
                   , t_bound_l, t_bound_r, t_bound_b, t_bound_u)
      call newgrd(ytn, yt, gtx_n, gty_n, gtx, gty, 3)
    end if
    if (j_dens == 1) then
      call bound_c(ycn, jc_in_l, jc_in_r, jc_in_b, jc_in_u &
                   , jc_l_s, jc_l_e, jc_r_s, jc_r_e, ic_b_s, ic_b_e, ic_u_s, ic_u_e &
                   , c_bound_l, c_bound_r, c_bound_b, c_bound_u)
      call newgrd(ycn, yc, gcx_n, gcy_n, gcx, gcy, 4)
    end if
!---------------------------------------------------
    call v1cal(yun, yvn)
    call v2cal
!-------------------------------------------------
    if (j_uadvec == 1) then
      call dcip0(yun, gux_n, guy_n, 1) !yun,v2--->yun
    else if (j_uadvec == 2) then
      call upwind(yun, gux_n, guy_n, 1) !yun,v2--->yun
    end if
    call bound_u(yun, yvn)   !bound_u(3)
    call gbound_u(gux_n, guy_n)
!-------------------------------------------------
    if (j_uadvec == 1) then
      call dcip0(yvn, gvx_n, gvy_n, 2)
    else
      call upwind(yvn, gvx_n, gvy_n, 2)
    end if
    call bound_v(yvn)
    call gbound_v(gvx_n, gvy_n)
!-------------------------------------------------
    if (j_temp == 1) then
      if (j_cadvec == 1) then
        call dcip0(ytn, gtx_n, gty_n, 3)
      else
        call upwind(ytn, gtx_n, gty_n, 3)
      end if
      call bound_t(ytn, j_in_l, j_in_r, j_in_b, j_in_u &
                   , j_l_s, j_l_e, j_r_s, j_r_e, i_b_s, i_b_e, i_u_s, i_u_e &
                   , t_bound_l, t_bound_r, t_bound_b, t_bound_u)
!       call gbound_t(gtx_n,gty_n)
    end if
!-------------------------------------------------
    if (j_dens == 1) then
      if (j_cadvec == 1) then
        call dcip0(ycn, gcx_n, gcy_n, 4)
      else
        call upwind(ycn, gcx_n, gcy_n, 4)
      end if
      call bound_c(ycn, jc_in_l, jc_in_r, jc_in_b, jc_in_u &
                   , jc_l_s, jc_l_e, jc_r_s, jc_r_e, ic_b_s, ic_b_e, ic_u_s, ic_u_e &
                   , c_bound_l, c_bound_r, c_bound_b, c_bound_u)
!       call gbound_t(gcx_n,gcy_n)
    end if
!--------------------------------------------------
    total_1 = 0.; total_2 = 0.
    do i = 1, nx
      do j = 1, ny + 1
        total_1 = total_1 + abs(yun(i, j) - yuo(i, j))
        total_2 = total_2 + abs(yun(i, j))
        yuo(i, j) = yun(i, j)
        yvo(i, j) = yvn(i, j)
      end do
    end do
    if (total_2 == 0.) then
      diff_ratio = 0.
    else
      diff_ratio = total_1/total_2
    end if
    if (diff_ratio < hloop_err) exit
  end do
!--------------------------------------------------
  call shift(yun, yu)
  call shift(yvn, yv)
  call shift(ypn, yp)
  call shift(gux_n, gux)
  call shift(guy_n, guy)
  call shift(gvx_n, gvx)
  call shift(gvy_n, gvy)
  if (j_temp == 1) then
    call shift(ytn, yt)
    call shift(gtx_n, gtx)
    call shift(gty_n, gty)
  end if
  if (j_dens == 1) then
    call shift(ycn, yc)
    call shift(gcx_n, gcx)
    call shift(gcy_n, gcy)
  end if
  call qcal(yun)
!--------------------------------------------------
  if (j_surf == 1 .and. time > stime_surf) then
    call hshift(h, hn, hs, eta)
  end if
  if (time > etime) goto 3000
  time = time + dt
  goto 2000
3000 continue
!     close(133)
end Program
