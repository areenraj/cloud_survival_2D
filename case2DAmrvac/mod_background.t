! Custom fluid background representing the ISM around falling cloud
module mod_background

  use mod_global_parameters
  implicit none

  ! Using a lookup table to speed up the process of allocating values at each time-step
  integer, parameter          :: nz_tab    = 50001
  double precision, parameter :: z_min_tab = 0.0d0     ! in kpc
  double precision, parameter :: z_max_tab = 100.0d0   ! in kpc

  ! Arrays to store lookup table values
  double precision :: z_tab(nz_tab)
  double precision :: rho_tab(nz_tab)
  double precision :: g_tab(nz_tab)
  double precision :: p_tab(nz_tab)

  ! Miller & Bregman flattened model
  double precision, parameter :: mb_n0   = 0.46d0    ! cm^-3
  double precision, parameter :: mb_Rc   = 0.42d0    ! kpc
  double precision, parameter :: mb_zc   = 0.26d0    ! kpc
  double precision, parameter :: mb_beta = 0.71d0

  ! Galactocentric radius for the simulation column (kpc)
  ! R_Sun is the Ferriere reference radius — fixed.
  double precision            :: R0       = 8.122d0   ! overridable from par file
  double precision, parameter :: R_Sun    = 8.5d0
  double precision            :: p_offset = 5.4d-13
  ! Ferriere–Miller halo blending regime
  double precision, parameter :: z_blend_lo = 1.0d0   ! kpc
  double precision, parameter :: z_blend_hi = 2.0d0   ! kpc

  ! Gravitational potential parameters - Bulge
  double precision, parameter :: M_b = 2.54d0      ! 10^10 M_sun
  double precision, parameter :: a_b  = 0.425d0    ! kpc

  ! Gravitional potential parameters - Halo
  double precision, parameter :: r_h  = 5.56d0     ! kpc
  double precision, parameter :: v_h  = 169.77d0   ! km/s

  double precision, parameter :: fmass_disk = 1.0d0
 
  ! Disk parameters for Miyamoto-Nagai Disks (Mass in 10^10 Msun and distance in kpc)

  ! Thin stellar disk — Model 3
  double precision, parameter :: b_thin  = 0.243d0     
  double precision, parameter :: M_thin1 =  2.282d0     
  double precision, parameter :: M_thin2 =  2.342d0
  double precision, parameter :: M_thin3 = -1.846d0
  double precision, parameter :: a_thin1 =  3.859d0    
  double precision, parameter :: a_thin2 =  9.052d0
  double precision, parameter :: a_thin3 =  3.107d0
 
  ! Thick stellar disk — Model 1
  double precision, parameter :: b_thick  = 0.776d0
  double precision, parameter :: M_thick1 =  0.061d0
  double precision, parameter :: M_thick2 =  4.080d0
  double precision, parameter :: M_thick3 = -3.521d0
  double precision, parameter :: a_thick1 =  0.993d0
  double precision, parameter :: a_thick2 =  6.555d0
  double precision, parameter :: a_thick3 =  7.651d0
 
  ! H I disk — Model 2
  double precision, parameter :: b_HI  = 0.168d0
  double precision, parameter :: M_HI1 =  2.217d0
  double precision, parameter :: M_HI2 =  2.350d0
  double precision, parameter :: M_HI3 = -3.303d0
  double precision, parameter :: a_HI1 =  9.021d0
  double precision, parameter :: a_HI2 =  9.143d0
  double precision, parameter :: a_HI3 =  7.758d0
 
  ! H2 disk — Model 3
  double precision, parameter :: b_H2  = 0.128d0
  double precision, parameter :: M_H21 =  1.005d0
  double precision, parameter :: M_H22 =  0.177d0
  double precision, parameter :: M_H23 = -0.907d0
  double precision, parameter :: a_H21 =  6.062d0
  double precision, parameter :: a_H22 =  3.141d0
  double precision, parameter :: a_H23 =  4.485d0

  ! step size for numerical solution of gravitational potential
  double precision, parameter :: dz_potential = 1.0d-3     ! 1 pc

  ! Par file namelist
  namelist /bg_list/ R0

contains

  subroutine bg_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, bg_list, end=111)
111    close(unitpar)
    end do
  end subroutine bg_params_read

  ! Build look-up tables
  subroutine init_background()
    integer          :: i
    double precision :: zk, dz

    do i = 1, nz_tab
       z_tab(i) = z_min_tab + (z_max_tab - z_min_tab) * dble(i-1) / dble(nz_tab-1)
    end do

    do i = 1, nz_tab
       zk = z_tab(i)
       call gas_density_cgs(zk, R0, rho_tab(i))
       call gravity_cgs(zk, R0, g_tab(i))
    end do

    p_tab(nz_tab) = p_offset
    do i = nz_tab-1, 1, -1
       dz = (z_tab(i+1) - z_tab(i)) * 1.0d3 * const_pc
       p_tab(i) = p_tab(i+1) &
                + 0.5d0 * (rho_tab(i)*g_tab(i) + rho_tab(i+1)*g_tab(i+1)) * dz
    end do

    do i = 1, nz_tab
       rho_tab(i) = rho_tab(i) / unit_density
       g_tab(i)   = g_tab(i)   / (unit_length / unit_time**2)
       p_tab(i)   = p_tab(i)   / unit_pressure
    end do

  end subroutine init_background


  ! Ferriere density profile - Molecular Hydrogen + CNM + WNM + WIM + HM

  subroutine ferriere(z, R, nH)
    double precision, intent(in)  :: z, R
    double precision, intent(out) :: nH
    double precision :: aR, H1, H2, H3, Hm, Hh
    double precision :: n_mol, n_cnm, n_wnm, n_wim, n_hm

    ! Radial Scaling Factor
    if (R > R_Sun) then
       aR = R / R_Sun
    else
       aR = 1.0d0
    end if

    ! Profile Parameters
    H1 = 0.127d0 * aR
    H2 = 0.318d0 * aR
    H3 = 0.403d0 * aR
    Hm = 0.081d0 * (R/R_Sun)**0.58d0

    ! Molecular 
    n_mol = 0.58d0 &
          * exp( -((R-4.5d0)**2 - (R_Sun-4.5d0)**2) / 2.9d0**2 ) &
          * (R/R_Sun)**(-0.58d0) &
          * exp( -(z/Hm)**2 )

    ! Cold neutral 
    n_cnm = (0.340d0 / aR**2) * ( &
            0.859d0 * exp( -(z/H1)**2 ) + &
            0.047d0 * exp( -(z/H2)**2 ) + &
            0.094d0 * exp( -z/H3 ) )

    ! Warm neutral 
    n_wnm = (0.226d0 / aR) * ( &
            (1.745d0 - 1.289d0/aR) * exp( -(z/H1)**2 ) + &
            (0.473d0 - 0.070d0/aR) * exp( -(z/H2)**2 ) + &
            (0.283d0 - 0.142d0/aR) * exp( -z/H3 ) )

    ! Warm ionized 
    n_wim = 0.0237d0 * exp( -(R**2 - R_Sun**2)/37.0d0**2 ) &
                    * exp( -z / 1.0d0 ) &
          + 0.0013d0 * exp( -((R-4.0d0)**2 - (R_Sun-4.0d0)**2) &
                            / (2.0d0)**2 ) &
                    * exp( -z / 0.150d0 )

    ! Hot medium (eq. 13) + scale height from eq. (52)
    Hh = 1.5d0 * (R/R_Sun)**1.65d0
    n_hm = 4.8d-4 * ( &
            0.12d0 * exp( -(R - R_Sun)/4.9d0 ) + &
            0.88d0 * exp( -((R-4.5d0)**2 - (R_Sun-4.5d0)**2) &
                          / 2.9d0**2 ) ) &
          * (R/R_Sun)**(-1.65d0) &
          * exp( -z/Hh )

    nH = max(n_mol + n_cnm + n_wnm + n_wim + n_hm, 0.0d0)

  end subroutine ferriere


  ! Miller & Bregman Flattened Beta-model

  subroutine miller(z, R, ne)
    double precision, intent(in)  :: z, R
    double precision, intent(out) :: ne

    ne = mb_n0 * (1.0d0 + (R/mb_Rc)**2 + (z/mb_zc)**2)**(-1.5d0*mb_beta) 
  end subroutine miller


  ! Total gas density

  subroutine gas_density_cgs(z, R, rho_cgs)
    double precision, intent(in)  :: z, R
    double precision, intent(out) :: rho_cgs
    double precision :: nF, nM

    call ferriere(z, R, nF)
    call miller(z, R, nM)

    rho_cgs = mp_cgs * (1.36d0 * nF + 0.62d0 * nM)
  end subroutine gas_density_cgs


  ! Miyamoto-Nagai Model - 1

  subroutine add_phi_mn1(R_cm, z_cm, M, a, b, phi_acc)
    double precision, intent(in)    :: R_cm, z_cm, M, a, b
    double precision, intent(inout) :: phi_acc
    double precision :: M_g, a_cm, b_cm, zeta, u, S

    M_g  = fmass_disk * M * 1.0d10 * const_msun
    a_cm = a * 1.0d3 * const_pc
    b_cm = b * 1.0d3 * const_pc

    zeta = sqrt(z_cm**2 + b_cm**2)
    u    = a_cm + zeta
    S    = R_cm**2 + u**2

    phi_acc = phi_acc - const_g * M_g / sqrt(S)
  end subroutine add_phi_mn1


  ! Miyamoto-Nagai Model - 2

  subroutine add_phi_mn2(R_cm, z_cm, M, a, b, phi_acc)
    double precision, intent(in)    :: R_cm, z_cm, M, a, b
    double precision, intent(inout) :: phi_acc
    double precision :: M_g, a_cm, b_cm, zeta, u, S

    M_g  = fmass_disk * M * 1.0d10 * const_msun
    a_cm = a * 1.0d3 * const_pc
    b_cm = b * 1.0d3 * const_pc

    zeta = sqrt(z_cm**2 + b_cm**2)
    u    = a_cm + zeta
    S    = R_cm**2 + u**2

    phi_acc = phi_acc - const_g * M_g / sqrt(S) &
            * (1.0d0 + a_cm * u / S)
  end subroutine add_phi_mn2


  ! Miyamoto-Nagai Model - 3

  subroutine add_phi_mn3(R_cm, z_cm, M, a, b, phi_acc)
    double precision, intent(in)    :: R_cm, z_cm, M, a, b
    double precision, intent(inout) :: phi_acc
    double precision :: M_g, a_cm, b_cm, zeta, u, S

    M_g  = fmass_disk * M * 1.0d10 * const_msun
    a_cm = a * 1.0d3 * const_pc
    b_cm = b * 1.0d3 * const_pc

    zeta = sqrt(z_cm**2 + b_cm**2)
    u    = a_cm + zeta
    S    = R_cm**2 + u**2

    phi_acc = phi_acc - const_g * M_g / sqrt(S) &
            * (1.0d0 + a_cm * u / S &
               - (a_cm**2 / 3.0d0) * (R_cm**2 - 2.0d0 * u**2) / S**2)
  end subroutine add_phi_mn3


  ! Total gravitational potential 
  ! Disks + Bulge + Halo.

  subroutine total_potential_cgs(z, R, phi)
    double precision, intent(in)  :: z, R
    double precision, intent(out) :: phi
    double precision :: R_cm, z_cm, r_sph
    double precision :: Mb_g, ab_cm, rh_cm, vh_cm
    double precision :: phi_bulge, phi_halo, phi_disks

    R_cm = R * 1.0d3 * const_pc
    z_cm = z * 1.0d3 * const_pc

    ! Bulge
    Mb_g  = M_b * 1.0d10 * const_msun
    ab_cm = a_b * 1.0d3 * const_pc
    r_sph = sqrt(R_cm**2 + z_cm**2)

    phi_bulge = -const_g * Mb_g / (r_sph + ab_cm)

    ! Halo
    rh_cm = r_h * 1.0d3 * const_pc
    vh_cm = v_h * 1.0d5   

    phi_halo = 0.5d0 * vh_cm**2 * log(R_cm**2 + z_cm**2 + rh_cm**2)

    ! Disks
    phi_disks = 0.0d0

    ! Thin stellar disk 
    call add_phi_mn3(R_cm, z_cm, M_thin1, a_thin1, b_thin, phi_disks)
    call add_phi_mn3(R_cm, z_cm, M_thin2, a_thin2, b_thin, phi_disks)
    call add_phi_mn3(R_cm, z_cm, M_thin3, a_thin3, b_thin, phi_disks)

    ! Thick stellar disk 
    call add_phi_mn1(R_cm, z_cm, M_thick1, a_thick1, b_thick, phi_disks)
    call add_phi_mn1(R_cm, z_cm, M_thick2, a_thick2, b_thick, phi_disks)
    call add_phi_mn1(R_cm, z_cm, M_thick3, a_thick3, b_thick, phi_disks)

    ! H I disk 
    call add_phi_mn2(R_cm, z_cm, M_HI1, a_HI1, b_HI, phi_disks)
    call add_phi_mn2(R_cm, z_cm, M_HI2, a_HI2, b_HI, phi_disks)
    call add_phi_mn2(R_cm, z_cm, M_HI3, a_HI3, b_HI, phi_disks)

    ! H2 disk 
    call add_phi_mn3(R_cm, z_cm, M_H21, a_H21, b_H2, phi_disks)
    call add_phi_mn3(R_cm, z_cm, M_H22, a_H22, b_H2, phi_disks)
    call add_phi_mn3(R_cm, z_cm, M_H23, a_H23, b_H2, phi_disks)

    phi = phi_bulge + phi_halo + phi_disks
  end subroutine total_potential_cgs


  ! Vertical gravitational acceleration
  ! Computed by central difference

  subroutine gravity_cgs(z, R, g_cgs)
    double precision, intent(in)  :: z, R
    double precision, intent(out) :: g_cgs
    double precision :: phi_plus, phi_minus, dz_cm
    double precision, parameter :: z_cut = 100.0d0

    dz_cm  = dz_potential * 1.0d3 * const_pc

    call total_potential_cgs(z + dz_potential, R, phi_plus)
    call total_potential_cgs(z - dz_potential, R, phi_minus)

    ! g_z = -dPhi/dz; Only magnitude
    g_cgs = abs( (phi_plus - phi_minus) / (2.0d0 * dz_cm) )

    ! Quadratic falloff beyond 100 kpc - from paper
    if (z > z_cut) then
       g_cgs = g_cgs * (z_cut / z)**2
    end if
  end subroutine gravity_cgs


  ! Final return function

  subroutine get_background(z_code, rho_code, p_code, g_code)
    double precision, intent(in)  :: z_code
    double precision, intent(out) :: rho_code, p_code, g_code
    double precision :: zq, frac, z
    integer          :: i

    z = abs(z_code) * unit_length / (1.0d3 * const_pc)

    zq = min( max(z, z_min_tab), z_max_tab * (1.0d0 - 1.0d-9) )               ! Make sure is in range
    i  = 1 + int( (zq - z_min_tab) / (z_max_tab - z_min_tab) * (nz_tab-1) )   ! Find index in table
    i  = min(max(i, 1), nz_tab - 1)

    frac     = (zq - z_tab(i)) / (z_tab(i+1) - z_tab(i))

    ! Interpolation
    rho_code = (1.0d0 - frac) * rho_tab(i) + frac * rho_tab(i+1)
    p_code   = (1.0d0 - frac) * p_tab(i)   + frac * p_tab(i+1)
    g_code   = (1.0d0 - frac) * g_tab(i)   + frac * g_tab(i+1)
  end subroutine get_background

end module mod_background

