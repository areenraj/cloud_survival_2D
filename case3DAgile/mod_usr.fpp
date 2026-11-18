#:mute
#:include "../src/mod_gpu_directives.fpp"
#:endmute

module mod_usr
  use mod_amrvac
  use mod_physics
  use mod_background

  implicit none

  ! Cloud initial position, size and velocity in pc and km/s
  double precision, parameter :: cloud_x   = 0.0d0
  double precision, parameter :: cloud_y   = 0.0d0
  double precision, parameter :: cloud_z   = 5.0d3
  double precision, parameter :: cloud_r   = 1.0d2
  double precision, parameter :: cloud_vz  = -1.0d2
  double precision, parameter :: cloud_nH  = 0.162d0  ! cm-3,cm-3

  double precision :: rho_cl
  double precision :: rho_cl_bulk !refinement threshold density
  ${GPU_DECLARE_CREATE('rho_cl_bulk')}$

contains

  subroutine usr_init()

    unit_length        = 3.0857d18   ! 1 pc
    unit_velocity      = 1.0d5       ! 1 km/s
    unit_numberdensity = 0.1d0       ! 0.1 cm-3,cm-3

    call set_coordinate_system("Cartesian_3D")

    usr_set_parameters => usr_setup_background
    usr_init_one_grid  => usr_initial_conditions

    call phys_activate()

  end subroutine usr_init
  
  subroutine usr_setup_background()
    call init_background()

    rho_cl      = mp_cgs * cloud_nH / unit_density
    rho_cl_bulk = 0.01d0 * rho_cl
    ${GPU_UPDATE_DEVICE('rho_cl_bulk')}$
  end subroutine usr_setup_background

  subroutine usr_initial_conditions(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndim)

    double precision :: rho_bg, p_bg, g_bg     ! background
    double precision :: rho_loc                ! cloud
    double precision :: dist, xx, yy, zz
    integer          :: ix1,ix2,ix3

    ! Get and allocate background
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       xx = x(ix1,ix2,ix3, 1)     ! horizontal coordinate (code = pc)
       yy = x(ix1,ix2,ix3, 2)     ! horizontal coordinate (code = pc)
       zz = x(ix1,ix2,ix3, 3)     ! vertical coordinate   (code = pc)

       call get_background(zz, rho_bg, p_bg, g_bg)

       dist = dsqrt((xx - cloud_x)**2 + (yy - cloud_y)**2 + (zz - cloud_z)**2)

          ! Sharp boundary
       if (dist < cloud_r) then
          rho_loc = rho_cl
       else
          rho_loc = rho_bg
       end if

       w(ix1,ix2,ix3, rho_)    = rho_loc
       w(ix1,ix2,ix3, mom(1))  = 0.0d0
       w(ix1,ix2,ix3, mom(2))  = 0.0d0
       w(ix1,ix2,ix3, p_)      = p_bg !pressure equilibrium

       if (dist < cloud_r) then
         w(ix1,ix2,ix3, tracer(1)) = rho_loc 
         w(ix1,ix2,ix3, mom(3))    = cloud_vz
       else
         w(ix1,ix2,ix3, tracer(1)) = 0.d0 !C = 0
         w(ix1,ix2,ix3, mom(3))    = 0.0d0
       end if

    enddo
    enddo
    enddo

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine usr_initial_conditions


  ! Gravity source term

  pure real(dp) function gravity_field(wCT, x, idim) result(field)
    ${GPU_ROUTINE_SEQ()}$
    real(dp), intent(in)        :: wCT(nw_phys)
    real(dp), intent(in)        :: x(1:ndim)
    integer, value, intent(in)  :: idim

    double precision :: rho_bg, p_bg, g_bg, zz

    field = 0.0d0
    if (idim /= 3) return

    zz = x(3)

    call get_background(zz, rho_bg, p_bg, g_bg)

    ! setting sign
    if (zz >= 0.0d0) then
       field = -g_bg
    else
       field =  g_bg
    end if

  end function gravity_field


  subroutine usr_refine_grid(level, qt, w, x, refineflag, coarsenflag, norefineflag, nocoarsenflag)
    ${GPU_ROUTINE_SEQ()}$
    integer, intent(in)          :: level
    double precision, intent(in) :: qt, w(1:nw), x(1:ndim)
    logical, intent(inout)       :: refineflag, coarsenflag, norefineflag, nocoarsenflag

    integer, parameter          :: level_bg = 3
    double precision, parameter :: C_min    = 1.0d-6

    if (qt <= 0.0d0 .and. &
        (x(1)-cloud_x)**2 + (x(2)-cloud_y)**2 + (x(3)-cloud_z)**2 < (1.2d0*cloud_r)**2) then
       refineflag  = .true.  ; norefineflag  = .false.
       coarsenflag = .false. ; nocoarsenflag = .true.
       return
    end if

    ! w is CONSERVED here, so tracer(1) holds rho*C
    if (w(tracer(1)) > C_min * w(rho_)) then
       if (w(rho_) > rho_cl_bulk) then
          refineflag  = .true.  ; norefineflag  = .false.
          coarsenflag = .false. ; nocoarsenflag = .true.
       else if (level > level_bg) then
          refineflag  = .false. ; norefineflag  = .true.
          coarsenflag = .true.  ; nocoarsenflag = .false.
       else
          refineflag  = .false. ; norefineflag  = .true.
          coarsenflag = .false. ; nocoarsenflag = .true.
       end if
    else
       refineflag  = .false.            ; norefineflag  = .true.
       coarsenflag = (level > level_bg) ; nocoarsenflag = .false.
    end if

  end subroutine usr_refine_grid

end module mod_usr


