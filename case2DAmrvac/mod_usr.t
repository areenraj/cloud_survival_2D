module mod_usr
  use mod_global_parameters
  use mod_usr_methods
  use mod_hd
  use mod_background

  implicit none

  ! Cloud initial position, size and velocity in pc and km/s
  double precision, parameter :: cloud_x   = 0.0d0    
  double precision, parameter :: cloud_z   = 5.0d3    
  double precision, parameter :: cloud_r   = 1.0d2    
  double precision, parameter :: cloud_vz  = -1.0d2   
  double precision, parameter :: cloud_nH  = 0.162d0  ! cm^-3

  double precision :: rho_cl 
  double precision :: rho_cl_bulk                 ! refinement threshold density

contains

  subroutine usr_init()

    unit_length        = 3.0857d18   ! 1 pc
    unit_velocity      = 1.0d5       ! 1 km/s
    unit_numberdensity = 0.1d0       ! 0.1 cm^-3

    call usr_params_read(par_files)
    call bg_params_read(par_files)

    usr_set_parameters => usr_setup_background
    usr_init_one_grid  => usr_initial_conditions
    usr_gravity        => usr_gravity_field
    usr_refine_grid    => usr_amr_refine

    call hd_activate()

  end subroutine usr_init

  subroutine usr_setup_background
    call init_background()

    rho_cl      = mp_cgs * cloud_nH / unit_density
    rho_cl_bulk = 0.1d0 * rho_cl
  end subroutine usr_setup_background

  subroutine usr_initial_conditions(ixI^L, ixO^L, w, x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    double precision :: rho_bg, p_bg, g_bg     ! background 
    double precision :: rho_loc          ! cloud 
    double precision :: dist, xx, zz
    integer          :: ix^D

    ! Get and allocate background
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       xx = x(ix^D, 1)     ! horizontal coordinate (code = pc)
       zz = x(ix^D, 2)     ! vertical coordinate   (code = pc)

       call get_background(zz, rho_bg, p_bg, g_bg)

       dist = dsqrt((xx - cloud_x)**2 + (zz - cloud_z)**2)

       ! Sharp boundary
       if (dist < cloud_r) then
          rho_loc = rho_cl
        else
          rho_loc = rho_bg
       end if

       w(ix^D, rho_)    = rho_loc
       w(ix^D, mom(1))  = 0.0d0    
       w(ix^D, p_)      = p_bg                         ! pressure equilibrium

       if (dist < cloud_r) then   
         w(ix^D, tracer(1)) = 1.0d0              
         w(ix^D, mom(2))    = cloud_vz
       else 
         w(ix^D, tracer(1)) = 0.d0                     ! C = 0              
         w(ix^D, mom(2))    = 0.0d0 
       end if

    {enddo\}

    call hd_to_conserved(ixI^L, ixO^L, w, x)

  end subroutine usr_initial_conditions


  ! Gravity source term

  subroutine usr_gravity_field(ixI^L, ixO^L, wCT, x, gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: wCT(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: gravity_field(ixI^S, ndim)

    double precision :: rho_bg, p_bg, g_bg, zz
    integer          :: ix^D

    gravity_field(ixO^S, 1) = 0.0d0

    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      zz = x(ix^D, 2)

      call get_background(zz, rho_bg, p_bg, g_bg)

      ! setting sign
      if (zz >= 0.0d0) then
        gravity_field(ix^D, 2) = -g_bg
      else
        gravity_field(ix^D, 2) =  g_bg
      end if
    {enddo\}

  end subroutine usr_gravity_field


  subroutine usr_amr_refine(igrid, level, ixI^L, ixO^L, qt, w, x, &
                            refine, coarsen)
    integer, intent(in)          :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S, nw), x(ixI^S, 1:ndim)
    integer, intent(inout)       :: refine, coarsen

    integer, parameter          :: level_bg = 2
    double precision, parameter :: C_min    = 1.0d-3
    logical                     :: has_tracer, has_bulk

    if (hd_n_tracer < 1) return

    has_tracer = any(w(ixO^S, tracer(1)) > C_min * w(ixO^S, rho_))
    has_bulk   = any(w(ixO^S, tracer(1)) > C_min * w(ixO^S, rho_) .and. &
                     w(ixO^S, rho_)      > rho_cl_bulk)

    if (has_bulk) then
       refine  = 1                           ! bulk cloud: force refinement
       coarsen = 0
    else if (has_tracer) then
       refine = -1                           ! diluted remnant: stop refining it
       if (level > level_bg) coarsen = 1     ! and release the cells
    end if
    ! no tracer at all -> refine = coarsen = 0 -> Lohner decides, untouched
  end subroutine usr_amr_refine
end module mod_usr

