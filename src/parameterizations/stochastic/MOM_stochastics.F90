!> Top-level module for the MOM6 ocean model in coupled mode.
module MOM_stochastics

! This file is part of MOM6. See LICENSE.md for the license.

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, update, and writing restart of stochastic physics. This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, time_type, post_data
use MOM_diag_mediator,       only : enable_averages, disable_averaging
use MOM_grid,                only : ocean_grid_type
use MOM_variables,           only : thermo_var_ptrs
use MOM_domains,             only : pass_var, pass_vector
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_error_handler,       only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, close_param_file, param_file_type
use mpp_domains_mod,         only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,         only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use MOM_domains,             only : root_PE, num_PEs
use MOM_coms,                only : Get_PElist
use MOM_EOS,                 only : calculate_density, EOS_domain
use stochastic_physics,      only : init_stochastic_physics_ocn, run_stochastic_physics_ocn

#include <MOM_memory.h>

implicit none ; private

public stochastics_init, update_stochastics,apply_skeb

!> This control structure holds parameters for the MOM_stochastics module
type, public:: stochastic_CS
  logical :: do_sppt         !< If true, stochastically perturb the diabatic
  logical :: do_skeb         !< If true, stochastically perturb the diabatic
  logical :: pert_epbl       !< If true, then randomly perturb the KE dissipation and genration terms
  integer :: id_sppt_wts  = -1 !< Diagnostic id for SPPT
  integer :: id_skeb_wts  = -1 !< Diagnostic id for SKEB
  integer :: id_skebu     = -1 !< Diagnostic id for SKEB
  integer :: id_skebv     = -1 !< Diagnostic id for SKEB
  integer :: id_diss      = -1 !< Diagnostic id for SKEB
  integer :: skeb_npass   = -1 !< number of passes of the 9-point smoother for the dissipation estimate
  integer :: id_psi       = -1 !< Diagnostic id for SPPT
  integer :: id_epbl1_wts = -1 !< Diagnostic id for epbl generation perturbation
  integer :: id_epbl2_wts = -1 !< Diagnostic id for epbl dissipation perturbation
  ! stochastic patterns
  real, allocatable :: sppt_wts(:,:)  !< Random pattern for ocean SPPT
                                     !! tendencies with a number between 0 and 2
  real, allocatable :: skeb_wts(:,:)  !< Random pattern for ocean SKEB
  real, allocatable :: epbl1_wts(:,:) !< Random pattern for K.E. generation
  real, allocatable :: epbl2_wts(:,:) !< Random pattern for K.E. dissipation
  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the
end type stochastic_CS

contains

!!   This subroutine initializes the stochastics physics control structure.
subroutine stochastics_init(dt, grid, GV, CS, param_file, diag, Time)
  real, intent(in)                       :: dt      !< time step [T ~> s]
  type(ocean_grid_type),   intent(in)    :: grid    !< horizontal grid information
  type(verticalGrid_type), intent(in)    :: GV      !< vertical grid structure
  type(stochastic_CS), pointer, intent(inout) :: CS !< stochastic control structure
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag    !< structure to regulate diagnostic output
  type(time_type), target                :: Time    !< model time

  ! Local variables
  integer, allocatable :: pelist(:) ! list of pes for this instance of the ocean
  integer :: mom_comm          ! list of pes for this instance of the ocean
  integer :: num_procs         ! number of processors to pass to stochastic physics
  integer :: iret              ! return code from stochastic physics
  integer :: pe_zero           !  root pe
  integer :: nx                ! number of x-points including halo
  integer :: ny                ! number of x-points including halo

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "ocean_stochastics_init"  ! This module's name.

  call callTree_enter("ocean_model_stochastic_init(), MOM_stochastics.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_stochastics_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%Time => Time
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! get number of processors and PE list for stochastic physics initialization
  call get_param(param_file, mdl, "DO_SPPT", CS%do_sppt, &
                 "If true, then stochastically perturb the thermodynamic "//&
                 "tendencies of T,S, amd h.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_SKEB", CS%do_skeb, &
                 "If true, then stochastically perturb the currents "//&
                 "using the stochastic kinetic energy backscatter schmem.",&
                 default=.false.)
  call get_param(param_file, mdl, "SKEB_NPASS", CS%skeb_npass, &
                 "number of passes of a 9-point smoother of the "//&
                 "dissipation estimate.", default=3)
  call get_param(param_file, mdl, "PERT_EPBL", CS%pert_epbl, &
                 "If true, then stochastically perturb the kinetic energy "//&
                 "production and dissipation terms.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)

  if (CS%do_sppt .OR. CS%pert_epbl .OR. CS%do_skeb) then
     num_procs = num_PEs()
     allocate(pelist(num_procs))
     call Get_PElist(pelist,commID = mom_comm)
     pe_zero = root_PE()
     nx = grid%ied - grid%isd + 1
     ny = grid%jed - grid%jsd + 1
     call init_stochastic_physics_ocn(dt,grid%geoLonT,grid%geoLatT,nx,ny,GV%ke, &
                                      CS%pert_epbl,CS%do_sppt,CS%do_skeb,pe_zero,mom_comm,iret)
     if (iret/=0)  then
         call MOM_error(FATAL, "call to init_stochastic_physics_ocn failed")
         return
     endif

     if (CS%do_sppt) allocate(CS%sppt_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
     if (CS%do_skeb) allocate(CS%skeb_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
     if (CS%pert_epbl) then
       allocate(CS%epbl1_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
       allocate(CS%epbl2_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
     endif
  endif
  CS%id_sppt_wts = register_diag_field('ocean_model', 'sppt_pattern', CS%diag%axesT1, Time, &
       'random pattern for sppt', 'None')
  CS%id_skeb_wts = register_diag_field('ocean_model', 'skeb_pattern', CS%diag%axesT1, Time, &
       'random pattern for skeb', 'None')
  CS%id_epbl1_wts = register_diag_field('ocean_model', 'epbl1_wts', CS%diag%axesT1, Time, &
      'random pattern for KE generation', 'None')
  CS%id_epbl2_wts = register_diag_field('ocean_model', 'epbl2_wts', CS%diag%axesT1, Time, &
      'random pattern for KE dissipation', 'None')
  CS%id_skebu = register_diag_field('ocean_model', 'skebu', CS%diag%axesTL, Time, &
       'zonal current perts', 'None')
  CS%id_skebv = register_diag_field('ocean_model', 'skebv', CS%diag%axesTL, Time, &
       'zonal current perts', 'None')
  CS%id_diss = register_diag_field('ocean_model', 'diss', CS%diag%axesTL, Time, &
       'dissipation', 'None')
  CS%id_psi  = register_diag_field('ocean_model', 'psi', CS%diag%axesTL, Time, &
       'stream function', 'None')

  if (CS%do_sppt .OR. CS%pert_epbl) &
    call MOM_mesg('            === COMPLETED MOM STOCHASTIC INITIALIZATION =====')

  call callTree_leave("ocean_model_init(")

end subroutine stochastics_init

!> update_ocean_model uses the forcing in Ice_ocean_boundary to advance the
!! ocean model's state from the input value of Ocean_state (which must be for
!! time time_start_update) for a time interval of Ocean_coupling_time_step,
!! returning the publicly visible ocean surface properties in Ocean_sfc and
!! storing the new ocean properties in Ocean_state.
subroutine update_stochastics(CS)
  type(stochastic_CS),      intent(inout) :: CS        !< diabatic control structure
  call callTree_enter("update_stochastics(), MOM_stochastics.F90")

! update stochastic physics patterns before running next time-step
  call run_stochastic_physics_ocn(CS%sppt_wts,CS%skeb_wts,CS%epbl1_wts,CS%epbl2_wts)

  return
end subroutine update_stochastics

subroutine apply_skeb(grid,GV,CS,uc,vc,thickness,tv,dt,Time_end)

  type(ocean_grid_type),      intent(in) :: grid      !< ocean grid structure
  type(verticalGrid_type),    intent(in) :: GV        !< ocean vertical grid
  type(stochastic_CS),        intent(in) :: CS        !< stochastic control structure

  real, dimension(SZIB_(grid),SZJ_(grid),SZK_(GV)), intent(inout) :: uc        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJB_(grid),SZK_(GV)), intent(inout) :: vc        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJ_(grid),SZK_(GV)),  intent(in)    :: thickness !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                            intent(in)    :: tv       !< points to thermodynamic fields
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
! locals

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: diss
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: psi_3d 
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: ustar
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: vstar
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: psi,diss_tmp

  real    :: shr,ten,tot,kh,con_skeb,con_diss
  integer :: i,j,k,iter
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  
  call callTree_enter("apply_skeb(), MOM_stochastics.F90")
  ALLOC_(diss(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke))
  ALLOC_(diss_tmp(grid%isd:grid%ied,grid%jsd:grid%jed))
  ALLOC_(psi_3d(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke))
  ALLOC_(psi(grid%isd:grid%ied,grid%jsd:grid%jed))
  ALLOC_(ustar(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke))
  ALLOC_(vstar(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke))

  con_diss=1.0
  con_skeb=1.0!0.001
!con_conv=1.0
! fill in halos with zeros
  DO k=1,GV%ke
     DO j=grid%jsd,grid%jed
        DO i=grid%isd,grid%ied
           diss(i,j,k)=0.0
        ENDDO
     ENDDO
  ENDDO
  DO j=grid%jsd,grid%jed
     DO i=grid%isd,grid%ied
        psi(i,j)=0.0
     ENDDO
  ENDDO
  
 !kh needs to be scaled

  kh=1!(120*111)**2
  DO k=1,GV%ke
     DO j=grid%jsc,grid%jec
        DO i=grid%isc,grid%iec
           ! Shear
           shr = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)*grid%IdxCv(i,J) + &
                 (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)* grid%IdyCu(I,j)
           ! Tension
           ten = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)* grid%IdyCv(i,J) + &
                 (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)* grid%IdxCu(I,j)
     
           tot = sqrt( shr**2 + ten**2 ) * grid%mask2dT(i,j)
           diss(i,j,k) = tot**3*kh*grid%areaT(i,j)!!**2
        ENDDO
     ENDDO
  ENDDO
  
  call pass_var(diss, grid%Domain)
  ! smooth dissipation x times
  do iter=1,CS%skeb_npass
     DO k=1,GV%ke
        !print*,'smooth diss',iter,k,minval(diss),maxval(diss)
        DO j=grid%jsc,grid%jec
           DO i=grid%isc,grid%iec
              diss_tmp(i,j)=(diss(i-1,j-1,k)+diss(i  ,j-1,k)+diss(i+1,j-1,k)+&
                             diss(i-1,j  ,k)+diss(i  ,j  ,k)+diss(i+1,j  ,k)+&
                             diss(i-1,j+1,k)+diss(i  ,j+1,k)+diss(i+1,j+1,k))/9.0* &
                             grid%mask2dT(i,j)
           ENDDO
        ENDDO
        DO j=grid%jsc,grid%jec
           DO i=grid%isc,grid%iec
              diss(i,j,k)=diss_tmp(i,j)
           ENDDO
        ENDDO
     ENDDO
     call pass_var(diss, grid%Domain)
  ENDDO

  !print*,'before psi', minval(diss),maxval(diss),minval(CS%skeb_wts(:,:)),maxval(CS%skeb_wts(:,:))
  DO k=1,GV%ke
     DO j=grid%jsc,grid%jec
        DO i=grid%isc,grid%iec
           psi(i,j) = con_skeb * sqrt( con_diss * diss(i,j,k)) *CS%skeb_wts(i,j)*dt
           psi_3d(i,j,k)=psi(i,j)
        ENDDO
     ENDDO
     call pass_var(psi, grid%Domain)
     DO j=grid%jsc,grid%jec
        DO i=grid%isc,grid%iec
            ustar(i,j,k) =   ( psi(i,j)-psi(i,j-1))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)* grid%IdyCu(i,j)
            vstar(i,j,k) = - ( psi(i,j)-psi(i-1,j))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)* grid%IdxCv(i,j)
            uc(I,j,k)=uc(I,j,k)+ustar(I,j,k)
            vc(i,J,k)=vc(i,J,k)+vstar(i,J,k)
       ENDDO
     ENDDO
  ENDDO
  call enable_averages(dt, Time_end, CS%diag)
  if (CS%id_diss > 0) then
     call post_data(CS%id_diss, diss(:,:,:), CS%diag)
  endif
  if (CS%id_skeb_wts > 0) then
     call post_data(CS%id_skeb_wts, CS%skeb_wts, CS%diag)
  endif
  if (CS%id_skebu > 0) then
     call post_data(CS%id_skebu, ustar(:,:,:), CS%diag)
  endif
  if (CS%id_skebv > 0) then
     call post_data(CS%id_skebv, vstar(:,:,:), CS%diag)
  endif
  if (CS%id_psi > 0) then
     call post_data(CS%id_psi, psi_3d(:,:,:), CS%diag)
  endif
  call disable_averaging(CS%diag)
  deallocate(diss)
  deallocate(diss_tmp)
  deallocate(psi_3d)
  deallocate(ustar)
  deallocate(vstar)
  deallocate(psi)
  return
end subroutine apply_skeb

end module MOM_stochastics

