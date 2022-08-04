module scatter

  implicit none
  private
  public :: compton_polarized, lorentz_trans, sample_scattered_photons, calc_chi, rotate_q_chi,&
            calc_thomson_diff, calc_kn_diff, update_stokes, lorentz_inverse

contains

  subroutine compton_polarized(k_f_pre, beta,  gamma, n, stokes_pre, k_f_post, stokes_post)
    use mtmod

    implicit none

    real(8), dimension(4), intent(in)  :: k_f_pre
    real(8), dimension(3), intent(in)  :: stokes_pre
    real(8), dimension(4), intent(out) :: k_f_post
    real(8), dimension(3), intent(out) :: stokes_post
    real(8), intent(in) :: beta, gamma
    real(8), dimension(4) :: k_e_pre, k_e_post
    real(8), dimension(3), intent(in) :: n

    real(8) :: q_eex_pre, u_eex_pre, mu_k, k0_e_pre
    real(8) :: chi_pre,chi_post,cos_chi_post
    real(8) :: q_ein_pre
    real(8) :: theta_pre, theta_post, phi_pre, phi_post
    real(8) :: cross_section_norm
    real(8) :: r
    real(8) :: n_k, n_k_post
    real(8) :: k0_e_post_dummy  

    q_eex_pre = stokes_pre(1)
    u_eex_pre = stokes_pre(2)
    
    call lorentz_trans(k_f_pre, k_e_pre, beta, gamma, n, k0_e_pre)

    ! ...
    ! Determine Scattered Photons by Rejection Method
    ! <-

    do
       
      call sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k,theta_pre, theta_post, phi_pre, phi_post)

      call calc_chi(k_e_pre, k_e_post, mu_k, chi_pre, &
           chi_post,theta_pre, theta_post, phi_pre, phi_post, cos_chi_post)

      call rotate_q_chi(q_eex_pre, u_eex_pre, q_ein_pre, chi_pre, chi_post)

      if (k0_e_pre .lt. 1.0d-4) then  ! Thomson Scattering Limit
        call calc_thomson_diff(k0_e_pre, q_ein_pre, mu_k, &
                                k0_e_post_dummy, cross_section_norm)
      else  ! Compton Scattering
        call calc_kn_diff(k0_e_pre, q_ein_pre, mu_k, &
                          k0_e_post_dummy, cross_section_norm)
     end if
     r = grnd()
           if (cross_section_norm .gt. r) exit  ! Rejection method
        end do
        
        call update_stokes(stokes_pre, k_e_pre, k_e_post, mu_k, chi_pre, chi_post, stokes_post)

        call lorentz_inverse(k_f_post, k_e_post, beta, gamma, n)

    return
  end subroutine compton_polarized

!===============================================================================
! Convert  photons' 4-momentum from laboratry frame to electron rest frame
!-------------------------------------------------------------------------------
! input
!  - k_f_pre : 4-momentum of photons in laboratry frame
!  - n : direnction of electron in laboratry frame
!  - beta : electron's velocity normalized by c   
! output
!   - k_e_pre : 4-momentum of photons in electron rest frame
!   - k0_e_pre : energy of photons in electron rest frame
!===============================================================================
  subroutine lorentz_trans(k_f_pre, k_e_pre, beta, gamma, n, k0_e_pre)

    implicit none

    real(8), dimension(4), intent(in) :: k_f_pre
    real(8), intent(in) :: beta,gamma
    real(8), dimension(3), intent(in) :: n
    real(8) :: n_k
    real(8), dimension(4), intent(out) :: k_e_pre
    real(8), intent(out) :: k0_e_pre

    
    n_k = n(1)*k_f_pre(2) + n(2)*k_f_pre(3) + n(3)*k_f_pre(4) 

    k_e_pre(2) = k_f_pre(2) + (gamma - 1.0d0)*n_k*n(1) - gamma*beta*k_f_pre(1)*n(1)
    k_e_pre(3) = k_f_pre(3) + (gamma - 1.0d0)*n_k*n(2) - gamma*beta*k_f_pre(1)*n(2)
    k_e_pre(4) = k_f_pre(4) + (gamma - 1.0d0)*n_k*n(3) - gamma*beta*k_f_pre(1)*n(3)
    k_e_pre(1) = sqrt(k_e_pre(2)**2 + k_e_pre(3)**2 + k_e_pre(4)**2)
    k0_e_pre  = k_e_pre(1)

    return
  end subroutine lorentz_trans
  
  !===============================================================================
  ! Determine scatterd photons' direction randomly and struct its 4-momentum
  !-------------------------------------------------------------------------------
  ! input
  !   - k0_e_pre : energy of pre-photons in Electron Rest Frame
  ! output
  !   - k_e_post : 4-momentum of scattered photons in Electron Rest Frame
  !   - mu_k     : scattering angle between pre- and post-photons
  !===============================================================================
  
  subroutine sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k, theta_pre, theta_post, phi_pre, phi_post)

    use mtmod

    implicit none

    real(8), intent(in)  ::  k0_e_pre
    real(8), dimension(4), intent(in) :: k_e_pre
    real(8), dimension(4), intent(out) :: k_e_post
    real(8), intent(out) :: mu_k

    real(8) :: g, phi, k0_e_post,r_
    real(8), intent(out) :: theta_pre, phi_pre, theta_post, phi_post

    real(8) :: pi
    pi = 4.0d0*atan(1.0d0)

    ! ...
    ! Determine the direction of scattered photons randomly
    ! and its energy with Copmton's formula
    ! to struct its 4-momentum in the electron rest frame
    ! <-

    theta_pre = acos(k_e_pre(4)/k_e_pre(1)) ! Initial theta

    if(k_e_pre(3) .gt. 0)then
       phi_pre = acos(k_e_pre(2)/(sqrt(k_e_pre(2)*k_e_pre(2) + k_e_pre(3)*k_e_pre(3))))
    else
       phi_pre = 2.0d0*pi - acos(k_e_pre(2)/(sqrt(k_e_pre(2)*k_e_pre(2) + k_e_pre(3)*k_e_pre(3))))
    end if ! Initial phi
    
    !call random_number(r_) ! use mt
    r_ = grnd()
    theta_post = pi*r_ !angle of scattering !3/11 with determination of mu check
    
    mu_k = cos(theta_post)

    ! call random_number(g) !use mt
    g = grnd()
     phi_post = 2.0d0*pi*g !determination of phi

     
  ! k0_e_post = k0_e_pre/( 1.0d0 + k0_e_pre*( 1.0d0 - mu_k )) !Compton 
    k0_e_post = k0_e_pre !Thomson
    k_e_post(1) = k0_e_post
    k_e_post(4) = k0_e_post*(sin(theta_post)*cos(phi_post)*sin(theta_pre) + cos(theta_pre)*cos(theta_post))
    k_e_post(2) = k0_e_post*(-sin(theta_post)*cos(phi_post)*cos(theta_pre)*cos(phi_pre)&
         +cos(theta_post)*sin(theta_pre)*cos(phi_pre) + sin(theta_post)*sin(phi_post)*sin(phi_pre))
    k_e_post(3) = k0_e_post*(-sin(theta_post)*cos(phi_post)*cos(theta_pre)*sin(phi_pre)&
         +cos(theta_post)*sin(theta_pre)*sin(phi_pre) - sin(theta_post)*sin(phi_post)*cos(phi_pre)) !3/11 check
    ! ->

    return
  end subroutine sample_scattered_photons

!===============================================================================
! Calculate Rotation Angle CHI (See Nagirner & Poutanen 1993, Davis et al. 2009)
!-------------------------------------------------------------------------------
! input
!   - k_e_pre  : 4-momentum of pre-photons in Electron Rest Frame
!   - k_e_post : 4-momentum of scattered photons in Electron Rest Frame
!   - mu_k     : scattering angle between pre- and post-photons
! output
!   - chi_pre   : Nagirner1993 (E2), Davis2009 (A5)
!   - cos_chi_post, chi_post :              (E3),           (A6)
!===============================================================================
!-------------------------------------------------------------------------------
! TODO(from Iwanaga)
!   - Is mu_k really need?
!     You can calculate cosine bitween pre- and post-photons.
!   - k@RAIKOU isn't unit direction vectors but 4-momentum.
!     There can be differnce from references such as Davis2009.
!      -> Alreay settled.
!         Made direction vectors from 4-momentum.
!-------------------------------------------------------------------------------
  subroutine calc_chi(k_e_pre, k_e_post, mu_k, &
                       chi_pre, chi_post, theta_pre, theta_post, phi_pre, phi_post, cos_chi_post)

    implicit none

    real(8), dimension(4), intent(in)  :: k_e_pre, k_e_post
    real(8), intent(in) :: mu_k, theta_pre, theta_post, phi_pre, phi_post
    real(8), intent(out) :: chi_pre,chi_post
    real(8), intent(out) :: cos_chi_post
    real(8) :: sin_chi_pre,cos_chi_pre, sin_chi_post

    integer :: i
    real(8), dimension(4) :: unit_k_e_pre, unit_k_e_post
    real(8) :: pi
    real(8) :: denominator,denominator_2
    pi = 4.0d0*atan(1.0d0)
    
    do i = 1,4
       unit_k_e_pre(i) = k_e_pre(i)/k_e_pre(1)
       unit_k_e_post(i) = k_e_post(i)/k_e_post(1) !3/11 cleaning
    end do

    denominator = sin(theta_post)*sqrt(1-unit_k_e_pre(4)*unit_k_e_pre(4))
    denominator_2 = sin(theta_post)*sqrt(1-unit_k_e_post(4)*unit_k_e_post(4))

    cos_chi_pre = (unit_k_e_post(4)-unit_k_e_pre(4)*mu_k)/denominator
    sin_chi_pre = (unit_k_e_pre(3)*unit_k_e_post(2)-unit_k_e_pre(2)*unit_k_e_post(3))/denominator ! cleanig and to 219

    if(cos_chi_pre .gt. 1.0d0)then
       cos_chi_pre = 1.0d0 !3/11 cleaning
    end if
    if(cos_chi_pre .lt. -1.0d0)then
       cos_chi_pre = -1.0d0
    end if

    if(sin_chi_pre .gt. 0)then
       chi_pre = acos(cos_chi_pre)
    else
       chi_pre = 2.0d0*pi - acos(cos_chi_pre)
    end if
    

    cos_chi_post = (unit_k_e_pre(4)-unit_k_e_post(4)*mu_k)/denominator_2
    sin_chi_post = (unit_k_e_post(3)*unit_k_e_pre(2)-unit_k_e_post(2)*unit_k_e_pre(3))/denominator_2

    if(cos_chi_post .gt. 1.0d0)then
       cos_chi_post = 1.0d0
    end if
    if(cos_chi_post .lt. -1.0d0)then
       cos_chi_post = -1.0d0
    end if  

    if(sin_chi_post .gt. 0)then
       chi_post = acos(cos_chi_post)
    else
       chi_post = 2.0d0*pi - acos(cos_chi_post)
    end if
    
    
    ! ->

    return
  end subroutine calc_chi

!===============================================================================
! Update photon's linear polarization q under rotating basis by L(chi)
! See Nagirner & Poutanen 1993 (C12) and Davis 2009 (A4)
!-------------------------------------------------------------------------------
! input
!   - q1      : linear polarization pararell/perpendicular
!               to the scattering plane before rotating chi
!   - u1      : 45 degrees linear polarization before rotating chi
!   - chi_pre : pre-rotation angle
!   - chi_post : post-rotation angle
! output
!   - q2      : linear polarization after rotating chi
!               q2 = q1*cos(2*chi) + u1*sin(2*chi)
!===============================================================================
  subroutine rotate_q_chi(q1, u1, q2, chi_pre, chi_post)

    implicit none

    real(8), intent(in)  :: q1, u1, chi_pre, chi_post
    real(8), intent(out) :: q2

    real(8) :: cos2chi, sin2chi

    cos2chi = cos(chi_pre)*cos(chi_pre) - sin(chi_pre)*sin(chi_pre)
    sin2chi = 2.0d0*sin(chi_pre)*cos(chi_pre)

    q2 = q1*cos2chi + u1*sin2chi

    return
  end subroutine rotate_q_chi

!===============================================================================
! Calculate Thomson differential cross section (normalized : max = 1)
!-------------------------------------------------------------------------------
! input
!   - k0_e_pre : energy of pre-photons in Electron Rest Frame
!   - q_ein_pre : linear polarization of pre-photons in e-internal frame(k)
!   - mu_k     : scattering angle between pre- and post-photons
! output
!   - k0_e_post : energy of scatterd photons in Electron Rest Frame
!   - cross_section_norm : Thomson differential cross section
!                          (normalized : max = 1)
!==============================================================================
  subroutine calc_thomson_diff(k0_e_pre, q_ein_pre, mu_k, &
                                k0_e_post, cross_section_norm)

    implicit none

    real(8), intent(in) :: k0_e_pre, q_ein_pre, mu_k
    real(8), intent(out) :: k0_e_post, cross_section_norm

    cross_section_norm = 0.5d0*( 2.0d0 - &
                          ( 1.0d0 + q_ein_pre )*( 1.0d0 - mu_k*mu_k ) )

    k0_e_post = k0_e_pre

    return
  end subroutine calc_thomson_diff

!===============================================================================
! Calculate Klein-Nishina differential cross section (normalized : max = 1)
!-------------------------------------------------------------------------------
! input
!   - k0_e_pre : energy of pre-photons in Electron Rest Frame
!   - q_ein_pre : linear polarization of pre-photons in e-internal frame(k)
!   - mu_k     : scattering angle between pre- and post-photons
! output
!   - k0_e_post : energy of scatterd photons in Electron Rest Frame
!   - cross_section_norm : Thomson differential cross section
!                          (normalized : max = 1)
!===============================================================================
  subroutine calc_kn_diff(k0_e_pre, q_ein_pre, mu_k, &
                          k0_e_post, cross_section_norm)

    implicit none

    real(8), intent(in) :: k0_e_pre, q_ein_pre, mu_k
    real(8), intent(out) :: k0_e_post, cross_section_norm

    ! real(8) :: k0_e_post_min, k0_e_post_max  ! safety
    real(8) :: k0_e_pre_inv, k0_e_post_inv
    ! k0_e_post_min = k0_e_pre/( 1.0d0 + 2.0d0*k0_e_pre )  ! safety
    k0_e_post = k0_e_pre                             !Thomson

    k0_e_pre_inv = 1.0d0/k0_e_pre
  ! k0_e_post = k0_e_pre/( 1.0d0 + k0_e_pre*( 1.0d0 - mu_k ) )
    k0_e_post_inv = 1.0d0/k0_e_post
    ! k0_e_post = max(k0_e_post, k0_e_post_min)  ! safety
    ! k0_e_post = min(k0_e_post, k0_e_post_max)  ! safety

    cross_section_norm = 0.5d0 * k0_e_post*k0_e_post*k0_e_pre_inv*k0_e_pre_inv &
                          *( k0_e_post*k0_e_pre_inv + k0_e_pre*k0_e_post_inv &
                          - ( 1.0d0 + q_ein_pre )*( 1.0d0 - mu_k*mu_k ) )

    return
  end subroutine calc_kn_diff

  !===============================================================================
  ! Update Stokes paraeters
  !-------------------------------------------------------------------------------
  ! input
  !   - k_e_pre : energy of pre-photons in Electron Rest Frame
  !   - k_e_post : energy of post-photons in Electron Rest Frame
  !   - stokes_pre : stokes parametes of pre-photons in e-external frame
  !   - mu_k     : scattering angle between pre- and post-photons
  !   - chi_pre : pre-rotation angle
  !   - chi_post : post-rotation angle
  ! output
  !   - stokes_post : stokes parametes of post-photons in e-external frame 
  !===============================================================================
  subroutine update_stokes(stokes_pre, k_e_pre, k_e_post, mu_k, chi_pre, chi_post, stokes_post)

    implicit none

    real(8), dimension(3), intent(in) :: stokes_pre
    real(8), dimension(4), intent(in) :: k_e_pre, k_e_post
    real(8),intent(in) :: mu_k, chi_pre, chi_post
    real(8), dimension(3), intent(out) :: stokes_post

    real(8) :: cos2chi_pre, sin2chi_pre, cos2chi_post, sin2chi_post
    real(8) :: s_a, s_b, s_c, s_d
    real(8) ::q_eex_pre, u_eex_pre, q_eex_post, u_eex_post, intensity
    real(8) :: k0_e_pre, k0_e_post

    q_eex_pre = stokes_pre(1)
    u_eex_pre = stokes_pre(2)

    k0_e_pre = k_e_pre(1)
    k0_e_post = k_e_post(1)

    cos2chi_pre = cos(chi_pre)*cos(chi_pre) - sin(chi_pre)*sin(chi_pre)
    sin2chi_pre = 2.0d0*sin(chi_pre)*cos(chi_pre)
    cos2chi_post = cos(chi_post)*cos(chi_post) - sin(chi_post)*sin(chi_post)
    sin2chi_post = 2.0d0*sin(chi_post)*cos(chi_post)

    s_a = mu_k*mu_k + 1.0d0
    s_b = k0_e_post/k0_e_pre + k0_e_pre/k0_e_post - 2.0d0
    s_c = (s_a - 2.0d0) ! probably -(s_a-2.0d0)...
    s_d = 2.0d0*mu_k


    intensity = s_a + s_b + s_c*(q_eex_pre*cos2chi_pre + u_eex_pre*sin2chi_pre)


    q_eex_post = s_c*cos2chi_post &
         + q_eex_pre*(s_a*cos2chi_pre*cos2chi_post + s_d*sin2chi_pre*sin2chi_post) &
         + u_eex_pre*(s_a*sin2chi_pre*cos2chi_post - s_d*cos2chi_pre*sin2chi_post)

    u_eex_post = s_c*sin2chi_post &
         + q_eex_pre*(s_a*cos2chi_pre*sin2chi_post - s_d*sin2chi_pre*cos2chi_post) &
         + u_eex_pre*(s_a*sin2chi_pre*sin2chi_post + s_d*cos2chi_pre*cos2chi_post)

    stokes_post(1) = q_eex_post/intensity
    stokes_post(2) = u_eex_post/intensity
    stokes_post(3) = stokes_pre(3)


    return
  end subroutine update_stokes 

  !===============================================================================
  ! Convert  photons' 4-momentum from electron rest frame to laboratry frame
  !-------------------------------------------------------------------------------
  ! input
  !  - k_e_post : 4-momentum of photons in electron rest frame
  !  - n : direnction of electron in laboratry frame
  !  - beta : electron's velocity normalized by c
  ! output
  !   - k_f_post : 4-momentum of photons in laboratry frame
  !===============================================================================
  subroutine lorentz_inverse(k_f_post, k_e_post, beta, gamma, n)

    implicit none

    real(8), dimension(4), intent(in) :: k_e_post
    real(8), dimension(3), intent(in) :: n
    real(8), intent(in) :: beta, gamma
    real(8), dimension(4), intent(out) :: k_f_post
    real(8) :: n_k_post

    n_k_post = n(1)*k_e_post(2) + n(2)*k_e_post(3) + n(3)*k_e_post(4)

    k_f_post(2) = k_e_post(2) + (gamma - 1.0d0)*n_k_post*n(1) + gamma*beta*k_e_post(1)*n(1)
    k_f_post(3) = k_e_post(3) + (gamma - 1.0d0)*n_k_post*n(2) + gamma*beta*k_e_post(1)*n(2)
    k_f_post(4) = k_e_post(4) + (gamma - 1.0d0)*n_k_post*n(3) + gamma*beta*k_e_post(1)*n(3)
    k_f_post(1) = sqrt(k_f_post(2)**2 + k_f_post(3)**2 + k_f_post(4)**2) 


    return
    end subroutine lorentz_inverse
    
end module scatter


  
 
