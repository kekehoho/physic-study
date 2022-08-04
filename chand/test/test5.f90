program main

  !  use KN_diff
  use scatter
  use mtmod

  implicit none
  
  integer :: iter, iter_max, fileNum,j,i
  character(50) :: fileName
  
  real(8) :: mu_k, q_in_post, k0_e_post_dummy, cross_section_norm, r, pi, theta
  real(8), dimension(4) :: k_e_post
  real(8), dimension(4) :: k_e_pre
  real(8) :: k,i_, cross_percent
  real(8) :: k0_e_pre
  real(8), dimension(3) :: stokes_pre
  real(8), dimension(3) :: stokes_post
  real(8), dimension(3) :: n
  real(8), dimension(7) :: q_in_pre
  real(8), dimension(11) :: theta_
  real(8) :: cos_chi_pre, sin_chi_pre, cos_chi_post, sin_chi_post, chi_pre, chi_post
  real(8) :: q_ein_pre, q_eex_pre, u_eex_pre, q_ein_post
  real(8) :: beta, gamma
  real(8) :: theta_post


  
  integer :: access

  iter_max = 100000 ! 100K
  !iter_max = 1000000  ! 1M
  !iter_max = 1000  ! 1K
  !iter_max = 10000 ! 10K

  pi = 4.0d0*atan(1.0d0)

  k_e_pre(1) =0.001d0
  k_e_pre(2) = k_e_pre(1)/sqrt(3.0d0)
  k_e_pre(3) = k_e_pre(1)/sqrt(3.0d0)
  k_e_pre(4) = k_e_pre(1)/sqrt(3.0d0)

  n(1:3) = (/1.0d0, 0.0d0, 0.0d0/)
  
  stokes_pre(1:3) = (/ 1.0d0, 0.0d0, 0.0d0/)
  
  q_in_pre(1:7) = (/ -0.8d0, -0.6d0, -0.4d0, 0.0d0, 0.4d0, 0.6d0, 0.8d0 /)
  
  theta_(1:11) = (/0.0d0, pi/10.0d0, pi/5.0d0, 3.0d0*pi/10.0d0, 4.0d0*pi/10.0d0, pi/2.0d0, &
       3.0d0*pi/5.0d0, 7.0d0*pi/10.0d0, 4.0d0*pi/5.0d0, 9.0d0*pi/10.0d0, pi/)
  
  beta = 0.0d0
  gamma = 1.0d0

  k0_e_pre  = k_e_pre(1)
  q_eex_pre = stokes_pre(1)
  u_eex_pre = stokes_pre(2)
  cross_percent = 0.0d0

  do i = 1,7
     q_ein_pre = q_in_pre(i)

      write(fileName, '("q^{e-in}"f11.6".csv")') q_ein_pre 
      fileName = trim(adjustl(fileName))

      fileNum  = i
      
      open(fileNum, file=fileName, status="replace")
      
      write(fileNum, '("iter_max = ",i6,", k_e_pre = ",f8.2,", q_exx_pre = ",f5.2,", u_exx_pre = ",f5.2)') &
        iter_max, k_e_pre(1),stokes_pre(1),stokes_pre(2) 
      write(fileNum, '( "theta, percentage ")')

      do j = 1,11

         k = 0.0d0
         theta = theta_(j)
      
      do iter = 1, iter_max
        do
           call sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k, theta_post, theta)

          ! if (k0_e_pre .lt. 1.0d-4) then  ! Thomson Scattering Limit
              call calc_thomson_diff(k0_e_pre, q_ein_pre, mu_k, &
                   k0_e_post_dummy, cross_section_norm)
          ! else  ! Compton Scattering
            !  call calc_kn_diff(k0_e_pre, q_ein_pre, mu_k, &
             !      k0_e_post_dummy, cross_section_norm)
           !end if
           
           k = k + 1.0d0

           ! call random_number(r)
           r = grnd()
                 if (cross_section_norm .gt. r) exit  ! Rejection method

!                 if(k .gt. 1000)then
!                    print '(f11.6, f11.6, f11.6, f11.6, f11.6)', k, theta, q_ein_pre, cross_section_norm, r
!                 end if
                 
              end do
              
           !   call update_stokes_parameters(stokes_pre, k_e_pre, k_e_post, mu_k, chi_pre, chi_post, stokes_post, q_ein_post)

              i_ = real(iter_max)
              cross_percent = i_/k
                 
    end do
    write(fileNum, '(f11.2,",",f11.6)') &
         theta_(j)*180/pi ,cross_percent
    
 end do
    
      close(fileNum)
   end do
   
   
end program main
