program main

  use scatter
  use mtmod

  implicit none

  integer :: iter, iter_max,fileNum,i
  character(50) :: fileName

  real(8) :: r,a
  real(8) :: k0_e_pre, mu_k, theta_pre, theta_post, phi_pre, phi_post
  real(8),dimension(4) :: k_e_pre, k_e_post

   integer :: access


  iter_max = 100

  do i = 1,5

  write(fileName, '("test_mt"i1".csv")')i
  fileName = trim(adjustl(fileName))

  fileNum  = 1

  open(fileNum, file=fileName, status="replace")

  write(fileNum, '( "iter, ramdom number")')

  k_e_pre(1) =0.001d0
  k_e_pre(2) = k_e_pre(1)/sqrt(3.0d0)
  k_e_pre(3) = k_e_pre(1)/sqrt(3.0d0)
  k_e_pre(4) = k_e_pre(1)/sqrt(3.0d0)

  k0_e_pre = k_e_pre(1)

  do iter = 1,iter_max

     r = grnd()
     a = grnd()

     call sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k, theta_pre, theta_post, phi_pre, phi_post)

     write(fileNum, '(i5,",",f11.6,",",f11.6,",",f11.6,",",f11.6)') &
          iter ,r, a, theta_post, phi_post

  end do

  close(fileNum)

  end do


end program main

  
