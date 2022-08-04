program main

  use scatter
  use mtmod  
  
  implicit none

    integer :: fileNum,iter,iter_max
    character(30) :: fileName

    real(8) :: z_total,x_total,y_total
    integer :: t,t_max, i,j,number
    real(8) :: mu_k, q_in_post, k0_e_post_dummy, cross_section_norm, r
    real(8), dimension(4) :: k_e_post, k_f_post
    real(8), dimension(4) :: k_e_pre, k_f_pre
    real(8) :: k0_e_pre
    real(8), dimension(3) :: stokes_pre
    real(8), dimension(3) :: stokes_post
    real(8), dimension(3) :: n
   ! real(16),dimension(100) :: num,mu_total
   ! real(16),dimension(100) :: q,u
   ! real(16),dimension(100) :: q_norm, u_norm
   ! real(16),dimension(100) :: mu_ave !bug? deback option check !FFLAGS = -O2 -g  -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
    real(8),dimension(100) :: num,mu_total
    real(8),dimension(100) :: q,u
    real(8),dimension(100) :: q_norm, u_norm
    real(8),dimension(100) :: mu_ave
    real(8) :: chi_pre, chi_post, cos_chi_post
    real(8) :: theta_pre, theta_post, phi_pre, phi_post
    real(8) :: q_ein_pre, q_eex_pre, u_eex_pre
    real(8) :: beta, gamma
    real(8) :: x,y,z,tau,mu,mu_
    real(8) :: pi,a,b,c,d
    real(8) :: p,o,q_i,u_i
    real(8) :: theta_e, phi_e
    real(8) :: depth,i_,j_
    real(8) :: theta_, phi_
    real(8) :: theta_s, s, mu_dummy
    real(8) :: count,number_

    integer :: access

    ! test1...use mtmod
    ! test2...use shokika
    ! test3...use use only thomson
    ! test4...change double precision
    
    t_max = 1000000 !!!3/10
    iter_max = 100 !3/10
    pi = 4.0d0*atan(1.0d0)
    depth = 6.0d0


    write(fileName, '("test_theta.csv")')
    fileName = trim(adjustl(fileName))

    fileNum  = 1

    open(fileNum, file=fileName, status="replace")

    write(fileNum, '("cos(i), q, u, number")')
    

    do i = 1,100
          q(i) = 0.0d0
          u(i) = 0.0d0
          mu_total(i) = 0.0d0
          num(i) = 0.0d0 ! 3/11 shokika
       end do
       

    do iter = 1, iter_max
     
       x_total = 0
       y_total = 0
       z_total = 0
   
       ! call random_number(s) !mtmod grnd
       s = grnd()
       theta_s = 2.0d0*pi*s

       stokes_pre(1) = cos(theta_s)
       stokes_pre(2) = sqrt(1 - stokes_pre(1)**2.0d0) !3/10 !test3

      ! stokes_pre(1) = 0.0d0
      ! stokes_pre(2) = 0.0d0 !nonpolarization
       stokes_pre(3) = 0.0d0 ! Initial condition
        mu_k = 1.0d0


       do t = 0, t_max
         ! call random_number(a)
         ! call random_number(b)
         ! call random_number(p)
         ! call random_number(o)

          a = grnd()
          b = grnd()
          p = grnd()
          o = grnd()

          theta_e = acos(p)
          phi_e = 2.0d0*pi*o

          n(1) = sin(theta_e)*cos(phi_e)
          n(2) = sin(theta_e)*sin(phi_e)
          n(3) = cos(theta_e)

          beta = b
          gamma = 1/sqrt(1 - beta*beta)
          tau = log(1/(1-a))
          
          if (t .eq. 0) then
            ! call random_number(c)
            ! call random_number(d)
             c = grnd()
             d = grnd()
         
             mu_= c
             !mu_ = 2.0d0*c-1.0d0
             theta_ = acos(mu_)
             phi_= 2.0d0*pi*d

             k_f_pre(1) = 0.001d0
             k_f_pre(2) = k_f_pre(1)*sin(theta_)*cos(phi_)
             k_f_pre(3) = k_f_pre(1)*sin(theta_)*sin(phi_)
             k_f_pre(4) = k_f_pre(1)*cos(theta_) ! Initial conditions

            ! k_f_pre(2) = 0.01d0*k_f_pre(1)
            ! k_f_pre(3) = 0.0d0
            ! k_f_pre(4) = sqrt(k_f_pre(1)**2.0d0-k_f_pre(2)**2.0d0)     
          
             k_f_post(1) = k_f_pre(1)
             k_f_post(2) = k_f_pre(2)
             k_f_post(3) = k_f_pre(3)
             k_f_post(4) = k_f_pre(4)

             stokes_post(1) = stokes_pre(1)
             stokes_post(2) = stokes_pre(2)         
             
          else

             q_eex_pre = stokes_pre(1)
             u_eex_pre = stokes_pre(2)

        

         !    call lorentz_trans(k_f_pre, k_e_pre, beta, gamma, n, k0_e_pre) ! lorentz

             k_e_pre(2) = k_f_pre(2)
             k_e_pre(3) = k_f_pre(3)
             k_e_pre(4) = k_f_pre(4)
             k_e_pre(1) = sqrt(k_e_pre(2)*k_e_pre(2)+k_e_pre(3)*k_e_pre(3)+k_e_pre(4)*k_e_pre(4)) !cleaning
             k0_e_pre = k_e_pre(1)
             

             
             do

                call sample_scattered_photons(k0_e_pre, k_e_pre, k_e_post, mu_k,theta_pre, theta_post, phi_pre, phi_post)

                if(phi_pre/=phi_pre)then
                   print '("chi=NaN")'
                end if
                

      
                
                call calc_chi(k_e_pre, k_e_post, mu_k, &
                     chi_pre, chi_post, theta_pre, theta_post, phi_pre, phi_post, cos_chi_post)

                call rotate_q_chi(q_eex_pre, u_eex_pre, q_ein_pre, chi_pre, chi_post)
                                           

               ! if (k0_e_pre .lt. 1.0d-4) then  ! Thomson Scattering Limit
                   call calc_thomson_diff(k0_e_pre, q_ein_pre, mu_k, &
                        k0_e_post_dummy, cross_section_norm) !3/11 test3
               ! else  ! Compton Scattering
                !   call calc_kn_diff(k0_e_pre, q_ein_pre, mu_k, &
                 !       k0_e_post_dummy, cross_section_norm)
               ! end if

                r = grnd()


                if (cross_section_norm .gt. r) exit !rejection


             end do
             

             call update_stokes(stokes_pre, k_e_pre, k_e_post, mu_k, chi_pre, chi_post, stokes_post)

          !   call lorentz_inverse(k_f_post, k_e_post, beta, gamma, n) ! lorentz method


             k_f_post(1) = k_e_post(1)
             k_f_post(2) = k_e_post(2)
             k_f_post(3) = k_e_post(3)
             k_f_post(4) = k_e_post(4)
             
          end if


          x = k_f_post(2)*tau/k_f_post(1)
          y = k_f_post(3)*tau/k_f_post(1)
          z = k_f_post(4)*tau/k_f_post(1)
          

       x_total = x_total + x
       y_total = y_total + y
       z_total = z_total + z

       mu = k_f_post(4)/sqrt(k_f_post(2)**2.0d0 + k_f_post(3)**2.0d0 + k_f_post(4)**2.0d0)

       stokes_pre(1) = stokes_post(1)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2))
       stokes_pre(2) = stokes_post(2)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2)) ! cleaning of stokes parameters

       if(stokes_post(1) .eq. 0.0d0 .and. stokes_post(2) .eq. 0.0d0)then
          stokes_pre(1) = 0.0d0
          stokes_pre(2) = 0.0d0
       end if
       

     !  stokes_pre(1) = stokes_post(1)
       !  stokes_pre(2) = sqrt(1-stokes_pre(1)*stokes_pre(1))  !another cleaning

       mu_in = (k_f_pre(2)*k_f_post(2) + k_f_pre(3)*k_f_post(3) + k_f_pre(4)*k_f_post(4))/(k_f_pre(1)*k_f_post(1))

       k_f_pre(1) = k_f_post(1)
       k_f_pre(2) = k_f_post(2)
       k_f_pre(3) = k_f_post(3)
       k_f_pre(4) = k_f_post(4)

       if(abs(z_total) .gt. depth) exit !Escape Conditions
     ! if(z_total .lt. 0) exit
      
    end do
    
    
    
    if(abs(z_total) .gt. depth)then
       do i = 1,100
          i_ = real(i)
          if(abs(mu) .gt. (i_-1)*0.01d0 .and. abs(mu) .lt. i_*0.01d0)then
          q(i) = q(i) + stokes_post(1)
          u(i) = u(i) + stokes_post(2)
          mu_total(i) = mu_total(i) + abs(mu)
          num(i) = num(i) + 1.0d0 ! 3/11 shokika
       end if
       
       end do         
    end if

    count = real(iter)
    do number = 1,10
       number_ = real(number)
       if (count .eq. number_*(10**7.0d0))then
          print '(i9)',iter
       end if
    end do
    
 end do
 
 
 do j = 1,100
    q_norm(j) = q(j)/num(j)
    u_norm(j) = u(j)/num(j)
    mu_ave(j) = mu_total(j)/num(j)
    write(fileNum, '(f15.6,",",f15.6,",",f15.6,",",f15.0)')&
         mu_ave(j),q_norm(j),u_norm(j),num(j)
 end do
 
 
        
 close(fileNum)


    
end program main
       

          

    

    
