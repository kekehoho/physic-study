program main

  use scatter
  use mtmod  
  
  implicit none

    integer :: fileNum,iter,iter_max
    character(50) :: fileName

    real(8) :: z_total,x_total,y_total
    integer :: t,t_max, i,j,number, m
    real(8) :: mu_k, q_in_post, k0_e_post_dummy, cross_section_norm, r
    real(8), dimension(4) :: k_e_post, k_f_post
    real(8), dimension(4) :: k_e_pre, k_f_pre
    real(8) :: k0_e_pre
    real(8), dimension(3) :: stokes_pre
    real(8), dimension(3) :: stokes_post
    real(8), dimension(3) :: n
    real(8),dimension(16,100) :: num,mu_total
    real(8),dimension(16,100) :: q,u
    real(8),dimension(16,100) :: q_norm, u_norm
    real(8),dimension(16,100) :: mu_ave
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
    real(8) :: theta_s, s
    real(8) :: count,number_
    real(8) :: section,random
    real(8) :: delta_tau, delta_section,tau_slab
    real(8) :: x_p,y_p,z_p,s_p
    real(8) :: z_delta,tent_section
    real(8) :: alpha,alpha_0
    real(8) :: x_delta,y_delta ! shuoud be fixed
    real(8) :: radius,radius_xy
    real(8) :: mu_pre,t_,t_total,t_ave
    integer :: times
    real(8) :: times_
    integer :: access

    !bug? deback option check !FFLAGS = -O2 -g  -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
    
     t_max = 100000000
   ! t_max = 2 ! test to be fixed 
    iter_max = 100000000 !should be fixed 
    pi = 4.0d0*atan(1.0d0)
    depth = 6.0d0 !*10.0d0**(163.0d0/40.0d0) ! shuold be fixed
    tau_slab = 6.0d0

    tent_section = 0.1d0
    alpha_0 = (tau_slab/(depth*(1.0d0-exp(-1.0d0))))!*10.0d0**(-133.0d0/40.0d0)


    ! write(fileName, '("test_mu_pre_depth_6_q.csv")')
    ! write(fileName, '("lv1_x_layer.csv")') 
    t_total = 0.0d0
    

    do i = 1,16
       do j = 1,100
          q(i,j) = 0.0d0
          u(i,j) = 0.0d0
          mu_total(i,j) = 0.0d0
          num(i,j) = 0.0d0 
       end do
    end do
    
       

    do iter = 1, iter_max
     
       x_total = 0.0d0
       y_total = 0.0d0
       z_total = 0.0d0
   
       s = grnd()
       theta_s = 2.0d0*pi*s

       stokes_pre(1) = cos(theta_s)
       stokes_pre(2) = sin(theta_s) 

      ! stokes_pre(1) = 0.0d0
      ! stokes_pre(2) = 0.0d0 !nonpolarization
       stokes_pre(3) = 0.0d0 ! Initial condition
        mu_k = 1.0d0


       do t = 0, t_max
        
         ! b = grnd()
         ! p = grnd()
         ! o = grnd()

         ! theta_e = acos(p)
         ! phi_e = 2.0d0*pi*o

         ! n(1) = sin(theta_e)*cos(phi_e)
         ! n(2) = sin(theta_e)*sin(phi_e)
         ! n(3) = cos(theta_e)

         ! beta = b
         ! gamma = 1/sqrt(1 - beta*beta)
          
          if (t .eq. 0) then

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
          
          section = 0.0d0

          do 
              random = grnd()

              alpha = alpha_0*exp(-abs(z_total)/depth)

             tau = tent_section*alpha

              if(exp(-tau) .lt. random)exit
             
                section = section + tent_section

                z = k_f_post(4)*tent_section/k_f_post(1) ! shuoud be fixed                                                                                        
               z_total = z_total + z ! shuoud be fixed 

               if(abs(z_total) .gt. depth)exit
          end do

        ! delta_tau = log(1/(1-random))
         delta_tau = log(1/random)
          
         delta_section = delta_tau/alpha

         s_p = section + delta_section
          
          x = k_f_post(2)*s_p/k_f_post(1)! shuoud be fixed
          y = k_f_post(3)*s_p/k_f_post(1)! shuoud be fixed 
          z_delta = k_f_post(4)*delta_section/k_f_post(1)! shuoud be fixed
          

       z_total = z_total + z_delta ! shuoud be fixed
       y_total = y_total + y
       x_total = x_total + x ! shuoud be fixed

       mu = k_f_post(4)/sqrt(k_f_post(2)**2.0d0 + k_f_post(3)**2.0d0 + k_f_post(4)**2.0d0)

       t_ = real(t)
       
       do times = 1,16
          times_ = real(times)

          if((times_-1.0d0)*5.0d0 .eq. t_)then
             do i = 1,100
                i_ = real(i)
          if(abs(mu) .ge. (i_- 1.0d0)*0.01d0 .and. abs(mu) .lt. i_*0.01d0)then
             q(times,i) = q(times,i) + stokes_post(1)
             u(times,i) = u(times,i) + stokes_post(2)
             mu_total(times,i) = mu_total(times,i) + abs(mu)
             num(times,i) = num(times,i) + 1.0d0
          end if
          
       end do
    end if
    
 end do
 
       

       if(abs(z_total) .gt. depth) exit !Escape Conditions
       
       stokes_pre(1) = stokes_post(1)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2))
       stokes_pre(2) = stokes_post(2)/sqrt(stokes_post(1)*stokes_post(1)+stokes_post(2)*stokes_post(2)) ! cleaning of stokes parameters

       if(stokes_post(1) .eq. 0.0d0 .and. stokes_post(2) .eq. 0.0d0)then
          stokes_pre(1) = 0.0d0
          stokes_pre(2) = 0.0d0
       end if
       

     !  stokes_pre(1) = stokes_post(1)
       !  stokes_pre(2) = sqrt(1-stokes_pre(1)*stokes_pre(1))  !another cleaning
       k_f_pre(1) = k_f_post(1)
       k_f_pre(2) = k_f_post(2)
       k_f_pre(3) = k_f_post(3)
       k_f_pre(4) = k_f_post(4)
     
    end do
    
    t_ = real(t)
    t_total = t_total + t_
    
 !      do i = 1,100
 !         i_ = real(i)
 !         if(abs(mu) .ge. (i_- 1.0d0)*0.01d0 .and. abs(mu) .lt. i_*0.01d0)then
 !         q(i) = q(i) + stokes_post(1)
 !         u(i) = u(i) + stokes_post(2)
 !         mu_total(i) = mu_total(i) + abs(mu)
 !         num(i) = num(i) + 1.0d0 
 !      end if
 !   end do
    


 count = real(iter)
    do number = 1,10
       number_ = real(number)
       if (count .eq. number_*(10**7.0d0))then
          print '(i9)',iter
       end if
    end do
 end do

 t_ave = t_total/iter_max

 print'(f11.2)',t_ave
 

 do j = 1,16
   
    write(fileName, '("test_",i2,"_times.csv")') (j-1)*5
    
    fileName = trim(adjustl(fileName))

    fileNum  = j

    open(fileNum, file=fileName, status="replace")

    write(fileNum, '("mu, q, u, N")')

    do i = 1,100
    
       q_norm(j,i) = q(j,i)/num(j,i)
       u_norm(j,i) = u(j,i)/num(j,i)
       mu_ave(j,i) = mu_total(j,i)/num(j,i)
       
    write(fileNum, '(f15.6,",",f15.6,",",f15.6,",",f15.0)')&
         mu_ave(j,i) ,q_norm(j,i),u_norm(j,i),num(j,i)
 end do
end do


 
 
! do j = 1,100
!    q_norm(j) = q(j)/num(j)
!    u_norm(j) = u(j)/num(j)
!    mu_ave(j) = mu_total(j)/num(j)
!    i_ = real(j)
!    write(fileNum, '(f15.6,",",f15.6,",",f15.6,",",f15.0)')&
!         mu_ave(j) ,q_norm(j),u_norm(j),num(j)
! end do
 
 
        
 close(fileNum)


    
end program main
       

          

    

    
