! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


module soap_turbo_radial

  use soap_turbo_functions

  contains

!**************************************************************************
!
! This function returns the normalization coefficients for the polynomial
! basis functions used to construct the orthonormal basis.
!
  function N_a(rcut, a)
    implicit none

    integer :: a, b
    real*8 :: rcut, N_a

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    N_a = dsqrt( rcut / dfloat(b) )
!   *********************************************

  return
  end function
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and the saparable Gaussian representation for the
! atomic sites.
!
  subroutine get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                     rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                     amplitude_scaling, nf, W, scaling_mode, mask, &
                                                     radial_enhancement, do_derivatives, do_central, &
                                                     central_weight, exp_coeff, exp_coeff_der)
!   Expansion coefficients using the polynomial basis with smooth filter

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: n, i, j, k, alpha_max_der
    real*8 :: I_n, I_np1, I_np2, pi, sq2, rj, N_n, N_np1, N_np2, C1, C2, dr, s2, sf2
    real*8 :: W(:,:)
    real*8 :: atom_sigma_f, rj_f
!   Results will be stored in exp_coeff, which is an array of dimension
!   (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_temp1(:), exp_coeff_temp2(:), exp_coeff_der_temp(:)
    logical, save :: print_basis = .false., print_message = .true.
    real*8 :: denom, der_sjf_rj, der_rjf_rj, amplitude_der, pref_f, der_pref_f

!   For this basis numerical instabilities in the orthonormal basis construction develop
!   above alpha_max = 7. These are relatively small up to alpha_max = 10 and become
!   catastrophic at alpha_max = 12.
    if( alpha_max > 10 )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "ERROR: Due to numerical instabilities in the basis construction for the poly3   <---- ERROR"
      write(*,*) "basis, it is strongly recommended not to exceed alpha_max = 7. For"
      write(*,*) "alpha_max > 10 the instabilities are too large to proceed. You can do your own"
      write(*,*) "testing for 7 < alpha_max < 11, the results might still be useful within that"
      write(*,*) "range. Note that the poly3gauss basis allows you to add one extra basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
      stop
    else if( alpha_max > 7 .and. print_message )then
      print_message = .false.
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "WARNING: Due to numerical instabilities in the basis construction for the poly3 <---- WARNING"
      write(*,*) "basis, it is strongly recommended not to exceed alpha_max = 7. For"
      write(*,*) "alpha_max > 10 this warning will turn into an error. You can do your own"
      write(*,*) "testing for 7 < alpha_max < 11, the results might still be useful within that"
      write(*,*) "range. Note that the poly3gauss basis allows you to add one extra basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
    end if

!   If the user requests derivatives, we need to get the expansion coefficients up to
!   alpha_max + 2
    if( do_derivatives )then
      alpha_max_der = alpha_max + 2
    else
      alpha_max_der = alpha_max
    end if
    allocate( exp_coeff_temp1(1:alpha_max_der) )
    allocate( exp_coeff_temp2(1:alpha_max_der) )
    allocate( exp_coeff_der_temp(1:alpha_max) )


!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

    pi = dacos(-1.d0)
    sq2 = dsqrt(2.d0)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if
!
!   **************** New basis ******************
!   Redefine all the distances by dividing them by rcut_hard
!   We do this to avoid numerical instability when the value
!   of alpha is too high
!
!    rcut_soft = rcut_soft_in
!    rcut_hard = rcut_hard_in
!    atom_sigma = atom_sigma_in
!    dr = rcut_hard - rcut_soft
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
    dr = 1.d0 - rcut_soft_in/rcut_hard_in
!   *********************************************
!


    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
!       Check if we need to do the central atom
        if( j == 1 .and. .not. do_central )then
          cycle
        end if
        if( rjs_in(k) < rcut_hard_in .and. mask(k) )then
          exp_coeff_temp1 = 0.d0
          exp_coeff_temp2 = 0.d0
          exp_coeff_der_temp = 0.d0
!   **************** New basis ******************
!          rj = rjs_in(k)
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
!         We leave this here because in the future atom_sigma will be rj dependent
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitude = 1.d0 / atom_sigma_scaled
              amplitude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              amplitude = 0.d0
              amplitude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                amplitude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              else
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                amplitude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( j == 1 )then
            amplitude = central_weight * amplitude
            amplitude_der = central_weight * amplitude_der
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitude_der = amplitude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            amplitude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            amplitude = amplitude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            amplitude_der = amplitude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            amplitude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            amplitude = amplitude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if       
!         We have the recursion series starting at n = 0, which means alpha = -2
!         However, we only need to save the expansion coefficients for alpha >= 1
!         This is I_-1
          I_n = 0.d0
          N_n = 1.d0
!         This is I_0
          N_np1 = N_a(rcut_hard, -2)
          I_np1 = dsqrt(pi/2.d0) * atom_sigma_scaled * ( derf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - &
                                                  derf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1
!         Speed up the computation of these coefficients
          if( rcut_hard_in == rcut_soft_in )then
            C1 = 0.d0
          else
            C1 = s2 / dr * dexp(-0.5d0 * (rcut_soft - rj)**2 / s2)
          end if
          C2 = s2 / rcut_hard * dexp(-0.5d0 * rj**2 / s2)
          do n = -1, alpha_max_der
            C1 = C1 * dr
            C2 = C2 * rcut_hard
            N_np2 = N_a(rcut_hard, n)
!           This is I_alpha
            I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                    - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 &
                    + C1 / N_np2 &
                    - C2 / N_np2
            if(n > 0)then
              exp_coeff_temp1(n) = I_np2
            end if
            N_n = N_np1
            N_np1 = N_np2
            I_n = I_np1
            I_np1 = I_np2
          end do
!         Compute the contribution to the derivative for this part (excludes the amplitude)
          if( do_derivatives )then
            do n = 1, alpha_max
              exp_coeff_der_temp(n) = (rj - rcut_hard)/s2 * ( atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled &
                                  - 1.d0 ) * exp_coeff_temp1(n) + &
                                 rcut_hard*N_a(rcut_hard, n+1)/s2/N_a(rcut_hard, n) * ( 2.d0 * atom_sigma_scaling &
                                  * (rj - rcut_hard) / atom_sigma_scaled - 1.d0 ) * exp_coeff_temp1(n+1) + &
                                 atom_sigma_scaling*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_scaled**3/ &
                                  N_a(rcut_hard, n) * exp_coeff_temp1(n+2)
            end do
          end if
!         If the atom is less than 4 standard deviations away from the soft cutoff, we add
!         also this correction to the integrals. This corresponds to a Gaussian filter. We
!         integrate between rcut_soft and rcut_hard in this case
!
!         This explicit ".true." or ".false." logical statement is there for debugging. For
!         regular code use it can be set to .false.
          if( .false. .or. (rcut_soft - rj) < 4.d0*atom_sigma_scaled )then
            atom_sigma_f = atom_sigma_scaled * dr / nf / dsqrt(s2 + dr**2/nf**2)
            rj_f = (s2 * rcut_soft + dr**2/nf**2*rj) / (s2 + dr**2/nf**2)
!           We leave this here because in the future atom_sigma will be rj dependent
            sf2 = atom_sigma_f**2
!           The products of two Gaussians is a Gaussian, but we need to add a prefactor
            pref_f = dexp( -0.5d0 * (rcut_soft-rj)**2 / ( s2 + dr**2/nf**2 ) )
!           We have the recursion series starting at n = 0, which means alpha = -2
!           However, we only need to save the expansion coefficients for alpha >= 1
!           This is I_-1
            I_n = 0.d0
            N_n = 1.d0
!           This is I_0
            N_np1 = N_a(rcut_hard, -2)
            I_np1 = dsqrt(pi/2.d0) * atom_sigma_f * ( derf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - &
                                                      derf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1
!           Speed up the computation of these coefficients
            C2 = sf2 / dr * dexp(-0.5d0 * (rcut_soft - rj_f)**2 / sf2)
            do n = -1, alpha_max_der
              C2 = C2 * dr
              N_np2 = N_a(rcut_hard, n)
!             This is I_alpha
              I_np2 = sf2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                      - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 &
                      - C2 / N_np2
              if(n > 0)then
                exp_coeff_temp2(n) = I_np2
              end if
              N_n = N_np1
              N_np1 = N_np2
              I_n = I_np1
              I_np1 = I_np2
            end do
!           Compute the contribution to the derivative for this part (excludes the amplitude)
            if( do_derivatives )then
              denom = s2 + dr**2/nf**2
              der_pref_f = pref_f * ( (rcut_soft - rj) / denom + (rcut_soft - rj)**2 / denom**2 * &
                                      atom_sigma_scaled * atom_sigma_scaling )
              der_rjf_rj = (2.d0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + dr**2/nf**2) / denom &
                           - (s2*rcut_soft + dr**2/nf**2 * rj) * 2.d0 * atom_sigma_scaled * &
                             atom_sigma_scaling / denom**2
              der_sjf_rj = atom_sigma_scaling * dr/nf / dsqrt(denom) * (1.d0 - atom_sigma_scaled**2/denom)
              do n = 1, alpha_max
                exp_coeff_der_temp(n) = exp_coeff_der_temp(n) + &
                                      pref_f * ( &
                                      (rj_f - rcut_hard)/sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f &
                                        - der_rjf_rj ) * exp_coeff_temp2(n) + &
                                      rcut_hard*N_a(rcut_hard, n+1)/sf2/N_a(rcut_hard, n) * ( 2.d0 * der_sjf_rj &
                                        * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj ) * exp_coeff_temp2(n+1) + &
                                      der_sjf_rj*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_f**3/ &
                                      N_a(rcut_hard, n) * exp_coeff_temp2(n+2) ) + &
                                      der_pref_f * &
                                      exp_coeff_temp2(n)
              end do
            end if
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          if( do_derivatives )then
            exp_coeff_der_temp(1:alpha_max) = amplitude * exp_coeff_der_temp(1:alpha_max) + amplitude_der * &
                                              (exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max))
            exp_coeff_der(1:alpha_max, k) = matmul( W, exp_coeff_der_temp(1:alpha_max) )
          end if
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max) )
        end if
      end do
    end do




!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   *********************************************

!   This is for debugging
    if( .false. )then
      open(10, file="radial_expansion_coefficients.dat", status="unknown", access="append")
      write(10,*) exp_coeff(1:alpha_max, 2)
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_derivatives.dat", status="unknown", access="append")
        write(10,*) exp_coeff_der(1:alpha_max, 2)
        close(10)
      end if
     end if

    deallocate( exp_coeff_temp1, exp_coeff_temp2, exp_coeff_der_temp )

  return
  end subroutine
!**************************************************************************










!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set augmented with a Gaussian function at the origin
! and the saparable Gaussian representation for the atomic sites.
!
  subroutine get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                          rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                          amplitude_scaling, nf, W, scaling_mode, mask, &
                                                          radial_enhancement, do_derivatives, exp_coeff, &
                                                          exp_coeff_der)
!   Expansion coefficients using the polynomial basis with smooth filter plus a Gaussian centered at the origin
!
!   TRY OUT: Check for very small numbers (that could lead to negative numbers) and then truncate them
!   to zero
!
    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives
    character(*), intent(in) :: scaling_mode
!
    integer :: n, i, j, k, alpha_max_der
    real*8 :: I_n, I_np1, I_np2, pi, sq2, rj, N_n, N_np1, N_np2, C1, C2, dr, s2, sf2
    real*8 :: W(:,:)
    real*8 :: atom_sigma_f, rj_f, N_gauss
!   Results will be stored in exp_coeff, which is an array of dimension
!   (n_sites, n_neigh_max, alpha_max)
    real*8 :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_temp1(:), exp_coeff_temp2(:), exp_coeff_der_temp(:)
    logical, save :: print_basis = .false., print_message = .true.
    real*8 :: denom, der_sjf_rj, der_rjf_rj, amplitude_der, pref_f, der_pref_f, sigma_star

!   For this basis numerical instabilities in the orthonormal basis construction develop
!   above alpha_max = 8. These are relatively small up to alpha_max = 11 and become
!   catastrophic at alpha_max = 13.
    if( alpha_max > 11 )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "ERROR: Due to numerical instabilities in the basis construction for the         <---- ERROR"
      write(*,*) "poly3gauss basis, it is strongly recommended not to exceed alpha_max = 8. For"
      write(*,*) "alpha_max > 11 the instabilities are too large to proceed. You can do your own"
      write(*,*) "testing for 8 < alpha_max < 12, the results might still be useful within that"
      write(*,*) "range. Note that the poly3 basis allows you to add one *less* basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
      stop
    else if( alpha_max > 8 .and. print_message )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "WARNING: Due to numerical instabilities in the basis construction for the       <---- WARNING"
      write(*,*) "poly3gauss basis, it is strongly recommended not to exceed alpha_max = 8. For"
      write(*,*) "alpha_max > 11 this warning will turn into an error. You can do your own"
      write(*,*) "testing for 8 < alpha_max < 12, the results might still be useful within that"
      write(*,*) "range. Note that the poly3 basis allows you to add one *less* basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
    end if

!   If the user requests derivatives, we need to get the expansion coefficients up to
!   alpha_max - 1 + 2. The "-1" is there because the Gaussian basis at the origin does not
!   participate in the calculation of the derivatives for the polynomial basis functions
    if( do_derivatives )then
      alpha_max_der = alpha_max + 2
    else
      alpha_max_der = alpha_max
    end if
    allocate( exp_coeff_temp1(1:alpha_max_der) )
    allocate( exp_coeff_temp2(1:alpha_max_der) )
    allocate( exp_coeff_der_temp(1:alpha_max) )

!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
!    if( .false. .and. print_basis )then
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      write(*,*) "g(x,s) = exp(-0.5*x**2/s**2) * sqrt(2./s) / pi**0.25 "
      do j = 1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max-1
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,"(A,F16.10,A,E16.8,A)",advance="no") "g(x,", atom_sigma_in, ") *", W(j,alpha_max), "+"
        write(*,*) "0."
      end do
    end if

    pi = dacos(-1.d0)
    sq2 = dsqrt(2.d0)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if
!
!   **************** New basis ******************
!   Redefine all the distances by dividing them by rcut_hard
!   We do this to avoid numerical instability when the value
!   of alpha is too high
!
!    rcut_soft = rcut_soft_in
!    rcut_hard = rcut_hard_in
!    atom_sigma = atom_sigma_in
!    dr = rcut_hard - rcut_soft
!    N_gauss = dsqrt(2.d0/atom_sigma_in) / pi**0.25
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
    dr = 1.d0 - rcut_soft_in/rcut_hard_in
    N_gauss = dsqrt(2.d0/atom_sigma) / pi**0.25d0
!   *********************************************
!
    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
!       IMPORTANT: for this basis, we skip i itself, which is neighbor number 1
        if( j == 1 )then
          cycle
        end if
        if( rjs_in(k) < rcut_hard_in .and. mask(k) )then
          exp_coeff_temp1 = 0.d0
          exp_coeff_temp2 = 0.d0
          exp_coeff_der_temp = 0.d0
!   **************** New basis ******************
!          rj = rjs_in(k)
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
!         We leave this here because in the future atom_sigma will be rj dependent
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitude = 1.d0 / atom_sigma_scaled
              amplitude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              amplitude = 0.d0
              amplitude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                amplitude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              else
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                amplitude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              end if
            end if
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitude_der = amplitude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            amplitude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            amplitude = amplitude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            amplitude_der = amplitude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            amplitude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            amplitude = amplitude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if
!         We have the recursion series starting at n = 0, which means alpha = -2
!         However, we only need to save the expansion coefficients for alpha >= 1
!         This is I_-1
          I_n = 0.d0
          N_n = 1.d0
!         This is I_0
          N_np1 = N_a(rcut_hard, -2)
          I_np1 = dsqrt(pi/2.d0) * atom_sigma_scaled * ( derf( (rcut_soft-rj)/sq2/atom_sigma_scaled ) - &
                                                  derf( (-rj)/sq2/atom_sigma_scaled ) ) / N_np1
!         Speed up the computation of these coefficients
          if( rcut_hard_in == rcut_soft_in )then
            C1 = 0.d0
          else
            C1 = s2 / dr * dexp(-0.5d0 * (rcut_soft - rj)**2 / s2)
          end if
          C2 = s2 / rcut_hard * dexp(-0.5d0 * rj**2 / s2)
!         This is different wrt the regular polynomial basis, we only go up to alpha_max-1
          do n = -1, alpha_max_der-1
            C1 = C1 * dr
            C2 = C2 * rcut_hard
            N_np2 = N_a(rcut_hard, n)
!           This is I_alpha
            I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                    - N_np1 * (rj - rcut_hard) / N_np2 * I_np1 &
                    + C1 / N_np2 &
                    - C2 / N_np2
            if(n > 0)then
              exp_coeff_temp1(n) = I_np2
            end if
            N_n = N_np1
            N_np1 = N_np2
            I_n = I_np1
            I_np1 = I_np2
          end do
!         Compute the contribution to the derivative for this part (excludes the amplitude)
          if( do_derivatives )then
            do n = 1, alpha_max-1
              exp_coeff_der_temp(n) = (rj - rcut_hard)/s2 * ( atom_sigma_scaling * (rj - rcut_hard) / atom_sigma_scaled &
                                  - 1.d0 ) * exp_coeff_temp1(n) + &
                                 rcut_hard*N_a(rcut_hard, n+1)/s2/N_a(rcut_hard, n) * ( 2.d0 * atom_sigma_scaling &
                                  * (rj - rcut_hard) / atom_sigma_scaled - 1.d0 ) * exp_coeff_temp1(n+1) + &
                                 atom_sigma_scaling*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_scaled**3/ &
                                  N_a(rcut_hard, n) * exp_coeff_temp1(n+2)
            end do
          end if
!         If the atom is less than 4 standard deviations away from the soft cutoff, we add
!         also this correction to the integrals. This corresponds to a Gaussian filter. We
!         integrate between rcut_soft and rcut_hard in this case
!
!         This explicit ".true." or ".false." logical statement is there for debugging. For
!         regular code use it can be set to .false.
          if( .false. .or. (rcut_soft - rj) < 4.d0*atom_sigma_scaled )then
            atom_sigma_f = atom_sigma_scaled * dr / nf / dsqrt(s2 + dr**2/nf**2)
            rj_f = (s2 * rcut_soft + dr**2/nf**2*rj) / (s2 + dr**2/nf**2)
!           We leave this here because in the future atom_sigma will be rj dependent
            sf2 = atom_sigma_f**2
!           The products of two Gaussians is a Gaussian, but we need to add a prefactor
            pref_f = dexp( -0.5d0 * (rcut_soft-rj)**2 / ( s2 + dr**2/nf**2 ) )
!           We have the recursion series starting at n = 0, which means alpha = -2
!           However, we only need to save the expansion coefficients for alpha >= 1
!           This is I_-1
            I_n = 0.d0
            N_n = 1.d0
!           This is I_0
            N_np1 = N_a(rcut_hard, -2)
            I_np1 = dsqrt(pi/2.d0) * atom_sigma_f * ( derf( (rcut_hard-rj_f)/sq2/atom_sigma_f ) - &
                                                      derf( (rcut_soft-rj_f)/sq2/atom_sigma_f ) ) / N_np1
!           Speed up the computation of these coefficients
            C2 = sf2 / dr * dexp(-0.5d0 * (rcut_soft - rj_f)**2 / sf2)
!           This is different wrt the regular polynomial basis, we only go up to alpha_max-1
            do n = -1, alpha_max_der-1
              C2 = C2 * dr
              N_np2 = N_a(rcut_hard, n)
!             This is I_alpha
              I_np2 = sf2 * dfloat(n+1) * N_n/ N_np2 * I_n &
                      - N_np1 * (rj_f - rcut_hard) / N_np2 * I_np1 &
                      - C2 / N_np2
              if(n > 0)then
                exp_coeff_temp2(n) = exp_coeff_temp2(n) + I_np2
              end if
              N_n = N_np1
              N_np1 = N_np2
              I_n = I_np1
              I_np1 = I_np2
            end do
!           Compute the contribution to the derivative for this part (excludes the amplitude)
            if( do_derivatives )then
              denom = s2 + dr**2/nf**2
              der_pref_f = pref_f * ( (rcut_soft - rj) / denom + (rcut_soft - rj)**2 / denom**2 * &
                                      atom_sigma_scaled * atom_sigma_scaling )
              der_rjf_rj = (2.d0*atom_sigma_scaled*rcut_soft*atom_sigma_scaling + dr**2/nf**2) / denom &
                           - (s2*rcut_soft + dr**2/nf**2 * rj) * 2.d0 * atom_sigma_scaled * &
                             atom_sigma_scaling / denom**2
              der_sjf_rj = atom_sigma_scaling * dr/nf / dsqrt(denom) * (1.d0 - atom_sigma_scaled**2/denom)
              do n = 1, alpha_max-1
                exp_coeff_der_temp(n) = exp_coeff_der_temp(n) + &
                                      pref_f * ( &
                                      (rj_f - rcut_hard)/sf2 * ( der_sjf_rj * (rj_f - rcut_hard) / atom_sigma_f &
                                        - der_rjf_rj ) * exp_coeff_temp2(n) + &
                                      rcut_hard*N_a(rcut_hard, n+1)/sf2/N_a(rcut_hard, n) * ( 2.d0 * der_sjf_rj &
                                        * (rj_f - rcut_hard) / atom_sigma_f - der_rjf_rj ) * exp_coeff_temp2(n+1) + &
                                      der_sjf_rj*rcut_hard**2*N_a(rcut_hard, n+2)/atom_sigma_f**3/ &
                                      N_a(rcut_hard, n) * exp_coeff_temp2(n+2) ) + &
                                      der_pref_f * &
                                      exp_coeff_temp2(n)
              end do
            end if
          end if
!         Now we obtain the overlap integral between the atomic density and the Gaussian centered at the origin
!         We assume that at the soft cutoff the Gaussian basis function is approx. zero, and we only
!         compute the overlap coefficient if the atom is close to the origin
!         Note that the sigma for the Gaussian basis function and the atom's sigma are not the same if there is
!         sigma scaling
!         We reset these to zero (temp2 will stay zero since it's the filter one, and we do not apply filter to
!         the overlap with the central Gaussian. This should be a good approximation if rcut_soft >= 4 * atom_sigma
          exp_coeff_temp1(alpha_max) = 0.d0
          exp_coeff_temp2(alpha_max) = 0.d0
          if( .false. .or. rj < 4.d0*(atom_sigma+atom_sigma_scaled) )then
              sigma_star = dsqrt(atom_sigma**2 + s2)
              exp_coeff_temp1(alpha_max) = dexp(- 0.5d0 * rj**2 / sigma_star**2 ) * dsqrt(pi/2.d0) * &
                                        atom_sigma_scaled*atom_sigma / sigma_star * ( 1.d0 &
                                        + derf(atom_sigma/atom_sigma_scaled*rj/sq2/sigma_star) )* N_gauss
            if( do_derivatives )then
              exp_coeff_der_temp(alpha_max) = ( rj**2 * atom_sigma_scaling / atom_sigma_scaled**3 - rj/sigma_star**2 + &
                                                atom_sigma_scaling*rj**2*atom_sigma**4/atom_sigma_scaled**3/sigma_star**4 + &
                                                atom_sigma_scaling*atom_sigma**2/atom_sigma_scaled/sigma_star**2 - &
                                                2.d0*rj**2*atom_sigma_scaling*atom_sigma**2/atom_sigma_scaled**3/sigma_star**2 &
                                              ) * exp_coeff_temp1(alpha_max) + &
                                              (1./s2 - 2.d0*rj*atom_sigma_scaling/atom_sigma_scaled**3) * s2 * atom_sigma**2 / &
                                                sigma_star**2 * dsqrt(2.d0/atom_sigma) / pi**0.25d0 * &
                                                dexp(-0.5d0 * rj**2 / sigma_star**2 * (1.d0 + atom_sigma**2 / s2) ) + &
                                              dsqrt(2.d0/atom_sigma) / pi**0.25d0 * dexp(-0.5d0 * rj**2 / sigma_star**2 * &
                                                (1.d0 + atom_sigma**2 / s2) ) * atom_sigma_scaling / atom_sigma_scaled * &
                                                rj*atom_sigma**4/sigma_star**4
            end if
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          if( do_derivatives )then
            exp_coeff_der_temp(1:alpha_max) = amplitude * exp_coeff_der_temp(1:alpha_max) + amplitude_der * &
                                              (exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max))
            exp_coeff_der(1:alpha_max, k) = matmul( W, exp_coeff_der_temp(1:alpha_max) )
          end if
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_temp1(1:alpha_max) + pref_f * exp_coeff_temp2(1:alpha_max) )
        end if
      end do
    end do
!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   ***********************************

!   This is for debugging
    if( .false. )then
      open(10, file="coefficients.dat", status="unknown", access="append")
      write(10,*) exp_coeff(1:alpha_max, 1)
      close(10)
      if( do_derivatives )then
        open(10, file="derivatives.dat", status="unknown", access="append")
        write(10,*) exp_coeff_der(1:alpha_max, 1)
        close(10)
      end if
    end if

    deallocate( exp_coeff_temp1, exp_coeff_temp2, exp_coeff_der_temp )

  return
  end subroutine
!**************************************************************************










!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3(alpha_max, S, W)

    implicit none

    integer :: alpha_max, i, j, info
    real*8, intent(inout) :: W(:,:), S(:,:)
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:)
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

    do i = 1, alpha_max
      S(i,i) = 1.d0
      do j = i+1, alpha_max
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
    S = matmul(U,S)
    S = matmul(S,VT)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv )

  return
  end subroutine
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set augmented with a central Gaussian function. It is
! also needed to transform between original basis and orthonormal basis.
! It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3gauss(alpha_max, atom_sigma_in, rcut_hard_in, S, W)

    implicit none

    integer :: alpha_max, i, j, info, n
    real*8, intent(inout) :: W(:,:), S(:,:)
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:)
    real*8 :: s2, I_n, N_n, N_np1, I_np1, N_np2, I_np2, C2, sq2, pi, atom_sigma, rcut_hard
    real*8, intent(in) :: rcut_hard_in, atom_sigma_in
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

!   These are the overlap integrals for the polynomial functions
    do i = 1, alpha_max-1
      S(i,i) = 1.d0
      do j = i+1, alpha_max-1
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

!   These are the overlap integrals between the Gaussian and the polynomials
!   See derivation of radial expansion coefficients to understand this code
!   **** New basis ****
    atom_sigma = atom_sigma_in/rcut_hard_in
    rcut_hard = 1.d0
!   *******************
    s2 = atom_sigma**2
    sq2 = dsqrt(2.d0)
    pi = dacos(-1.d0)
    I_n = 0.d0
    N_n = 1.d0
    N_np1 = N_a(rcut_hard, -2)
    I_np1 = dsqrt(pi/2.d0) * atom_sigma * derf( rcut_hard/sq2/atom_sigma ) / N_np1
    C2 = s2 / rcut_hard
    do n = -1, alpha_max-1
      C2 = C2 * rcut_hard
      N_np2 = N_a(rcut_hard, n)
      I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
              + N_np1 * rcut_hard / N_np2 * I_np1 &
              - C2 / N_np2
      if(n > 0)then
!       Include the normalization factor of the Gaussian
        S(alpha_max, n) = I_np2 * sq2 / dsqrt(atom_sigma) / pi**0.25d0
        S(n, alpha_max) = S(alpha_max, n)
      end if
      N_n = N_np1
      N_np1 = N_np2
      I_n = I_np1
      I_np1 = I_np2
    end do
    S(alpha_max, alpha_max) = 1.d0

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
    S = matmul(U,S)
    S = matmul(S,VT)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv )

  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. Unlike the original function, it does not
! require blas/lapack to work and instead relies on pretabulated values.
!
  subroutine get_orthonormalization_matrix_poly3_tabulated(alpha_max, S, W)

    implicit none

    integer :: alpha_max
    real*8, intent(inout) :: W(:,:), S(:,:)

    select case(alpha_max)
      case(:0)
        write(*,*) "Bad value of alpha_max"
        stop
      case(1)
        W = reshape([1.0000000000000000], shape(W))
        S = reshape([1.0000000000000000], shape(S))
      case(2)
        W = reshape([5.9999999999999920, -5.2915026221291726, -5.2915026221291726, 5.9999999999999920 &
                    ], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.99215674164922152, 1.0000000000000000 &
                    ], shape(S))
      case(3)
        W = reshape([16.532871534059733, -29.078310649176263, 13.308493852760057, -29.078310649176263,  &
                    65.256761373659685, -36.000096455630107, 13.308493852760057, -36.000096455630107,  &
                    23.492063480194044], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.99215674164922152,  &
                    1.0000000000000000, 0.99498743710661997, 0.97499604304356913, 0.99498743710661997,  &
                    1.0000000000000000], shape(S))
      case(4)
        W = reshape([31.619473929163473, -82.899603631697133, 76.876069386827055, -24.858289195076921,  &
                    -82.899603631697133, 278.76362329353776, -309.81360154834925, 114.16667789081102,  &
                    76.876069386827055, -309.81360154834925, 401.66898510515369, -168.42692378140524,  &
                    -24.858289195076921, 114.16667789081102, -168.42692378140524, 79.877446522677999], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.99215674164922152, 1.0000000000000000, 0.99498743710661997, 0.98333216603563356,  &
                    0.97499604304356913, 0.99498743710661997, 1.0000000000000000, 0.99652172859178323,  &
                    0.95393920141694566, 0.98333216603563356, 0.99652172859178323, 1.0000000000000000 &
                    ], shape(S))
      case(5)
        W = reshape([48.701840798683406, -167.01751270577688, 232.29114353203897, -152.17733678177117,  &
                    38.937949480633065, -167.01751270577688, 745.79170301284648, -1251.5498958156577,  &
                    936.75356599200450, -263.84749322361466, 232.29114353203897, -1251.5498958156577,  &
                    2447.7952379389576, -2072.0455670547381, 643.96375585768624, -152.17733678177117,  &
                    936.75356599200450, -2072.0455670547381, 1959.7497439691933, -672.11823364760107,  &
                    38.937949480633065, -263.84749322361466, 643.96375585768624, -672.11823364760107,  &
                    253.84463194408610], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.99215674164922152, 1.0000000000000000, 0.99498743710661997,  &
                    0.98333216603563356, 0.96824583655185414, 0.97499604304356913, 0.99498743710661997,  &
                    1.0000000000000000, 0.99652172859178323, 0.98809481374347152, 0.95393920141694566,  &
                    0.98333216603563356, 0.99652172859178323, 1.0000000000000000, 0.99744571741206722,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000], shape(S))
      case(6)
        W = reshape([65.344980865540876, -270.96483658090318, 495.43665447034391, -487.28339480337212,  &
                    252.44181841585481, -54.246284684757562, -270.96483658090318, 1493.2157334161463,  &
                    -3338.9403306650820, 3783.7603589943619, -2167.7466416209732, 500.79512237094178,  &
                    495.43665447034391, -3338.9403306650820, 8760.0319024860873, -11249.630736045141,  &
                    7104.9320746203475, -1771.4461242822615, -487.28339480337212, 3783.7603589943619,  &
                    -11249.630736045141, 16097.101050349149, -11150.857982565711, 3007.1993994548748,  &
                    252.44181841585481, -2167.7466416209732, 7104.9320746203475, -11150.857982565711,  &
                    8419.3910413187914, -2457.9617037326980, -54.246284684757562, 500.79512237094178,  &
                    -1771.4461242822615, 3007.1993994548748, -2457.9617037326980, 776.42840231397884], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.99215674164922152, 1.0000000000000000,  &
                    0.99498743710661997, 0.98333216603563356, 0.96824583655185414, 0.95148591360407542,  &
                    0.97499604304356913, 0.99498743710661997, 1.0000000000000000, 0.99652172859178323,  &
                    0.98809481374347152, 0.97677102365552460, 0.95393920141694566, 0.98333216603563356,  &
                    0.99652172859178323, 1.0000000000000000, 0.99744571741206722, 0.99107124982123374,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.90905934288630952, 0.95148591360407542,  &
                    0.97677102365552460, 0.99107124982123374, 0.99804496391695696, 1.0000000000000000 &
                    ], shape(S))
      case(7)
        W = reshape([80.112087744414481, -380.83753196337295, 846.68851603934706, -1097.4660724104394,  &
                    854.00203553098982, -371.15373926820087, 69.379320671488060, -380.83753196337295,  &
                    2460.4826695826819, -6806.5629280860976, 10292.638409627154, -8930.9770566136904,  &
                    4194.6906070515997, -829.33608862620451, 846.68851603934706, -6806.5629280860976,  &
                    22347.744248873638, -38562.685975984692, 37023.176288616021, -18793.587621101400,  &
                    3945.6376780219634, -1097.4660724104394, 10292.638409627154, -38562.685975984692,  &
                    74347.213817597105, -78269.402026367752, 42877.575814051277, -9587.7268556762156,  &
                    854.00203553098982, -8930.9770566136904, 37023.176288616021, -78269.402026367752,  &
                    89507.134104252254, -52778.316794353603, 12594.783311461077, -371.15373926820087,  &
                    4194.6906070515997, -18793.587621101400, 42877.575814051277, -52778.316794353603,  &
                    33381.624954509396, -8510.6924230665845, 69.379320671488060, -829.33608862620451,  &
                    3945.6376780219634, -9587.7268556762156, 12594.783311461077, -8510.6924230665845,  &
                    2318.7297663139648], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.88712019959006128, 0.99215674164922152,  &
                    1.0000000000000000, 0.99498743710661997, 0.98333216603563356, 0.96824583655185414,  &
                    0.95148591360407542, 0.93404977361585861, 0.97499604304356913, 0.99498743710661997,  &
                    1.0000000000000000, 0.99652172859178323, 0.98809481374347152, 0.97677102365552460,  &
                    0.96378881965339736, 0.95393920141694566, 0.98333216603563356, 0.99652172859178323,  &
                    1.0000000000000000, 0.99744571741206722, 0.99107124982123374, 0.98226460284385697,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.99305547153730200, 0.90905934288630952,  &
                    0.95148591360407542, 0.97677102365552460, 0.99107124982123374, 0.99804496391695696,  &
                    1.0000000000000000, 0.99845559753396829, 0.88712019959006128, 0.93404977361585861,  &
                    0.96378881965339736, 0.98226460284385697, 0.99305547153730200, 0.99845559753396829,  &
                    1.0000000000000000], shape(S))
      case(8)
        W = reshape([92.539024757823469, -485.62295892374391, 1245.3098138485498, -1970.1315498890267,  &
                    2024.2063495743746, -1321.8047757365766, 499.34651114023347, -83.121761910171273,  &
                    -485.62295892374391, 3539.9434708348890, -11511.708266971060, 21569.895229194808,  &
                    -24986.945984915284, 17767.671463504408, -7134.3810386646874, 1241.2366889962939,  &
                    1245.3098138485498, -11511.708266971060, 45076.344590889283, -97380.089133774367,  &
                    125682.09286798064, -97005.646137706557, 41462.330953470708, -7568.2374319869477,  &
                    -1970.1315498890267, 21569.895229194808, -97380.089133774367, 236558.01704791249,  &
                    -335931.86898086715, 280178.31947277795, -127521.98045884802, 24497.997672243826,  &
                    2024.2063495743746, -24986.945984915284, 125682.09286798064, -335931.86898086715,  &
                    518562.62266775075, -464898.16962265247, 225191.00178569887, -45642.664557060503,  &
                    -1321.8047757365766, 17767.671463504408, -97005.646137706557, 280178.31947277795,  &
                    -464898.16962265247, 445491.35949089238, -229343.07661604194, 49131.675673804595,  &
                    499.34651114023347, -7134.3810386646874, 41462.330953470708, -127521.98045884802,  &
                    225191.00178569887, -229343.07661604194, 125237.47296520101, -28390.563739949452,  &
                    -83.121761910171273, 1241.2366889962939, -7568.2374319869477, 24497.997672243826,  &
                    -45642.664557060503, 49131.675673804595, -28390.563739949452, 6814.4475643740907], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.88712019959006128, 0.86602540378443860,  &
                    0.99215674164922152, 1.0000000000000000, 0.99498743710661997, 0.98333216603563356,  &
                    0.96824583655185414, 0.95148591360407542, 0.93404977361585861, 0.91651513899116799,  &
                    0.97499604304356913, 0.99498743710661997, 1.0000000000000000, 0.99652172859178323,  &
                    0.98809481374347152, 0.97677102365552460, 0.96378881965339736, 0.94991775959816649,  &
                    0.95393920141694566, 0.98333216603563356, 0.99652172859178323, 1.0000000000000000,  &
                    0.99744571741206722, 0.99107124982123374, 0.98226460284385697, 0.97192421422695907,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.99305547153730200, 0.98601329718326935,  &
                    0.90905934288630952, 0.95148591360407542, 0.97677102365552460, 0.99107124982123374,  &
                    0.99804496391695696, 1.0000000000000000, 0.99845559753396829, 0.99444440145743085,  &
                    0.88712019959006128, 0.93404977361585861, 0.96378881965339736, 0.98226460284385697,  &
                    0.99305547153730200, 0.99845559753396829, 1.0000000000000000, 0.99874921777190884,  &
                    0.86602540378443860, 0.91651513899116799, 0.94991775959816649, 0.97192421422695907,  &
                    0.98601329718326935, 0.99444440145743085, 0.99874921777190884, 1.0000000000000000 &
                    ], shape(S))
      case(9)
        W = reshape([102.73996023131654, -579.28008477226797, 1649.6329933089637, -3023.3784007460226,  &
                    3801.6152070879762, -3281.1404539213454, 1861.8484048267730, -625.49270608226868,  &
                    94.172577622820413, -579.28008477226797, 4626.3192102608891, -17024.980936140477,  &
                    37570.528484501476, -53968.531249935149, 51220.175635520885, -31105.502641416926,  &
                    10973.451092982890, -1712.0993861885026, 1649.6329933089637, -17024.980936140477,  &
                    76648.552042974829, -197491.45685094723, 319059.39758025052, -330888.18983548041,  &
                    214867.87633591256, -79759.988337705232, 12939.555536288681, -3023.3784007460226,  &
                    37570.528484501476, -197491.45685094723, 577420.51966726186, -1032806.0301229986,  &
                    1161968.8118061803, -805283.01678246027, 314926.85573368822, -53282.713938116933,  &
                    3801.6152070879762, -53968.531249935149, 319059.39758025052, -1032806.0301229986,  &
                    2015313.2813086673, -2441081.2560675377, 1800654.0587436187, -742277.93893866800,  &
                    131305.73650891884, -3281.1404539213454, 51220.175635520885, -330888.18983548041,  &
                    1161968.8118061803, -2441081.2560675377, 3159490.1849562805, -2472607.0610428583,  &
                    1074278.1958661408, -199099.57483824741, 1861.8484048267730, -31105.502641416926,  &
                    214867.87633591256, -805283.01678246027, 1800654.0587436187, -2472607.0610428583,  &
                    2045686.8262865648, -936138.89101480751, 182064.24577301860, -625.49270608226868,  &
                    10973.451092982890, -79759.988337705232, 314926.85573368822, -742277.93893866800,  &
                    1074278.1958661408, -936138.89101480751, 450714.89061671734, -92090.959139143888,  &
                    94.172577622820413, -1712.0993861885026, 12939.555536288681, -53282.713938116933,  &
                    131305.73650891884, -199099.57483824741, 182064.24577301860, -92090.959139143888,  &
                    19782.408284583180], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.88712019959006128, 0.86602540378443860,  &
                    0.84590516936330140, 0.99215674164922152, 1.0000000000000000, 0.99498743710661997,  &
                    0.98333216603563356, 0.96824583655185414, 0.95148591360407542, 0.93404977361585861,  &
                    0.91651513899116799, 0.89921841062113494, 0.97499604304356913, 0.99498743710661997,  &
                    1.0000000000000000, 0.99652172859178323, 0.98809481374347152, 0.97677102365552460,  &
                    0.96378881965339736, 0.94991775959816649, 0.93564551297569798, 0.95393920141694566,  &
                    0.98333216603563356, 0.99652172859178323, 1.0000000000000000, 0.99744571741206722,  &
                    0.99107124982123374, 0.98226460284385697, 0.97192421422695907, 0.96064535921058791,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.99305547153730200, 0.98601329718326935,  &
                    0.97758819057930046, 0.90905934288630952, 0.95148591360407542, 0.97677102365552460,  &
                    0.99107124982123374, 0.99804496391695696, 1.0000000000000000, 0.99845559753396829,  &
                    0.99444440145743085, 0.98868599666425949, 0.88712019959006128, 0.93404977361585861,  &
                    0.96378881965339736, 0.98226460284385697, 0.99305547153730200, 0.99845559753396829,  &
                    1.0000000000000000, 0.99874921777190884, 0.99545452192223205, 0.86602540378443860,  &
                    0.91651513899116799, 0.94991775959816649, 0.97192421422695907, 0.98601329718326935,  &
                    0.99444440145743085, 0.99874921777190884, 1.0000000000000000, 0.99896640799254144,  &
                    0.84590516936330140, 0.89921841062113494, 0.93564551297569798, 0.96064535921058791,  &
                    0.97758819057930046, 0.98868599666425949, 0.99545452192223205, 0.99896640799254144,  &
                    1.0000000000000000], shape(S))
      case(10)
        W = reshape([110.85899098970923, -656.32724208730247, 2002.0561806262745, -4030.0430894678625,  &
                    5747.1082633813776, -5876.5673849153736, 4220.3047498562373, -2019.2279545723384,  &
                    576.83178827809900, -74.279463179713289, -656.32724208730247, 5580.5224445491212,  &
                    -22331.612043043937, 54995.143701985420, -91095.538134775474, 104230.69173430053,  &
                    -81616.918643421988, 41841.558030521599, -12662.655255751506, 1715.2096007504856,  &
                    2002.0561806262745, -22331.612043043937, 110954.23971315096, -323698.56673909881,  &
                    611814.99361075275, -776207.90529295243, 659769.38292835664, -361528.21759789030,  &
                    115656.71780639346, -16430.693998954230, -4030.0430894678625, 54995.143701985420,  &
                    -323698.56673909881, 1086072.8812834970, -2301554.0722570983, 3206429.9814929799,  &
                    -2943668.1793391723, 1719782.1674985890, -580745.53209077963, 86416.331611407732,  &
                    5747.1082633813776, -91095.538134775474, 611814.99361075275, -2301554.0722570983,  &
                    5380578.4815178504, -8152091.3190899873, 8040206.2774941670, -4994941.7070384044,  &
                    1778403.4588345746, -277067.37753186666, -5876.5673849153736, 104230.69173430053,  &
                    -776207.90529295243, 3206429.9814929799, -8152091.3190899873, 13309104.674707124,  &
                    -14025079.427185630, 9238846.3749270160, -3464500.8524018270, 565144.53718240885,  &
                    4220.3047498562373, -81616.918643421988, 659769.38292835664, -2943668.1793391723,  &
                    8040206.2774941670, -14025079.427185630, 15706401.821684649, -10938310.914757373,  &
                    4315299.1100801602, -737221.24287861609, -2019.2279545723384, 41841.558030521599,  &
                    -361528.21759789030, 1719782.1674985890, -4994941.7070384044, 9238846.3749270160,  &
                    -10938310.914757373, 8029013.8383061187, -3328353.6668580547, 595670.14588087215,  &
                    576.83178827809900, -12662.655255751506, 115656.71780639346, -580745.53209077963,  &
                    1778403.4588345746, -3464500.8524018270, 4315299.1100801602, -3328353.6668580547,  &
                    1447833.8533988285, -271507.14094415621, -74.279463179713289, 1715.2096007504856,  &
                    -16430.693998954230, 86416.331611407732, -277067.37753186666, 565144.53718240885,  &
                    -737221.24287861609, 595670.14588087215, -271507.14094415621, 53355.279479472047], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.88712019959006128, 0.86602540378443860,  &
                    0.84590516936330140, 0.82679728470768454, 0.99215674164922152, 1.0000000000000000,  &
                    0.99498743710661997, 0.98333216603563356, 0.96824583655185414, 0.95148591360407542,  &
                    0.93404977361585861, 0.91651513899116799, 0.89921841062113494, 0.88235294117647056,  &
                    0.97499604304356913, 0.99498743710661997, 1.0000000000000000, 0.99652172859178323,  &
                    0.98809481374347152, 0.97677102365552460, 0.96378881965339736, 0.94991775959816649,  &
                    0.93564551297569798, 0.92128466398761111, 0.95393920141694566, 0.98333216603563356,  &
                    0.99652172859178323, 1.0000000000000000, 0.99744571741206722, 0.99107124982123374,  &
                    0.98226460284385697, 0.97192421422695907, 0.96064535921058791, 0.94882928301683922,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.99305547153730200, 0.98601329718326935,  &
                    0.97758819057930046, 0.96824583655185426, 0.90905934288630952, 0.95148591360407542,  &
                    0.97677102365552460, 0.99107124982123374, 0.99804496391695696, 1.0000000000000000,  &
                    0.99845559753396829, 0.99444440145743085, 0.98868599666425949, 0.98169181562325258,  &
                    0.88712019959006128, 0.93404977361585861, 0.96378881965339736, 0.98226460284385697,  &
                    0.99305547153730200, 0.99845559753396829, 1.0000000000000000, 0.99874921777190884,  &
                    0.99545452192223205, 0.99065885080469862, 0.86602540378443860, 0.91651513899116799,  &
                    0.94991775959816649, 0.97192421422695907, 0.98601329718326935, 0.99444440145743085,  &
                    0.99874921777190884, 1.0000000000000000, 0.99896640799254144, 0.99621210759909562,  &
                    0.84590516936330140, 0.89921841062113494, 0.93564551297569798, 0.96064535921058791,  &
                    0.97758819057930046, 0.98868599666425949, 0.99545452192223205, 0.99896640799254144,  &
                    1.0000000000000000, 0.99913156735681652, 0.82679728470768454, 0.88235294117647056,  &
                    0.92128466398761111, 0.94882928301683922, 0.96824583655185426, 0.98169181562325258,  &
                    0.99065885080469862, 0.99621210759909562, 0.99913156735681652, 1.0000000000000000 &
                    ], shape(S))
      case(11)
        W = reshape([117.21235214772904, -714.35806297792158, 2245.1643321280740, -4610.3933792100061,  &
                    6511.4009202362895, -6150.8302224825920, 3391.5937860140148, -436.01742880721537,  &
                    -738.92275647372799, 484.31715587976799, -98.454105559226775, -714.35806297792158,  &
                    6269.8952732304706, -25875.069939432633, 65110.259187545853, -107732.04536419440,  &
                    117390.86907425059, -78379.654709546943, 23457.309598262604, 5607.1440832265725,  &
                    -6684.8214165052414, 1550.5418761193073, 2245.1643321280740, -25875.069939432633,  &
                    132599.60954916809, -396361.31134008802, 756826.91023577540, -944113.58367922041,  &
                    746475.35184183146, -329819.32141719793, 40631.572575033351, 26062.477344374791,  &
                    -8671.4089555561513, -4610.3933792100061, 65110.259187545853, -396361.31134008802,  &
                    1372463.1491815776, -2988502.0270795608, 4244631.4027743703, -3916425.4633714431,  &
                    2232511.7325943527, -677804.31608534034, 55391.291470506068, 13595.784351342649,  &
                    6511.4009202362895, -107732.04536419440, 756826.91023577540, -2988502.0270795608,  &
                    7374201.5953415921, -11897330.128776217, 12685063.180279443, -8765094.4377386160,  &
                    3697769.2054832382, -829684.91928406083, 67971.552751104493, -6150.8302224825920,  &
                    117390.86907425059, -944113.58367922041, 4244631.4027743703, -11897330.128776217,  &
                    21851270.467823889, -26754279.226484302, 21633158.682814363, -11077651.903835505,  &
                    3246049.5859478437, -412975.14465735806, 3391.5937860140148, -78379.654709546943,  &
                    746475.35184183146, -3916425.4633714431, 12685063.180279443, -26754279.226484302,  &
                    37494848.639549591, -34689165.786953673, 20377863.358143974, -6891923.6598007092,  &
                    1022531.8799494690, -436.01742880721537, 23457.309598262604, -329819.32141719793,  &
                    2232511.7325943527, -8765094.4377386160, 21633158.682814363, -34689165.786953673,  &
                    36151306.746357322, -23647318.572119914, 8826384.4847903624, -1434984.5992588799,  &
                    -738.92275647372799, 5607.1440832265725, 40631.572575033351, -677804.31608534034,  &
                    3697769.2054832382, -11077651.903835505, 20377863.358143974, -23647318.572119914,  &
                    16913806.476479985, -6818990.6947852084, 1186826.9821223991, 484.31715587976799,  &
                    -6684.8214165052414, 26062.477344374791, 55391.291470506068, -829684.91928406083,  &
                    3246049.5859478437, -6891923.6598007092, 8826384.4847903624, -6818990.6947852084,  &
                    2933523.2297482090, -540611.16396153928, -98.454105559226775, 1550.5418761193073,  &
                    -8671.4089555561513, 13595.784351342649, 67971.552751104493, -412975.14465735806,  &
                    1022531.8799494690, -1434984.5992588799, 1186826.9821223991, -540611.16396153928,  &
                    104864.79597292397], shape(W))
        S = reshape([1.0000000000000000, 0.99215674164922152, 0.97499604304356913, 0.95393920141694566,  &
                    0.93154097872359987, 0.90905934288630952, 0.88712019959006128, 0.86602540378443860,  &
                    0.84590516936330140, 0.82679728470768454, 0.80868982852161886, 0.99215674164922152,  &
                    1.0000000000000000, 0.99498743710661997, 0.98333216603563356, 0.96824583655185414,  &
                    0.95148591360407542, 0.93404977361585861, 0.91651513899116799, 0.89921841062113494,  &
                    0.88235294117647056, 0.86602540378443871, 0.97499604304356913, 0.99498743710661997,  &
                    1.0000000000000000, 0.99652172859178323, 0.98809481374347152, 0.97677102365552460,  &
                    0.96378881965339736, 0.94991775959816649, 0.93564551297569798, 0.92128466398761111,  &
                    0.90703620734810975, 0.95393920141694566, 0.98333216603563356, 0.99652172859178323,  &
                    1.0000000000000000, 0.99744571741206722, 0.99107124982123374, 0.98226460284385697,  &
                    0.97192421422695907, 0.96064535921058791, 0.94882928301683922, 0.93674969975975964,  &
                    0.93154097872359987, 0.96824583655185414, 0.98809481374347152, 0.99744571741206722,  &
                    1.0000000000000000, 0.99804496391695696, 0.99305547153730200, 0.98601329718326935,  &
                    0.97758819057930046, 0.96824583655185426, 0.95831484749990992, 0.90905934288630952,  &
                    0.95148591360407542, 0.97677102365552460, 0.99107124982123374, 0.99804496391695696,  &
                    1.0000000000000000, 0.99845559753396829, 0.99444440145743085, 0.98868599666425949,  &
                    0.98169181562325258, 0.97383114934675230, 0.88712019959006128, 0.93404977361585861,  &
                    0.96378881965339736, 0.98226460284385697, 0.99305547153730200, 0.99845559753396829,  &
                    1.0000000000000000, 0.99874921777190884, 0.99545452192223205, 0.99065885080469862,  &
                    0.98476101329618471, 0.86602540378443860, 0.91651513899116799, 0.94991775959816649,  &
                    0.97192421422695907, 0.98601329718326935, 0.99444440145743085, 0.99874921777190884,  &
                    1.0000000000000000, 0.99896640799254144, 0.99621210759909562, 0.99215674164922152,  &
                    0.84590516936330140, 0.89921841062113494, 0.93564551297569798, 0.96064535921058791,  &
                    0.97758819057930046, 0.98868599666425949, 0.99545452192223205, 0.99896640799254144,  &
                    1.0000000000000000, 0.99913156735681652, 0.99679486355016889, 0.82679728470768454,  &
                    0.88235294117647056, 0.92128466398761111, 0.94882928301683922, 0.96824583655185426,  &
                    0.98169181562325258, 0.99065885080469862, 0.99621210759909562, 0.99913156735681652,  &
                    1.0000000000000000, 0.99926008128973698, 0.80868982852161886, 0.86602540378443871,  &
                    0.90703620734810975, 0.93674969975975964, 0.95831484749990992, 0.97383114934675230,  &
                    0.98476101329618471, 0.99215674164922152, 0.99679486355016889, 0.99926008128973698,  &
                    1.0000000000000000], shape(S))
      case(12:)
        write(*,*) "Bad value of alpha_max"
        stop
    end select

  return
  end subroutine
!**************************************************************************


end module soap_turbo_radial
