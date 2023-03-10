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
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:), Sc(:,:)
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( Sc(1:alpha_max, 1:alpha_max) )
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
!    S = matmul(U,S)
!    S = matmul(S,VT)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, U, alpha_max, S, alpha_max, 0.d0, Sc, alpha_max)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, Sc, alpha_max, VT, alpha_max, 0.d0, S, alpha_max)
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

    deallocate( Sb, U, Vt, svd, work, ipiv, Sc )

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
    real*8, allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:), Sc(:,:)
    real*8 :: s2, I_n, N_n, N_np1, I_np1, N_np2, I_np2, C2, sq2, pi, atom_sigma, rcut_hard
    real*8, intent(in) :: rcut_hard_in, atom_sigma_in
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( Sc(1:alpha_max, 1:alpha_max) )
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
!    S = matmul(U,S)
!    S = matmul(S,VT)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, U, alpha_max, S, alpha_max, 0.d0, Sc, alpha_max)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, Sc, alpha_max, VT, alpha_max, 0.d0, S, alpha_max)
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

    deallocate( Sb, U, Vt, svd, work, ipiv, Sc )

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
        W = reshape([1.0000000000000000d0], shape(W))
        S = reshape([1.0000000000000000d0], shape(S))
      case(2)
        W = reshape([5.9999999999999920d0, -5.2915026221291726d0, -5.2915026221291726d0,  &
                    5.9999999999999920d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0], shape(S))
      case(3)
        W = reshape([16.532871534059733d0, -29.078310649176263d0, 13.308493852760057d0, -29.078310649176263d0,  &
                    65.256761373659685d0, -36.000096455630107d0, 13.308493852760057d0, -36.000096455630107d0,  &
                    23.492063480194044d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0], shape(S))
      case(4)
        W = reshape([31.619473929163473d0, -82.899603631697133d0, 76.876069386827055d0, -24.858289195076921d0,  &
                    -82.899603631697133d0, 278.76362329353776d0, -309.81360154834925d0, 114.16667789081102d0,  &
                    76.876069386827055d0, -309.81360154834925d0, 401.66898510515369d0, -168.42692378140524d0,  &
                    -24.858289195076921d0, 114.16667789081102d0, -168.42692378140524d0, 79.877446522677999d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0], shape(S))
      case(5)
        W = reshape([48.701840798683406d0, -167.01751270577688d0, 232.29114353203897d0, -152.17733678177117d0,  &
                    38.937949480633065d0, -167.01751270577688d0, 745.79170301284648d0, -1251.5498958156577d0,  &
                    936.75356599200450d0, -263.84749322361466d0, 232.29114353203897d0, -1251.5498958156577d0,  &
                    2447.7952379389576d0, -2072.0455670547381d0, 643.96375585768624d0, -152.17733678177117d0,  &
                    936.75356599200450d0, -2072.0455670547381d0, 1959.7497439691933d0, -672.11823364760107d0,  &
                    38.937949480633065d0, -263.84749322361466d0, 643.96375585768624d0, -672.11823364760107d0,  &
                    253.84463194408610d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0], shape(S))
      case(6)
        W = reshape([65.344980865540876d0, -270.96483658090318d0, 495.43665447034391d0, -487.28339480337212d0,  &
                    252.44181841585481d0, -54.246284684757562d0, -270.96483658090318d0, 1493.2157334161463d0,  &
                    -3338.9403306650820d0, 3783.7603589943619d0, -2167.7466416209732d0, 500.79512237094178d0,  &
                    495.43665447034391d0, -3338.9403306650820d0, 8760.0319024860873d0, -11249.630736045141d0,  &
                    7104.9320746203475d0, -1771.4461242822615d0, -487.28339480337212d0, 3783.7603589943619d0,  &
                    -11249.630736045141d0, 16097.101050349149d0, -11150.857982565711d0, 3007.1993994548748d0,  &
                    252.44181841585481d0, -2167.7466416209732d0, 7104.9320746203475d0, -11150.857982565711d0,  &
                    8419.3910413187914d0, -2457.9617037326980d0, -54.246284684757562d0, 500.79512237094178d0,  &
                    -1771.4461242822615d0, 3007.1993994548748d0, -2457.9617037326980d0, 776.42840231397884d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0], shape(S))
      case(7)
        W = reshape([80.112087744414481d0, -380.83753196337295d0, 846.68851603934706d0, -1097.4660724104394d0,  &
                    854.00203553098982d0, -371.15373926820087d0, 69.379320671488060d0, -380.83753196337295d0,  &
                    2460.4826695826819d0, -6806.5629280860976d0, 10292.638409627154d0, -8930.9770566136904d0,  &
                    4194.6906070515997d0, -829.33608862620451d0, 846.68851603934706d0, -6806.5629280860976d0,  &
                    22347.744248873638d0, -38562.685975984692d0, 37023.176288616021d0, -18793.587621101400d0,  &
                    3945.6376780219634d0, -1097.4660724104394d0, 10292.638409627154d0, -38562.685975984692d0,  &
                    74347.213817597105d0, -78269.402026367752d0, 42877.575814051277d0, -9587.7268556762156d0,  &
                    854.00203553098982d0, -8930.9770566136904d0, 37023.176288616021d0, -78269.402026367752d0,  &
                    89507.134104252254d0, -52778.316794353603d0, 12594.783311461077d0, -371.15373926820087d0,  &
                    4194.6906070515997d0, -18793.587621101400d0, 42877.575814051277d0, -52778.316794353603d0,  &
                    33381.624954509396d0, -8510.6924230665845d0, 69.379320671488060d0, -829.33608862620451d0,  &
                    3945.6376780219634d0, -9587.7268556762156d0, 12594.783311461077d0, -8510.6924230665845d0,  &
                    2318.7297663139648d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0], shape(S))
      case(8)
        W = reshape([92.539024757823469d0, -485.62295892374391d0, 1245.3098138485498d0, -1970.1315498890267d0,  &
                    2024.2063495743746d0, -1321.8047757365766d0, 499.34651114023347d0, -83.121761910171273d0,  &
                    -485.62295892374391d0, 3539.9434708348890d0, -11511.708266971060d0, 21569.895229194808d0,  &
                    -24986.945984915284d0, 17767.671463504408d0, -7134.3810386646874d0, 1241.2366889962939d0,  &
                    1245.3098138485498d0, -11511.708266971060d0, 45076.344590889283d0, -97380.089133774367d0,  &
                    125682.09286798064d0, -97005.646137706557d0, 41462.330953470708d0, -7568.2374319869477d0,  &
                    -1970.1315498890267d0, 21569.895229194808d0, -97380.089133774367d0, 236558.01704791249d0,  &
                    -335931.86898086715d0, 280178.31947277795d0, -127521.98045884802d0, 24497.997672243826d0,  &
                    2024.2063495743746d0, -24986.945984915284d0, 125682.09286798064d0, -335931.86898086715d0,  &
                    518562.62266775075d0, -464898.16962265247d0, 225191.00178569887d0, -45642.664557060503d0,  &
                    -1321.8047757365766d0, 17767.671463504408d0, -97005.646137706557d0, 280178.31947277795d0,  &
                    -464898.16962265247d0, 445491.35949089238d0, -229343.07661604194d0, 49131.675673804595d0,  &
                    499.34651114023347d0, -7134.3810386646874d0, 41462.330953470708d0, -127521.98045884802d0,  &
                    225191.00178569887d0, -229343.07661604194d0, 125237.47296520101d0, -28390.563739949452d0,  &
                    -83.121761910171273d0, 1241.2366889962939d0, -7568.2374319869477d0, 24497.997672243826d0,  &
                    -45642.664557060503d0, 49131.675673804595d0, -28390.563739949452d0, 6814.4475643740907d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0,  &
                    0.91651513899116799d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0,  &
                    0.97192421422695907d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0,  &
                    0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0,  &
                    0.99444440145743085d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0,  &
                    0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0,  &
                    1.0000000000000000d0], shape(S))
      case(9)
        W = reshape([102.73996023131654d0, -579.28008477226797d0, 1649.6329933089637d0, -3023.3784007460226d0,  &
                    3801.6152070879762d0, -3281.1404539213454d0, 1861.8484048267730d0, -625.49270608226868d0,  &
                    94.172577622820413d0, -579.28008477226797d0, 4626.3192102608891d0, -17024.980936140477d0,  &
                    37570.528484501476d0, -53968.531249935149d0, 51220.175635520885d0, -31105.502641416926d0,  &
                    10973.451092982890d0, -1712.0993861885026d0, 1649.6329933089637d0, -17024.980936140477d0,  &
                    76648.552042974829d0, -197491.45685094723d0, 319059.39758025052d0, -330888.18983548041d0,  &
                    214867.87633591256d0, -79759.988337705232d0, 12939.555536288681d0, -3023.3784007460226d0,  &
                    37570.528484501476d0, -197491.45685094723d0, 577420.51966726186d0, -1032806.0301229986d0,  &
                    1161968.8118061803d0, -805283.01678246027d0, 314926.85573368822d0, -53282.713938116933d0,  &
                    3801.6152070879762d0, -53968.531249935149d0, 319059.39758025052d0, -1032806.0301229986d0,  &
                    2015313.2813086673d0, -2441081.2560675377d0, 1800654.0587436187d0, -742277.93893866800d0,  &
                    131305.73650891884d0, -3281.1404539213454d0, 51220.175635520885d0, -330888.18983548041d0,  &
                    1161968.8118061803d0, -2441081.2560675377d0, 3159490.1849562805d0, -2472607.0610428583d0,  &
                    1074278.1958661408d0, -199099.57483824741d0, 1861.8484048267730d0, -31105.502641416926d0,  &
                    214867.87633591256d0, -805283.01678246027d0, 1800654.0587436187d0, -2472607.0610428583d0,  &
                    2045686.8262865648d0, -936138.89101480751d0, 182064.24577301860d0, -625.49270608226868d0,  &
                    10973.451092982890d0, -79759.988337705232d0, 314926.85573368822d0, -742277.93893866800d0,  &
                    1074278.1958661408d0, -936138.89101480751d0, 450714.89061671734d0, -92090.959139143888d0,  &
                    94.172577622820413d0, -1712.0993861885026d0, 12939.555536288681d0, -53282.713938116933d0,  &
                    131305.73650891884d0, -199099.57483824741d0, 182064.24577301860d0, -92090.959139143888d0,  &
                    19782.408284583180d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0,  &
                    0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0,  &
                    0.96064535921058791d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.90905934288630952d0, 0.95148591360407542d0,  &
                    0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0,  &
                    0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0,  &
                    0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0,  &
                    0.99896640799254144d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0], shape(S))
      case(10)
        W = reshape([110.85899098970923d0, -656.32724208730247d0, 2002.0561806262745d0, -4030.0430894678625d0,  &
                    5747.1082633813776d0, -5876.5673849153736d0, 4220.3047498562373d0, -2019.2279545723384d0,  &
                    576.83178827809900d0, -74.279463179713289d0, -656.32724208730247d0, 5580.5224445491212d0,  &
                    -22331.612043043937d0, 54995.143701985420d0, -91095.538134775474d0, 104230.69173430053d0,  &
                    -81616.918643421988d0, 41841.558030521599d0, -12662.655255751506d0, 1715.2096007504856d0,  &
                    2002.0561806262745d0, -22331.612043043937d0, 110954.23971315096d0, -323698.56673909881d0,  &
                    611814.99361075275d0, -776207.90529295243d0, 659769.38292835664d0, -361528.21759789030d0,  &
                    115656.71780639346d0, -16430.693998954230d0, -4030.0430894678625d0, 54995.143701985420d0,  &
                    -323698.56673909881d0, 1086072.8812834970d0, -2301554.0722570983d0, 3206429.9814929799d0,  &
                    -2943668.1793391723d0, 1719782.1674985890d0, -580745.53209077963d0, 86416.331611407732d0,  &
                    5747.1082633813776d0, -91095.538134775474d0, 611814.99361075275d0, -2301554.0722570983d0,  &
                    5380578.4815178504d0, -8152091.3190899873d0, 8040206.2774941670d0, -4994941.7070384044d0,  &
                    1778403.4588345746d0, -277067.37753186666d0, -5876.5673849153736d0, 104230.69173430053d0,  &
                    -776207.90529295243d0, 3206429.9814929799d0, -8152091.3190899873d0, 13309104.674707124d0,  &
                    -14025079.427185630d0, 9238846.3749270160d0, -3464500.8524018270d0, 565144.53718240885d0,  &
                    4220.3047498562373d0, -81616.918643421988d0, 659769.38292835664d0, -2943668.1793391723d0,  &
                    8040206.2774941670d0, -14025079.427185630d0, 15706401.821684649d0, -10938310.914757373d0,  &
                    4315299.1100801602d0, -737221.24287861609d0, -2019.2279545723384d0, 41841.558030521599d0,  &
                    -361528.21759789030d0, 1719782.1674985890d0, -4994941.7070384044d0, 9238846.3749270160d0,  &
                    -10938310.914757373d0, 8029013.8383061187d0, -3328353.6668580547d0, 595670.14588087215d0,  &
                    576.83178827809900d0, -12662.655255751506d0, 115656.71780639346d0, -580745.53209077963d0,  &
                    1778403.4588345746d0, -3464500.8524018270d0, 4315299.1100801602d0, -3328353.6668580547d0,  &
                    1447833.8533988285d0, -271507.14094415621d0, -74.279463179713289d0, 1715.2096007504856d0,  &
                    -16430.693998954230d0, 86416.331611407732d0, -277067.37753186666d0, 565144.53718240885d0,  &
                    -737221.24287861609d0, 595670.14588087215d0, -271507.14094415621d0, 53355.279479472047d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0,  &
                    0.88235294117647056d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.93564551297569798d0, 0.92128466398761111d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0,  &
                    0.94882928301683922d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0,  &
                    0.98169181562325258d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.99545452192223205d0, 0.99065885080469862d0, 0.86602540378443860d0,  &
                    0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0,  &
                    0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0,  &
                    0.99621210759909562d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.82679728470768454d0,  &
                    0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0, 0.96824583655185426d0,  &
                    0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0, 0.99913156735681652d0,  &
                    1.0000000000000000d0], shape(S))
      case(11)
        W = reshape([117.21235214772904d0, -714.35806297792158d0, 2245.1643321280740d0, -4610.3933792100061d0,  &
                    6511.4009202362895d0, -6150.8302224825920d0, 3391.5937860140148d0, -436.01742880721537d0,  &
                    -738.92275647372799d0, 484.31715587976799d0, -98.454105559226775d0, -714.35806297792158d0,  &
                    6269.8952732304706d0, -25875.069939432633d0, 65110.259187545853d0, -107732.04536419440d0,  &
                    117390.86907425059d0, -78379.654709546943d0, 23457.309598262604d0, 5607.1440832265725d0,  &
                    -6684.8214165052414d0, 1550.5418761193073d0, 2245.1643321280740d0, -25875.069939432633d0,  &
                    132599.60954916809d0, -396361.31134008802d0, 756826.91023577540d0, -944113.58367922041d0,  &
                    746475.35184183146d0, -329819.32141719793d0, 40631.572575033351d0, 26062.477344374791d0,  &
                    -8671.4089555561513d0, -4610.3933792100061d0, 65110.259187545853d0, -396361.31134008802d0,  &
                    1372463.1491815776d0, -2988502.0270795608d0, 4244631.4027743703d0, -3916425.4633714431d0,  &
                    2232511.7325943527d0, -677804.31608534034d0, 55391.291470506068d0, 13595.784351342649d0,  &
                    6511.4009202362895d0, -107732.04536419440d0, 756826.91023577540d0, -2988502.0270795608d0,  &
                    7374201.5953415921d0, -11897330.128776217d0, 12685063.180279443d0, -8765094.4377386160d0,  &
                    3697769.2054832382d0, -829684.91928406083d0, 67971.552751104493d0, -6150.8302224825920d0,  &
                    117390.86907425059d0, -944113.58367922041d0, 4244631.4027743703d0, -11897330.128776217d0,  &
                    21851270.467823889d0, -26754279.226484302d0, 21633158.682814363d0, -11077651.903835505d0,  &
                    3246049.5859478437d0, -412975.14465735806d0, 3391.5937860140148d0, -78379.654709546943d0,  &
                    746475.35184183146d0, -3916425.4633714431d0, 12685063.180279443d0, -26754279.226484302d0,  &
                    37494848.639549591d0, -34689165.786953673d0, 20377863.358143974d0, -6891923.6598007092d0,  &
                    1022531.8799494690d0, -436.01742880721537d0, 23457.309598262604d0, -329819.32141719793d0,  &
                    2232511.7325943527d0, -8765094.4377386160d0, 21633158.682814363d0, -34689165.786953673d0,  &
                    36151306.746357322d0, -23647318.572119914d0, 8826384.4847903624d0, -1434984.5992588799d0,  &
                    -738.92275647372799d0, 5607.1440832265725d0, 40631.572575033351d0, -677804.31608534034d0,  &
                    3697769.2054832382d0, -11077651.903835505d0, 20377863.358143974d0, -23647318.572119914d0,  &
                    16913806.476479985d0, -6818990.6947852084d0, 1186826.9821223991d0, 484.31715587976799d0,  &
                    -6684.8214165052414d0, 26062.477344374791d0, 55391.291470506068d0, -829684.91928406083d0,  &
                    3246049.5859478437d0, -6891923.6598007092d0, 8826384.4847903624d0, -6818990.6947852084d0,  &
                    2933523.2297482090d0, -540611.16396153928d0, -98.454105559226775d0, 1550.5418761193073d0,  &
                    -8671.4089555561513d0, 13595.784351342649d0, 67971.552751104493d0, -412975.14465735806d0,  &
                    1022531.8799494690d0, -1434984.5992588799d0, 1186826.9821223991d0, -540611.16396153928d0,  &
                    104864.79597292397d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.80868982852161886d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0,  &
                    0.89921841062113494d0, 0.88235294117647056d0, 0.86602540378443871d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.92128466398761111d0, 0.90703620734810975d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0, 0.94882928301683922d0,  &
                    0.93674969975975964d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.95831484749990992d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0,  &
                    0.98868599666425949d0, 0.98169181562325258d0, 0.97383114934675230d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.99065885080469862d0, 0.98476101329618471d0, 0.86602540378443860d0, 0.91651513899116799d0,  &
                    0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0,  &
                    0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0, 0.99621210759909562d0,  &
                    0.99215674164922152d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.99679486355016889d0,  &
                    0.82679728470768454d0, 0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0,  &
                    0.96824583655185426d0, 0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0,  &
                    0.99913156735681652d0, 1.0000000000000000d0, 0.99926008128973698d0, 0.80868982852161886d0,  &
                    0.86602540378443871d0, 0.90703620734810975d0, 0.93674969975975964d0, 0.95831484749990992d0,  &
                    0.97383114934675230d0, 0.98476101329618471d0, 0.99215674164922152d0, 0.99679486355016889d0,  &
                    0.99926008128973698d0, 1.0000000000000000d0], shape(S))
      case(12:)
        write(*,*) "Bad value of alpha_max"
        stop
    end select

  return
  end subroutine
!**************************************************************************


end module soap_turbo_radial
