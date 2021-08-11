! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2021, Miguel A. Caro
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

module soap_turbo_desc

  use soap_turbo_radial
  use soap_turbo_angular

  contains

!**************************************************************************
  subroutine get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                      thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                      atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, &
                      amplitude_scaling, radial_enhancement, central_weight, basis, scaling_mode, do_timing, &
                      do_derivatives, compress_soap, compress_soap_indices, soap, soap_cart_der)

  implicit none

!-------------------
! Input variables
  real*8, intent(in) :: rjs(:), thetas(:), phis(:)
  real*8, intent(in) :: amplitude_scaling(:), atom_sigma_r_scaling(:), atom_sigma_t(:), atom_sigma_t_scaling(:)
  real*8, intent(in) :: central_weight(:), atom_sigma_r(:), global_scaling(:)
  real*8, intent(in) :: nf(:), rcut_hard(:), rcut_soft(:)

  integer, intent(in) :: n_species, radial_enhancement, species(:,:), species_multiplicity(:)
  integer, intent(in) :: n_sites, n_neigh(:), l_max, n_atom_pairs, alpha_max(:), compress_soap_indices(:)

  logical, intent(in) :: do_derivatives, do_timing, mask(:,:), compress_soap

  character(*), intent(in) :: basis, scaling_mode

! Output variables
  real*8, intent(inout) :: soap(:,:), soap_cart_der(:,:,:)
!-------------------


!-------------------
! Internal variables
  complex*16, allocatable :: angular_exp_coeff(:,:), cnk(:,:,:)
  complex*16, allocatable :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:), cnk_rad_der(:,:,:)
  complex*16, allocatable :: cnk_azi_der(:,:,:), angular_exp_coeff_pol_der(:,:), cnk_pol_der(:,:,:)
  complex*16, allocatable :: eimphi(:), prefm(:), eimphi_rad_der(:)

  real*8, allocatable, save :: W(:,:), S(:,:), multiplicity_array(:)
  real*8, allocatable :: soap_rad_der(:,:), sqrt_dot_p(:), soap_azi_der(:,:)
  real*8, allocatable :: W_temp(:,:), S_temp(:,:)
  real*8, allocatable :: radial_exp_coeff(:,:), soap_pol_der(:,:)
  real*8, allocatable :: preflm(:), plm_array(:), prefl(:), fact_array(:), prefl_rad_der(:)
  real*8, allocatable :: radial_exp_coeff_der(:,:), cnk_slice(:)
  real*8 :: amplitude, multiplicity, pi, rcut_max
  real*8 :: radial_time, angular_time, coeff_time, time3, total_time, soap_time, time1, time2, compress_time, &
            memory_time, basis_time

  integer, allocatable :: i_beg(:), i_end(:), starting_index(:)
  integer, save :: n_max_prev
  integer :: k_max, n_max
  integer :: i, counter, j, k, n_soap, k2, k3, n, l, m, np, counter2
  logical, allocatable :: do_central(:), skip_soap_component(:,:,:)
  logical, save :: recompute_basis = .true.
!-------------------

  if( do_timing )then
    call cpu_time(time3)
  end if

!-------------------
! Constants
  pi = dacos(-1.d0)
!-------------------



  if( do_timing )then
    call cpu_time(time1)
  end if

! Check if we need to expand the central atom
  allocate( do_central(1:n_species) )
  do i = 1, n_species
    if( central_weight(i) /= 0.d0 )then
      do_central(i) = .true.
    else
      do_central(i) = .false.
    end if
  end do


  rcut_max = 0.d0
  do i = 1, n_species
    if( rcut_hard(i) > rcut_max )then
      rcut_max = rcut_hard(i)
    end if
  end do


! This assigns begin and end indices to the components of the radial basis
  allocate( i_beg(1:n_species) )
  allocate( i_end(1:n_species) )
  i_beg(1) = 1
  i_end(1) = alpha_max(1)
  do i = 2, n_species
    i_beg(i) = i_end(i-1) + 1
    i_end(i) = i_beg(i) + alpha_max(i) - 1
  end do


! Do some array allocation for the SOAP expansion part
  k_max = 1 + l_max*(l_max+1)/2 + l_max
  allocate( preflm(1:k_max) )
  allocate( plm_array(1:k_max) )
  allocate( eimphi(1:k_max) )
  allocate( prefl(0:l_max) )
  allocate( prefm(0:l_max) )
  allocate( fact_array(1:l_max) )
  call get_preflm(preflm, l_max)
  if( do_derivatives )then
    allocate( prefl_rad_der(0:l_max) )
    allocate( eimphi_rad_der(1:k_max) )
  end if


  if( do_timing )then
    call cpu_time(time2)
    memory_time = time2 - time1
    time1 = time2
  end if


! This is to build the radial basis
  n_max = 0
  do i = 1, n_species
    n_max = n_max + alpha_max(i)
  end do
  if( n_max_prev /= n_max )then
    n_max_prev = n_max
    recompute_basis = .true.
  end if
  if( recompute_basis )then
    if( allocated(W) .or. allocated(S) )then
      deallocate(W, S)
    end if
    allocate( W(1:n_max, 1:n_max) )
    allocate( S(1:n_max, 1:n_max) )
    W = 0.d0
    S = 0.d0
!   This is done per species. Each loop iteration modifies the slice of the W and S matrices that
!   corresponds to the species in question. The radial basis functions for species A are always
!   assumed orthogonal to the basis functions for species B. W and S are therefore block diagonal.
    do i = 1, n_species
!     We pass these temp arrays with the right size because the Lapack/Blas routines internally fail
!     if the memory is not contiguous
      allocate( S_temp(1:alpha_max(i), 1:alpha_max(i)) )
      allocate( W_temp(1:alpha_max(i), 1:alpha_max(i)) )
      S_temp = 0.d0
      W_temp = 0.d0
      if( basis == "poly3gauss" )then
        call get_orthonormalization_matrix_poly3gauss(alpha_max(i), atom_sigma_r(i), rcut_hard(i), S_temp, W_temp)
      else if( basis == "poly3" )then
        call get_orthonormalization_matrix_poly3(alpha_max(i), S_temp, W_temp)
      end if
      S(i_beg(i):i_end(i), i_beg(i):i_end(i)) = S_temp
      W(i_beg(i):i_end(i), i_beg(i):i_end(i)) = W_temp
      deallocate( S_temp, W_temp )
    end do
  end if

  if( do_timing )then
    call cpu_time(time2)
    basis_time = time2 - time1
  end if


! This is for the expansion coefficients and the soap vectors
  allocate( radial_exp_coeff(1:n_max, 1:n_atom_pairs) )
  allocate( angular_exp_coeff(1:k_max, 1:n_atom_pairs) )
  angular_exp_coeff = 0.d0
  allocate( cnk( 1:k_max, 1:n_max, 1:n_sites) )
  cnk = 0.d0
! Handle SOAP compression here
  allocate( skip_soap_component(0:l_max, 1:n_max, 1:n_max) )
  skip_soap_component = .false.
  if( compress_soap )then
    if( do_timing )then
      call cpu_time(time1)
    end if
    n_soap = size(compress_soap_indices)
    skip_soap_component = .true.
    counter = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          counter = counter + 1
          do i = 1, n_soap
            if( compress_soap_indices(i) == counter )then
              skip_soap_component(l, np, n) = .false.
              exit
            end if
          end do
        end do
      end do
    end do
    if( do_timing )then
      call cpu_time(time2)
      compress_time = time2 - time1
    end if
  else
    n_soap = n_max*(n_max+1)/2 * (l_max+1)
  end if

  if( do_timing )then
    call cpu_time(time1)
  end if

  allocate( sqrt_dot_p(1:n_sites) )
  sqrt_dot_p = 0.d0
  if( do_derivatives )then
    allocate( radial_exp_coeff_der(1:n_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_rad_der(1:k_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_azi_der(1:k_max, 1:n_atom_pairs) )
    allocate( angular_exp_coeff_pol_der(1:k_max, 1:n_atom_pairs) )
    allocate( cnk_rad_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( cnk_azi_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( cnk_pol_der( 1:k_max, 1:n_max, 1:n_atom_pairs) )
    allocate( soap_rad_der(1:n_soap, 1:n_atom_pairs) )
    allocate( soap_azi_der(1:n_soap, 1:n_atom_pairs) )
    allocate( soap_pol_der(1:n_soap, 1:n_atom_pairs) )
    radial_exp_coeff_der = 0.d0
    angular_exp_coeff_rad_der = 0.d0
    angular_exp_coeff_azi_der = 0.d0
    angular_exp_coeff_pol_der = 0.d0
    cnk_rad_der = 0.d0
    soap_rad_der = 0.d0
    cnk_azi_der = 0.d0
    soap_azi_der = 0.d0
    cnk_pol_der = 0.d0
    soap_pol_der = 0.d0
  else
!   We need this dummy variable defined here. Note the decreased range for the second index
    allocate( radial_exp_coeff_der(1:n_max, 1:1) )
  end if


  if( do_timing )then
    call cpu_time(time2)
    memory_time = memory_time + time2 - time1
    time1 = time2
  end if


  do i = 1, n_species
    if( basis == "poly3gauss" )then
      call get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                        rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                        amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                        scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                        radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                        radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    else if( basis == "poly3" )then
      call get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                   rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                   amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                   scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                   do_central(i), central_weight(i), &
                                                   radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                   radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    end if
  end do



  if( do_timing )then
    call cpu_time(time2)
    radial_time = time2 - time1
    time1 = time2
  end if




! For the angular expansion the masking works differently, since we do not have a species-augmented basis as in the
! radial expansion part.
  call get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, atom_sigma_t, atom_sigma_t_scaling, &
                                          rcut_max, l_max, eimphi, preflm, plm_array, prefl, prefm, &
                                          fact_array, mask, n_species, eimphi_rad_der, &
                                          do_derivatives, prefl_rad_der, angular_exp_coeff, angular_exp_coeff_rad_der, &
                                          angular_exp_coeff_azi_der, angular_exp_coeff_pol_der )
  if( do_timing )then
    call cpu_time(time2)
    angular_time = time2 - time1
  end if
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"





!! For debugging (only gfortran)
!  if( .false. )then
!    write(*,*) "# Gaussian centered at ", rjs(4)
!    write(*,"(A)",advance="no") "rho(x) = "
!    do i = 1, alpha_max
!      write(*,"(A,I0,A,E16.8,A)",advance="no") "p", i, "(x) *", radial_exp_coeff(i,4), "+"
!    end do
!    write(*,*) "0."
!  end if




!  write(*,*) "Obtaining full expansion coefficients..."
!  write(*,*)
  if( do_timing )then
    call cpu_time(time1)
  end if

  if( allocated( starting_index ) )deallocate( starting_index )
  if( allocated( cnk_slice ) )deallocate( cnk_slice )
  allocate( starting_index(1:n_sites) )
  allocate( cnk_slice(1:1 + l_max*(l_max+1)/2 + l_max) )
  k2 = 0
  do i = 1, n_sites
    starting_index(i) = k2
    do j = 1, n_neigh(i)
      k2 = k2 + 1
    end do
  end do

!$omp parallel do private(i, k2, j, n, l, m, k, amplitude, cnk_slice) schedule(static,1)
  do i = 1, n_sites
    k2 = starting_index(i)
    do j = 1, n_neigh(i)
      k2 = k2 + 1
      do n = 1, n_max
        cnk_slice = 0.d0
        do l = 0, l_max
          do m = 0, l
            k = 1 + l*(l+1)/2 + m
!           It is messy with the prefactor in spherical harmonics but we need to be sure because of the central atom below
!            cnk(k, n, i) = cnk(k, n, i) + 4.d0*pi * radial_exp_coeff(n, k2) * angular_exp_coeff(k, k2)
            cnk_slice(k) = cnk_slice(k) + 4.d0*pi * radial_exp_coeff(n, k2) * angular_exp_coeff(k, k2)
          end do
        end do
        cnk(:, n, i) = cnk_slice(:)
      end do
    end do
    do k = 1, species_multiplicity(i)
      j = species(k, i)
      if( basis == "poly3gauss" .and. central_weight(j) /= 0.d0 )then
        if( radial_enhancement == 1 )then
          amplitude = dsqrt(2.d0/pi) * atom_sigma_r(j) / rcut_hard(j)
        else if( radial_enhancement == 2 )then
          amplitude = atom_sigma_r(j)**2 / rcut_hard(j)**2
        else
          amplitude = 1.d0
        end if
        cnk(1, i_beg(j):i_end(j), i) = cnk(1, i_beg(j):i_end(j), i) + &
                                       amplitude * central_weight(j) * dsqrt(4.d0*pi) * pi**0.25d0 * &
                                       dsqrt(atom_sigma_r(j) / 2.d0) * &
                                       rcut_hard(j)**3 / atom_sigma_t(j)**2 / atom_sigma_r(j) * &
                                       matmul(W(i_beg(j):i_end(j),i_beg(j):i_end(j)), S(i_beg(j):i_end(j), i_end(j)) )
      end if
    end do
  end do
!$omp end parallel do

! Do derivatives
  if( do_derivatives )then
! Get the radial derivatives:
    call get_derivatives(radial_exp_coeff, angular_exp_coeff, radial_exp_coeff_der, &
                         angular_exp_coeff_rad_der, angular_exp_coeff_azi_der, angular_exp_coeff_pol_der, &
                         n_sites, n_max, l_max, n_neigh, rjs, rcut_max, cnk_rad_der, cnk_azi_der, cnk_pol_der )
  end if

  if( do_timing )then
    call cpu_time(time2)
    coeff_time = time2 - time1
  end if







!  write(*,*) "Building SOAP vectors..."

  if( do_timing )then
    call cpu_time(time1)
  end if

! Create the multiplicity array
  if( recompute_basis )then
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          if( skip_soap_component(l, np, n) )cycle
          do m = 0, l
            counter2 = counter2 + 1
          end do
        end do
      end do
    end do
    if( allocated(multiplicity_array) )then
      deallocate(multiplicity_array)
    end if
    allocate( multiplicity_array(1:counter2) )
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          if( skip_soap_component(l, np, n) )cycle
          do m = 0, l
            counter2 = counter2 + 1
            multiplicity = 1.d0
            if( n /= np )then
              multiplicity = multiplicity * dsqrt(2.d0)
            end if
            if( m > 0 )then
              multiplicity = multiplicity * 2.d0
            end if
            multiplicity_array(counter2) = multiplicity
          end do
        end do
      end do
    end do
  end if
  recompute_basis = .false.


  do i = 1, n_sites
    counter = 0
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          if( skip_soap_component(l, np, n) )cycle
          counter = counter+1
          do m = 0, l
            k = 1 + l*(l+1)/2 + m
            counter2 = counter2 + 1
!            multiplicity = 1.d0
!            if( n /= np )then
!              multiplicity = multiplicity * dsqrt(2.d0)
!            end if
!            if( m > 0 )then
!              multiplicity = multiplicity * 2.d0
!            end if
            multiplicity = multiplicity_array(counter2)
            soap(counter, i) = soap(counter, i) + multiplicity * real(cnk(k, n, i) * conjg(cnk(k, np, i)))
          end do
        end do
      end do
    end do
    sqrt_dot_p(i) = dsqrt(dot_product(soap(1:n_soap, i), soap(1:n_soap, i)))
!   This is to avoid NaNs when the SOAP sphere is empty
    if( sqrt_dot_p(i) < 1.d-5 )then
      sqrt_dot_p(i) = 1.d0
    end if
  end do


  if( do_derivatives )then
!   Derivatives of the SOAP descriptor in spherical coordinates
!****************************
! Uncomment for detailed timing check

! call cpu_time(time1)
!****************************
    k2 = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        counter = 0
        counter2 = 0
        do n = 1, n_max
          do np = n, n_max
            do l = 0, l_max
              if( skip_soap_component(l, np, n) )cycle
              counter = counter+1
              do m = 0, l
                k = 1 + l*(l+1)/2 + m
                counter2 = counter2 + 1
!                multiplicity = 1.d0
!               These ifs here are slow. I can probably get a 20% speedup in this part of the code (which is the current
!               bottleneck) by optimizing this bit...
!                if( n /= np )then
!                  multiplicity = multiplicity * dsqrt(2.d0)
!                end if
!                if( m > 0 )then
!                  multiplicity = multiplicity * 2.d0
!                end if
                multiplicity = multiplicity_array(counter2)
                soap_rad_der(counter, k2) = soap_rad_der(counter, k2) + multiplicity * real( cnk_rad_der(k, n, k2) * &
                                            conjg(cnk(k, np, i)) + cnk(k, n, i) * conjg(cnk_rad_der(k, np, k2)) )
                soap_azi_der(counter, k2) = soap_azi_der(counter, k2) + multiplicity * real( cnk_azi_der(k, n, k2) * &
                                            conjg(cnk(k, np, i)) + cnk(k, n, i) * conjg(cnk_azi_der(k, np, k2)) )
                soap_pol_der(counter, k2) = soap_pol_der(counter, k2) + multiplicity * real( cnk_pol_der(k, n, k2) * &
                                            conjg(cnk(k, np, i)) + cnk(k, n, i) * conjg(cnk_pol_der(k, np, k2)) )
              end do
            end do
          end do
        end do
!****************************
! Uncomment for detailed timing check
!
!      end do
!    end do
! call cpu_time(time2)
! write(*,*) time2-time1
! call cpu_time(time1)
!    k2 = 0
!    do i = 1, n_sites
!      do j = 1, n_neigh(i)
!        k2 = k2 + 1
!****************************
        soap_rad_der(1:n_soap, k2) = soap_rad_der(1:n_soap, k2) / sqrt_dot_p(i) - &
                                     soap(1:n_soap, i) / sqrt_dot_p(i)**3 * &
                                     dot_product( soap(1:n_soap, i), soap_rad_der(1:n_soap, k2) )
        soap_azi_der(1:n_soap, k2) = soap_azi_der(1:n_soap, k2) / sqrt_dot_p(i) - &
                                     soap(1:n_soap, i) / sqrt_dot_p(i)**3 * &
                                     dot_product( soap(1:n_soap, i), soap_azi_der(1:n_soap, k2) )
        soap_pol_der(1:n_soap, k2) = soap_pol_der(1:n_soap, k2) / sqrt_dot_p(i) - &
                                     soap(1:n_soap, i) / sqrt_dot_p(i)**3 * &
                                     dot_product( soap(1:n_soap, i), soap_pol_der(1:n_soap, k2) )
!       Transform to Cartesian
        if( j == 1 )then
          k3 = k2
        else
          soap_cart_der(1, 1:n_soap, k2) = dsin(thetas(k2)) * dcos(phis(k2)) * soap_rad_der(1:n_soap, k2) - &
                                           dcos(thetas(k2)) * dcos(phis(k2)) / rjs(k2) * soap_pol_der(1:n_soap, k2) - &
                                           dsin(phis(k2)) / rjs(k2) * soap_azi_der(1:n_soap, k2)
          soap_cart_der(2, 1:n_soap, k2) = dsin(thetas(k2)) * dsin(phis(k2)) * soap_rad_der(1:n_soap, k2) - &
                                           dcos(thetas(k2)) * dsin(phis(k2)) / rjs(k2) * soap_pol_der(1:n_soap, k2) + &
                                           dcos(phis(k2)) / rjs(k2) * soap_azi_der(1:n_soap, k2)
          soap_cart_der(3, 1:n_soap, k2) = dcos(thetas(k2)) * soap_rad_der(1:n_soap, k2) + &
                                           dsin(thetas(k2)) / rjs(k2) * soap_pol_der(1:n_soap, k2)
!         MAKE SURE THAT THIS IS CORRECT FOR THE CENTRAL ATOM DERIVATIVES
          soap_cart_der(1, 1:n_soap, k3) = soap_cart_der(1, 1:n_soap, k3) - soap_cart_der(1, 1:n_soap, k2)
          soap_cart_der(2, 1:n_soap, k3) = soap_cart_der(2, 1:n_soap, k3) - soap_cart_der(2, 1:n_soap, k2)
          soap_cart_der(3, 1:n_soap, k3) = soap_cart_der(3, 1:n_soap, k3) - soap_cart_der(3, 1:n_soap, k2)
        end if
      end do
    end do
!****************************
! Uncomment for detailed timing check
!
! call cpu_time(time2)
! write(*,*) time2-time1
!****************************
  end if

! Now we normalize the soap vectors:
  do i = 1, n_sites
    soap(1:n_soap, i) = soap(1:n_soap, i) / sqrt_dot_p(i)
  end do


! This is for debugging
 if( .false. )then
!   open(unit=10, file="soap_desc.dat", status="unknown", access="append")
   open(unit=10, file="soap_desc.dat", status="unknown")
   write(10, *) soap(1:n_soap, 1)
   close(10)
   open(unit=10, file="soap_rad_der.dat", status="unknown", access="append")
   write(10, *) soap_rad_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_azi_der.dat", status="unknown", access="append")
   write(10, *) soap_azi_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_pol_der.dat", status="unknown", access="append")
   write(10, *) soap_pol_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_cart_der.dat", status="unknown")
   do i = 1, 5
     write(10, *) soap_cart_der(1, 1:n_soap, i)
     write(10, *) soap_cart_der(2, 1:n_soap, i)
     write(10, *) soap_cart_der(3, 1:n_soap, i)
     write(10, *)
   end do
   close(10)
!   open(unit=10, file="norm.dat", status="unknown", access="append")
!   write(10, *) sqrt_dot_p(1)
!   close(10)

  end if

!  call cpu_time(time2)
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"


  deallocate( eimphi, preflm, plm_array, prefl, prefm, fact_array, radial_exp_coeff, angular_exp_coeff, cnk, &
              i_beg, i_end, do_central, sqrt_dot_p )
  if( do_derivatives )then
    deallocate( radial_exp_coeff_der, angular_exp_coeff_rad_der, soap_rad_der, soap_azi_der, soap_pol_der, &
                cnk_rad_der, cnk_azi_der, cnk_pol_der, eimphi_rad_der, angular_exp_coeff_azi_der, &
                prefl_rad_der, angular_exp_coeff_pol_der )
  else
    deallocate( radial_exp_coeff_der )
  end if

  if( do_timing )then
    call cpu_time(time2)
    soap_time = time2-time1
    total_time = time2 - time3
    write(*,*)'                                       |'
    write(*,*)'SOAP timings:                          |'
    write(*,*)'                                       |'
    write(*,'(A, F8.3, A)') '  *) Radial expansion: ', radial_time, ' seconds |'
    write(*,'(A, F6.3, A)') '  *) Radial basis build: ', basis_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Angular expansion: ', angular_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Expansion coeffs.: ', coeff_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) SOAP vector build: ', soap_time, ' seconds |'
    if( compress_soap )then
      write(*,'(A, F8.3, A)') '  *) SOAP compression: ', compress_time, ' seconds |'
    end if
    write(*,'(A, F7.3, A)') '  *) Memory allocation: ', memory_time, ' seconds |'
    write(*,'(A, F19.3, A)') '  *) Total: ', total_time, ' seconds |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if

  end subroutine get_soap
!**************************************************************************











!**************************************************************************
  subroutine get_derivatives(radial_exp_coeff, angular_exp_coeff, radial_exp_coeff_der, &
                             angular_exp_coeff_rad_der, angular_exp_coeff_azi_der, angular_exp_coeff_pol_der, &
                             n_sites, n_max, l_max, n_neigh, rjs, cutoff, cnk_rad_der, cnk_azi_der, cnk_pol_der )

    implicit none

    real*8, intent(in) :: radial_exp_coeff(:,:), radial_exp_coeff_der(:,:), rjs(:), cutoff
    complex*16, intent(in) :: angular_exp_coeff(:,:)
    complex*16 :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:), angular_exp_coeff_pol_der(:,:)
    complex*16, intent(out) :: cnk_rad_der(:,:,:), cnk_azi_der(:,:,:), cnk_pol_der(:,:,:)
    integer, intent(in) :: n_sites, n_max, l_max, n_neigh(:)
    integer :: k2, i, j, n, m, k, l
    real*8 :: pi

    pi = dacos(-1.d0)


    k2 = 0
    do i = 1, n_sites
!     We could skip j = 1, which is the central atom
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        if( rjs(k2) < cutoff )then
          do n = 1, n_max
            do l = 0, l_max
              do m = 0, l
                k = 1 + l*(l+1)/2 + m
!               Radial derivatives:
                cnk_rad_der(k, n, k2) = 4.d0*pi * ( angular_exp_coeff(k, k2) * radial_exp_coeff_der(n, k2) + &
                                         angular_exp_coeff_rad_der(k, k2) * radial_exp_coeff(n, k2) )
!               Azimuthal angle derivatives:
                cnk_azi_der(k, n, k2) = 4.d0*pi * angular_exp_coeff_azi_der(k, k2) * radial_exp_coeff(n, k2)
!               Polar angle derivatives:
                cnk_pol_der(k, n, k2) = 4.d0*pi * angular_exp_coeff_pol_der(k, k2) * radial_exp_coeff(n, k2)
              end do
            end do
          end do
        end if
      end do
    end do

    return
  end subroutine
!**************************************************************************







!**************************************************************************
  subroutine assign_species_multiplicity(max_species_multiplicity, species_types_soap, xyz_species, &
                                         xyz_species_supercell, n_species_soap, all_atoms, which_atom, &
                                         indices, n_neigh, n_atom_pairs, neighbors_list, mask_species, &
                                         species, species_multiplicity, species_multiplicity_supercell )
!                                         species, species_multiplicity )

    implicit none

!   Input variables
    integer, intent(in) :: max_species_multiplicity, n_species_soap, which_atom, indices(1:3), &
                           n_neigh(:), n_atom_pairs, neighbors_list(:)
    character*8, intent(in) :: species_types_soap(:), xyz_species(:), xyz_species_supercell(:)
    logical, intent(in) :: all_atoms

!   Output variables
    integer, allocatable, intent(out) :: species(:,:), species_multiplicity(:), &
                                         species_multiplicity_supercell(:)
!    integer, allocatable :: species_multiplicity_supercell(:)
    logical, allocatable, intent(out) :: mask_species(:,:)

!   Internal variables
    integer, allocatable :: species_supercell(:,:)
    integer :: n_sites, n_sites_supercell, i, j, k, i2, j2, k2, ijunk, counter


    n_sites = size(xyz_species)
    n_sites_supercell = size(xyz_species_supercell)

    allocate( species(1:max_species_multiplicity, 1:n_sites) )
    allocate( species_multiplicity(1:n_sites) )
    species = 0
    species_multiplicity = 0

    do i = 1, n_sites
      do j = 1, n_species_soap
        if( xyz_species(i) == species_types_soap(j) )then
          species_multiplicity(i) = species_multiplicity(i) + 1
          species(species_multiplicity(i), i) = j
        end if
      end do
    end do

    if( n_sites_supercell > n_sites )then
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = 0
      species_multiplicity_supercell = 0
!      counter = 0
!      do i2 = 1, indices(1)
!        do j2 = 1, indices(2)
!          do k2 = 1, indices(3)
!            do i = 1, n_sites
!              counter = counter + 1
!              species_supercell(:, counter) = species(:, i)
!            end do
!          end do
!        end do
!      end do
      do i = 1, n_sites_supercell
        do j = 1, n_species_soap
          if( xyz_species_supercell(i) == species_types_soap(j) )then
            species_multiplicity_supercell(i) = species_multiplicity_supercell(i) + 1
            species_supercell(species_multiplicity_supercell(i), i) = j
          end if
        end do
      end do
    else
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = species
      species_multiplicity_supercell = species_multiplicity
    end if

!   This is perhaps not the most efficient way to select only one atom, fix in the future <----- FIX THIS
    if( .not. all_atoms )then
      n_sites = 1
      deallocate( species )
      allocate( species(1:max_species_multiplicity, 1:n_sites) )
      species(1:max_species_multiplicity, 1) = species_supercell(1:max_species_multiplicity, which_atom)
      species_supercell(1:max_species_multiplicity, which_atom) = species_supercell(1:max_species_multiplicity, 1)
      species_supercell(1:max_species_multiplicity, 1) = species(1:max_species_multiplicity, 1)

      ijunk = species_multiplicity(which_atom)
      deallocate( species_multiplicity )
      allocate( species_multiplicity(1:n_sites) )
      species_multiplicity(1) = ijunk
    end if

    allocate( mask_species(1:n_atom_pairs, 1:n_species_soap) )
    mask_species = .false.

    k2 = 0
    do i = 1, n_sites
      do k = 1, n_neigh(i)
        k2 = k2 + 1
        j = neighbors_list(k2)
        if( k == 1 )then
          do i2 = 1, species_multiplicity_supercell(i)
            mask_species(k2, species_supercell(i2, j)) = .true.
          end do
        else
          do i2 = 1, species_multiplicity_supercell(i)
            mask_species(k2, species_supercell(i2, j)) = .true.
          end do
        end if
      end do
    end do

!    deallocate( species_supercell, species_multiplicity_supercell )
    deallocate( species_supercell )


  end subroutine
!**************************************************************************


end module soap_turbo_desc
