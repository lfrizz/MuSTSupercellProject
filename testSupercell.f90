program testSupercell
    use SupercellModule
    !use, intrinsic :: iso_fortran_env, only: dp => real64
    !use KindParamModule, only : IntKind, RealKind, CmplxKind
    !use MathParamModule, only : ZERO, HALF, ONE, TWO, PI
    !use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom
    !use ChemElementModule, only : getDebyeTemperature, getAtomicMass, getAtomicRadius
    !use global_arrays
    implicit none

    !real(kind=RealKind) :: hbar_au = 1.0d0           !hbar in au
    !real(kind=RealKind) :: kB_au = 3.1668115663d-6    !Boltzmann constant in Hartree/K
    !real(kind=RealKind) :: kB_Ry = Boltzmann   !in Rydberg/K
    !real(kind=RealKind) :: electron_mass = 0.5d0     !m_e = 1/2 in au
    !real(kind=RealKind) :: sqrt_2pi = sqrt(2.0d0 * pi)
    real(kind=RealKind), save :: temperature
    !integer(kind=IntKind) :: ios, num_atoms, num_species, i, pos
    !real(kind=RealKind) :: avg_debye_temp, avg_mass, avg_radius

    write(*, "(A)", advance='no') "Enter temperature in Kelvin:"
    read(*,*) temperature
    print *, "temperature: ", temperature

    !call read_position(num_atoms, num_species)
    !call calculate_avg_properties(num_atoms, num_species, avg_debye_temp, avg_mass)
    !call calculate_species_msd(num_atoms, num_species, temperature, avg_debye_temp, avg_mass)
    !call generate_atomic_displacements_with_weights(num_atoms, num_species)

    call supercell_calculation(temperature)

    print *, "Calculation completed successfully!"

end program testSupercell