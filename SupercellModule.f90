!   ***************************************************************  
!   Subroutine for calculating the atomic displacements of an alloy using
!   the supercell approach

!   ***************************************************************
    module global_arrays
        implicit none
        character(len=10), allocatable, save :: species_names(:)
        character(len=10), allocatable, save :: species_list(:)
        integer, allocatable, save :: species_counts(:)
        integer, allocatable, save :: species_ids(:)
        real, allocatable, save :: x(:), y(:), z(:)
        real, allocatable, save :: msd_displacement(:)
        real, allocatable, save :: displacements(:,:), displacement_weights(:)
    end module global_arrays

    !program main
    !    use global_arrays
    !    implicit none
    !    call supercell_calculation()
    !end program main

    module SupercellModule
        use global_arrays
        implicit none

    contains

        subroutine supercell_calculation(temperature)

            use, intrinsic :: iso_fortran_env, only: dp => real64
            use KindParamModule, only : IntKind, RealKind, CmplxKind
            use MathParamModule, only : ZERO, HALF, ONE, TWO, PI
            use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom
            use ChemElementModule, only : getDebyeTemperature, getAtomicMass, getAtomicRadius
            use global_arrays

            implicit none

            real(kind=RealKind) :: hbar_au = 1.0d0           !hbar in au
            real(kind=RealKind) :: kB_au = 3.1668115663d-6    !Boltzmann constant in Hartree/K
            real(kind=RealKind) :: kB_Ry = Boltzmann   !in Rydberg/K
            real(kind=RealKind) :: electron_mass = 0.5d0     !m_e = 1/2 in au
            real(kind=RealKind) :: sqrt_2pi = sqrt(2.0d0 * pi)
            !real(kind=RealKind), save :: temperature
            real, intent(in) :: temperature
            integer(kind=IntKind) :: ios, num_atoms, num_species, i, pos
            real(kind=RealKind) :: avg_debye_temp, avg_mass, avg_radius
            real(kind=RealKind) :: temp_8

            temp_8 = real(temperature, kind=RealKind)

            call read_position(num_atoms, num_species)
            call calculate_avg_properties(num_atoms, num_species, avg_debye_temp, avg_mass)
            call calculate_species_msd(num_atoms, num_species, temp_8, avg_debye_temp, avg_mass)
            call generate_atomic_displacements_with_weights(num_atoms, num_species)

        end subroutine supercell_calculation

        function calculate_debye_function(x_val) result(phi)
            use KindParamModule, only : RealKind
            use MathParamModule, only : PI
            use global_arrays
            implicit none

            real(kind=RealKind), intent(in) :: x_val
            real(kind=RealKind) :: phi
            integer, parameter :: n_points = 100
            real(kind=RealKind) :: integral, t, dt
            integer :: i

            t=0.001d0
            integral = 0.0d0
            dt = real(x_val, kind=8)/real(n_points, kind=8)
            print *, "x_val: ", real(x_val, kind=8)
            print *, "n_points: ", real(n_points, kind=8)

            do i=1, n_points-1
                integral = integral + t**3 / (exp(t) - 1.0d0) * dt
                t = t + dt
            end do
            phi = 3.0d0 / x_val**3 * integral
        end function calculate_debye_function

        subroutine read_position(num_atoms, num_species)
            use global_arrays
            implicit none

            integer, intent(out) :: num_atoms, num_species
            integer :: i, j, ios, colon_pos, pos
            character(len=256) :: line, species_name_str
            logical :: found
            logical :: in_positions
            integer :: line_count, species_start_line

            integer, parameter :: unit = 10

            open(unit=unit, file="position.dat", status="old", action="read", iostat=ios)
            if (ios /= 0) then
                print *, "ERROR: Cannot open position.dat file"
                stop
            end if

            num_atoms = 0
            num_species = 0
            species_start_line = 0
            line_count = 0

            print *, "Reading position.dat file..."

            ! FIRST PASS: Read entire file to find counts and species
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line_count = line_count + 1
                
                ! Look for "# Number of medium atoms:"
                if (index(line, "# Number of medium atoms:") > 0) then
                    colon_pos = index(line, ":")
                    read(line(colon_pos+1:), *) num_atoms
                    print *, "Found num_atoms =", num_atoms
                end if
                
                ! Look for "# Number of medium atom type:"
                if (index(line, "# Number of medium atom type:") > 0) then
                    colon_pos = index(line, ":")
                    read(line(colon_pos+1:), *) num_species
                    print *, "Found num_species =", num_species
                end if
                
                ! Find where species information starts (after the medium atom type line)
                if (index(line, "# Number of medium atom type:") > 0) then
                    species_start_line = line_count + 1
                end if
            end do

            if (num_atoms == 0 .or. num_species == 0) then
                print *, "ERROR: Could not find atom/species information in position.dat"
                close(unit)
                stop
            end if

            ! Allocate arrays
            allocate(species_names(num_atoms))
            allocate(x(num_atoms), y(num_atoms), z(num_atoms))
            allocate(species_list(num_species))
            allocate(species_counts(num_species))
            allocate(species_ids(num_atoms))

            rewind(unit)

            ! SECOND PASS: Read species list and counts
            line_count = 0
            i = 0
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line_count = line_count + 1
                
                ! Read species lines
                if (line_count >= species_start_line .and. line_count < species_start_line + num_species) then
                    i = i + 1
                    ! Extract species name - format: "# Number of  Li :       64"
                    pos = index(line, "of ") + 3
                    colon_pos = index(line, ":")
                    species_name_str = adjustl(line(pos:colon_pos-1))
                    species_list(i) = trim(species_name_str)
                    
                    ! Extract the count after colon
                    read(line(colon_pos+1:), *) species_counts(i)
                    
                    !print *, "Species", i, ":", trim(species_list(i)), "count:", species_counts(i)
                end if
            end do

            rewind(unit)

            ! THIRD PASS: Read atom positions
            ! Skip to the position data - find the line after the second separator
            line_count = 0
            i = 0
            in_positions = .false.
            
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line_count = line_count + 1
                
                ! Skip empty lines
                if (trim(line) == "") cycle
                
                ! Check if we've found the separator line before positions
                if (.not. in_positions) then
                    if (index(line, "==============================================================") > 0) then
                        ! Found a separator, now start looking for positions after the next line
                        ! Read the next line (which should be the lattice vectors or empty)
                        read(unit, '(A)', iostat=ios) line
                        line_count = line_count + 1
                        if (ios /= 0) exit
                        
                        ! Read the three lattice vector lines
                        do j = 1, 3
                            read(unit, '(A)', iostat=ios) line
                            line_count = line_count + 1
                            if (ios /= 0) exit
                        end do
                        
                        ! Read the next separator line
                        read(unit, '(A)', iostat=ios) line
                        line_count = line_count + 1
                        if (ios /= 0) exit
                        
                        ! Now we should be at the position data
                        in_positions = .true.
                        cycle
                    end if
                end if
                
                ! Read positions
                if (in_positions) then
                    ! Skip comment lines and empty lines
                    if (line(1:1) == "#" .or. trim(line) == "") cycle
                    
                    i = i + 1
                    if (i <= num_atoms) then
                        read(line, *, iostat=ios) species_names(i), x(i), y(i), z(i)
                        if (ios /= 0) then
                            print *, "ERROR reading position line:", trim(line)
                            cycle
                        end if
                        
                        !if (i <= 5) then
                        !    print *, "Read atom", i, ":", trim(species_names(i)), x(i), y(i), z(i)
                        !end if
                    end if
                end if
            end do

            close(unit)

            if (i == 0) then
                print *, "ERROR: No positions were read from position.dat"
                print *, "Please check the file format"
                stop
            end if

            if (i /= num_atoms) then
                print *, "WARNING: Read", i, "atoms but expected", num_atoms
                num_atoms = i  ! Adjust to what was actually read
            end if

            ! CREATE ENUMERATION
            do i = 1, num_atoms
                found = .false.
                
                do j = 1, num_species
                    if (trim(species_names(i)) == trim(species_list(j))) then
                        species_ids(i) = j
                        found = .true.
                        exit
                    end if
                end do
                
                if (.not. found) then
                    print *, "Error: Unknown species ", trim(species_names(i))
                    print *, "Available species:", (trim(species_list(j)), j=1, num_species)
                    stop
                end if
            end do

            ! OUTPUT
            print *, "Total atoms:", num_atoms
            print *, "Number of species:", num_species
            
            print *, "Species mapping:"
            do i = 1, num_species
                print *, i, " -> ", trim(species_list(i))
            end do
            
            !print *, "First few atoms (name, id, x, y, z):"
            !do i = 1, min(5, num_atoms)
            !    print *, trim(species_names(i)), species_ids(i), x(i), y(i), z(i)
            !end do

        end subroutine read_position

        subroutine calculate_avg_properties(num_atoms, num_species, avg_debye_temp, avg_mass)
            use KindParamModule, only : RealKind
            use ChemElementModule, only : getDebyeTemperature, getAtomicMass
            use global_arrays
            implicit none

            integer, intent(in) :: num_atoms, num_species
            real(kind=RealKind), intent(out) :: avg_debye_temp, avg_mass
            integer :: i
            real :: total_debye, total_mass

            total_debye = 0.0d0
            total_mass = 0.0d0

            do i = 1, num_species
                total_debye = total_debye + getDebyeTemperature(species_list(i)) * species_counts(i)
                total_mass = total_mass + getAtomicMass(species_list(i)) * 1822.89 * species_counts(i)
                !print *, "Atomic mass(i): ", getAtomicMass(species_list(i)) * 1822.89
                !print *, "species counts(i): ", species_counts(i)
            end do

            avg_debye_temp = total_debye / num_atoms
            avg_mass = total_mass / num_atoms

            write(*,*)
            write(*,*) "Alloy Properties:"
            write(*,'(A,F10.2,A)') "  Average Debye temperature: ", avg_debye_temp, " K"
            write(*,'(A,F10.2,A)') "  Average mass: ", avg_mass, " m_e"
            write(*,'(A,I8)') "  Total atoms: ", num_atoms

        end subroutine calculate_avg_properties

        subroutine calculate_species_msd(num_atoms, num_species, temperature, avg_debye_temp, avg_mass)
            use ChemElementModule, only : getDebyeTemperature, getAtomicMass, getAtomicRadius
            use KindParamModule, only: RealKind
            use MathParamModule, only : PI
            use PhysParamModule, only : Boltzmann
            use global_arrays
            implicit none

            integer, intent(in) :: num_atoms, num_species
            real(kind=RealKind), intent(in) :: temperature, avg_debye_temp, avg_mass
            integer :: i
            real(kind=RealKind) :: debye_func, msd_thermal

            allocate(msd_displacement(num_species))

            do i = 1, num_species
                debye_func = calculate_debye_function(getDebyeTemperature(species_list(i)) / temperature)
                
                msd_thermal = (3.0d0)/(4.0d0 * PI**2*getAtomicMass(species_list(i)) * 1822.89 * Boltzmann/2 * getDebyeTemperature(species_list(i))) * (debye_func / (getDebyeTemperature(species_list(i)) / temperature ) + 0.25d0)
                msd_displacement(i) = msd_thermal
            end do

        end subroutine calculate_species_msd

        subroutine generate_atomic_displacements_with_weights(num_atoms, num_species)
            use global_arrays
            use MathParamModule, only : PI
            implicit none

            integer, intent(in) :: num_atoms, num_species
            real :: u1, u2, r, theta
            real :: gauss1, gauss2, gauss3, sigma_thermal, sigma_chemical
            real :: displacement_thermal(3), displacement_chemical(3)
            real :: weight_thermal, weight_chemical, total_weight, weight_count
            real :: sqrt_2pi = sqrt(2.0d0 * PI)
            integer :: i, species_idx

            allocate(displacements(3, num_atoms))
            allocate(displacement_weights(num_atoms))

            weight_count = 0

            do i=1, num_atoms
                species_idx = species_ids(i)
                sigma_thermal = sqrt(msd_displacement(species_idx)/3.0d0)
                sigma_chemical = sqrt(msd_displacement(species_idx)/3.0d0)
                if((i < num_atoms-1 .AND. species_names(i) /= species_names(i+1)) .OR. i==num_atoms - 1) then
                    print *, "species name: ", species_names(i)
                    print *, "sigma_thermal: ", sigma_thermal
                end if

                call box_muller_pair(gauss1, gauss2)
                call box_muller_pair(gauss3, r)

                displacement_thermal(1) = gauss1 * sigma_thermal
                displacement_thermal(2) = gauss2 * sigma_thermal
                displacement_thermal(3) = gauss3 * sigma_thermal


                weight_thermal = (1.0d0/(sqrt_2pi*sigma_thermal))**3 * &
                            exp(-0.5d0 * dot_product(displacement_thermal, displacement_thermal) / &
                                sigma_thermal**2)

                call box_muller_pair(gauss1, gauss2)
                call box_muller_pair(gauss3, r)
                displacement_chemical(1) = gauss1 * sigma_chemical
                displacement_chemical(2) = gauss2 * sigma_chemical
                displacement_chemical(3) = gauss3 * sigma_chemical

                weight_chemical = (1.0d0/(sqrt_2pi*sigma_chemical))**3 * &
                                exp(-0.5d0 * dot_product(displacement_chemical, displacement_chemical) / &
                                    sigma_chemical**2)

                displacements(:, i) = displacement_thermal + displacement_chemical
                
                total_weight = weight_thermal + weight_chemical
                displacement_weights(i) = total_weight
                weight_count = weight_count + total_weight
            end do
            displacement_weights = displacement_weights / weight_count
            print *, "Sum of weights: ", sum(displacement_weights)

            call write_outfile(num_atoms, num_species)

        end subroutine generate_atomic_displacements_with_weights

        subroutine box_muller(gauss1, gauss2)
            use MathParamModule, only : PI
            implicit none
            
            real, intent(out) :: gauss1, gauss2
            real :: u1, u2, r, theta
            
            ! Generate uniform random numbers in (0,1)
            call random_number(u1)
            call random_number(u2)
            
            ! Ensure we don't get exactly 0 or 1
            u1 = max(min(u1, 0.9999999999d0), 0.0000000001d0)
            u2 = max(min(u2, 0.9999999999d0), 0.0000000001d0)
            
            ! Box-Muller transform
            r = sqrt(-2.0d0 * log(u1))
            theta = 2.0d0 * PI * u2
            gauss1 = r * cos(theta)
            gauss2 = r * sin(theta)
            
        end subroutine box_muller

        subroutine box_muller_pair(gauss1, gauss2)
            use global_arrays
            use MathParamModule, only : PI
            implicit none
            
            real, intent(out) :: gauss1, gauss2
            real :: u1, u2, r, theta
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, 1.0d-10)
            u2 = max(u2, 1.0d-10)

            r = sqrt(-2.0d0 * log(u1))
            theta = 2.0d0 * PI * u2
            gauss1 = r * cos(theta)
            gauss2 = r * sin(theta)
        end subroutine box_muller_pair

        subroutine write_outfile(num_atoms, num_species)
            use global_arrays, only : displacements, displacement_weights, species_list, species_counts, species_names, x, y, z
            implicit none

            integer, intent(in) :: num_atoms, num_species
            integer :: i, ios

            open(newunit=ios, file="position_displacement.dat", status="replace", action="write")
            write(ios, '(A)') "# Units:: atomic units"
            write(ios, '(A)') "  3.50000000"
            write(ios, '(A)') "# If Angstroms is needed, comment out the line above and uncomment the following line"
            write(ios, '(A)') "#   1.85209500"
            write(ios, '(A)') "# =============================================================="
            write(ios, '(A)') "        4.00000000000      0.00000000000      0.00000000000"
            write(ios, '(A)') "        0.00000000000      4.00000000000      0.00000000000"
            write(ios, '(A)') "        0.00000000000      0.00000000000      4.00000000000"
            write(ios, '(A)') "# =============================================================="
            write(ios, '(A)') "# Number of clusters:        0"
            write(ios, '(A,I8)') "# Number of medium atoms:          ", num_atoms
            write(ios, '(A,I8)') "# Number of medium atom type:        ", num_species
            do i=1, num_species
                write(ios, '(A,A,A,I8)') "# Number of  ", trim(species_list(i)), ":       ", species_counts(i)
            end do
            write(ios, '(A)') "# =============================================================="
            
            ! Write original positions plus displacements
            do i=1, num_atoms
                write(ios, '(A,3(2X,F12.8))') trim(species_names(i)), x(i) + displacements(1,i), y(i) + displacements(2,i), z(i) + displacements(3,i)
            end do
            close(ios)

            open(newunit=ios, file="displacement.dat", status="replace", action="write")
            ! Write just displacements
            do i=1, num_atoms
                write(ios, '(A,2X,F24.20)') trim(species_names(i)), (displacements(1,i)**2 + displacements(2,i)**2 + displacements(3,i)**2)**(0.5d0)/3 !* displacement_weights(i) !note - result with weight gives wrong answer
            end do
            close(ios)
        end subroutine write_outfile

    end module SupercellModule