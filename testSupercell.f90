program testSupercell
    use SupercellModule
    implicit none

    real(kind=RealKind), save :: temperature

    write(*, "(A)", advance='no') "Enter temperature in Kelvin:"
    read(*,*) temperature
    print *, "temperature: ", temperature

    call supercell_calculation(temperature)

    print *, "Calculation completed successfully!"

end program testSupercell
