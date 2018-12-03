module theis_solution_wrapper

    use iso_c_binding, only: c_double, c_int
    use theis_solution, only: analyticaltheis

    implicit none

    contains

        subroutine c_theis_solution()

end module theis_solution_wrapper