program lifetimes
    use mod_lammps_reader
    use iso_fortran_env, only: int64, real64
    implicit none

    type(lammps_reader) :: reader
    integer(int64) :: natoms, step, old_step=-1, delta=-1, lifetime
    integer, dimension(:), allocatable :: curr_status
    integer(int64), dimension(:), allocatable :: prev_change

    integer(int64), dimension(5,10000000) :: histogram
    integer, parameter :: H2=1, H=2, D2=3, D=4, HD=5
    integer :: t, i, new_status, idx
    logical :: success
    character(len=50) :: infile

    ! infile(:) = " "
    ! infile = "600K/combined.bin"

    call get_command_argument(1, infile)

    histogram(:,:) = 0

    call reader%open_file(trim(infile))
    call reader%read_step()

    natoms = reader%header%num_atoms
    delta = reader%next_header%step - reader%header%step

    allocate(curr_status(natoms), prev_change(natoms))

    prev_change(:) = reader%header%step
    do i = 1, natoms
        curr_status(i) = get_status(reader%values(:,i))
    end do

    do while (reader%has_next_step)
        ! write(*,*) reader%header%step
        call reader%read_step(success)
        if ((.not. success) .or. reader%header%step==old_step) cycle

        do i = 1, natoms
            new_status = get_status(reader%values(:,i))
            if (new_status /= curr_status(i)) then
                lifetime = reader%header%step - prev_change(i)
                idx = min(lifetime/delta, size(histogram,2)+1)
                histogram(curr_status(i),idx) = histogram(curr_status(i), idx) + 1

                prev_change(i) = reader%header%step
                curr_status(i) = new_status
            end if
        end do
    end do

    do t = 1, size(histogram,2)+1
        write(*,*) t*delta, histogram(:,t)
    end do



contains
    integer function get_status(values) result(s)
        real(real64), intent(in) :: values(:)

        !  id type v_mytype   x   y   z c_coord_H c_coord_D c_coord_Au c_coord_Pd c_coord_PdAu c_coord_Au1 c_coord_Pd1 c_coord_PdAu1 c_coord_Au2 c_coord_Pd2 c_coord_PdAu2
        ! [0]  [1]      [2] [3] [4] [5]       [6]       [7]        [8]        [9]         [10]        [11]        [12]          [13]
        integer :: mytype, coordH, coordD, coordAu, coordPd, coordPdAu
        ! real(real64) :: z

        mytype = values(3)
        coordH = values(7)
        coordD = values(8)
        coordAu = values(9)
        coordPd = values(10)
        coordPdAu = values(11)

        if (mytype==3) then
            if (coordH > 0) then
                s = H2
            else if (coordH == 0 .and. coordD==1) then
                s = HD
            else
                s = H
            end if
        else
            if (coordD > 0) then
                s = D2
            else if (coordH == 1 .and. coordD==0) then
                s = HD
            else
                s = D
            end if
        end if

    end function
end program lifetimes
