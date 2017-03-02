! ==============================================================================
! wrf_storm_module 
!
! Module contains routines for returning wind and
! pressure fields based on data files from WRF output.
!
! Both NetCDF and ASCII files will be supported.
! TODO: Add NetCDF support.
!
! The ASCII data files are expected to be of the form:
! *lat.dat
! *lon.dat
! *u10.dat
! *v10.dat
! *pmsl.dat
! where * is specified in setrun.py or empty by default.
! 
! ==============================================================================
module wrf_storm_module

    implicit none
    save

    ! WRF storm type definition
    ! Specified wind & pressure field 
    type wrf_storm_type
        ! Size of spatial grids,
        !  corresponds to lengths of lat/lon arrays
        integer :: num_lats
        integer :: num_lons

        ! Location of storm field values
        ! longitude and latitude arrays
        ! start from SW corner
        real(kind=8), allocatable :: lat(:)
        real(kind=8), allocatable :: lon(:)

        ! We keep two time snapshots in memory 
        !  for interpolation, with t_next > t > t_prev
        real(kind=8) :: t_next, t_prev
        real(kind=8), allocatable :: u_next(:,:)
        real(kind=8), allocatable :: v_next(:,:)
        real(kind=8), allocatable :: p_next(:,:)
        real(kind=8), allocatable :: u_prev(:,:)
        real(kind=8), allocatable :: v_prev(:,:)
        real(kind=8), allocatable :: p_prev(:,:)
        ! The values will be updated incrementally:
        !  when t>t_next, replace older states with newer states
        !  and read in next snapshot.
        
        ! These will be used during the interpolation step
        !  for times between t_prev and t_next
        real(kind=8) :: t
        real(kind=8), allocatable :: u(:,:)
        real(kind=8), allocatable :: v(:,:)
        real(kind=8), allocatable :: p(:,:)

        ! The estimate center of the storm 
        !  is also stored for interpolation
        integer :: eye_next(2)
        integer :: eye_prev(2)

        ! Keep track of how many snapshots have been read
        integer :: last_storm_index

        ! for compatibility with parent storm module (?)
        ! Storm physics
        real(kind=8) :: ambient_pressure = 101.3d3 ! 101300 Pascals
        real(kind=8) :: rho_air = 1.3d0

        ! Store the storm data file for repeated reading
        character(len=4096) :: data_path_root

    end type wrf_storm_type

    ! Internal tracking variables for storm
    logical, private :: DEBUG = .true. 
    !logical, private :: DEBUG = .false. 

    ! Tolerance for floating point inequalities
    real(kind=8), parameter :: eps = 1.0e-8 


contains

    ! Setup routine for the WRF storm model
    ! Open data files, get parameters, allocate memory,
    !  and read in first two time snapshots of storm data
    subroutine set_wrf_storm(storm_data_path_root, storm, log_unit)

        use geoclaw_module, only: coordinate_system
        use amr_module, only: t0

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path_root
        type(wrf_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: l_file = 701
        integer :: i, j, io_status, num_lats, num_lons
        character(len=4096) :: storm_data_path

        ! Reading buffer variables
        character(len=100) :: dummy_read
        integer :: readsize

        ! Storm type only works on lat-long coordinate systems
        if (coordinate_system /= 2) then
            stop "explicit storm type does only works on lat-long coordinates."
        endif

        ! We need to count two things:
        !   number of latitude coords (ny)
        !   number of longitude coords (nx)
        ! We don't need number of time snapshots (nt)
        !   since we will load snapshots only as needed
        ! The datafiles for lat & lon contain nx*ny values
        ! The datafiles for u, v, p contain nx*ny*nt values

        ! Open latitudes data file
        if (present(storm_data_path_root)) then
            storm%data_path_root = storm_data_path_root
        else
            storm%data_path_root = './'
        endif
        storm_data_path = trim(storm%data_path_root) // "lat.dat"
        print *,'Reading latitudes data file ',trim(storm_data_path)
        open(unit=l_file,file=storm_data_path,status='old', &
                action='read',iostat=io_status)
        if (io_status /= 0) then
            print *, "Error opening latitudes data file. status = ", io_status
            stop 
        endif            

        ! Count number of data columns
        num_lons = 0
        do ! incremental read of lat file
            ! Read a chunk of the first line
            read(l_file,'(A)', ADVANCE='NO', SIZE=readsize, iostat=io_status) dummy_read
            ! Count number of entries in first record (using . as proxy)
            do i = 1, readsize ! equal to len(dummy_read) unless EOR reached
                if ( dummy_read(i:i) == '.'  ) Then
                    num_lons = num_lons + 1
                endif
            end do
            ! Exit loop if we ran into EOR or EOF
            if (io_status /= 0) then
                exit
            endif
        end do
        storm%num_lons = num_lons

        ! Count number of data rows
        rewind(l_file) 
        num_lats = 0
        do
            read (l_file, *, iostat=io_status)
            ! Exit loop if we ran into an error or we reached the end of the file
            if (io_status /= 0) then
                exit
            else
                num_lats = num_lats + 1
            endif
        end do
        rewind(l_file) ! closed later
        storm%num_lats = num_lats

        ! save lat & lon coords, times, and u/v/p fields
        allocate(storm%lat(num_lats))
        allocate(storm%lon(num_lons))
        allocate(storm%u_prev(num_lons,num_lats))
        allocate(storm%v_prev(num_lons,num_lats))
        allocate(storm%p_prev(num_lons,num_lats))
        allocate(storm%u_next(num_lons,num_lats))
        allocate(storm%v_next(num_lons,num_lats))
        allocate(storm%p_next(num_lons,num_lats))
        allocate(storm%u(num_lons,num_lats))
        allocate(storm%v(num_lons,num_lats))
        allocate(storm%p(num_lons,num_lats))

        ! Read in latitude coords
        do i = 1, num_lats
            read (l_file, *, iostat=io_status) storm%lat(i)
            if (io_status /= 0) exit
        end do
        close(l_file)

        ! Open longitudes data file
        storm_data_path = trim(storm%data_path_root) // "lon.dat"
        print *,'Reading longitudes data file ',trim(storm_data_path)
        open(unit=l_file,file=storm_data_path,status='old', &
                action='read',iostat=io_status)
        if (io_status /= 0) then
            print *, "Error opening longitudes data file. status = ", io_status
            stop 
        endif            
        ! Read in longitude coords. Just need first line.
        read (l_file, *, iostat=io_status) storm%lon
        close(l_file)

        ! This is used to speed up searching for correct storm data
        !  (using ASCII datafiles)
        storm%last_storm_index = 0

        ! Read in the first storm data snapshot as 'next'
        !  and increment storm%last_storm_index to 1
        call read_wrf_storm(storm,t0)

        ! Check if starting time of simulation
        !  is before the first storm data snapshot
        if (t0 < storm%t_next - eps) then
            print *, "Simulation start time precedes storm data. Using clear skies."
            if (DEBUG) print *, "t0=", t0, "first storm t:",storm%t_next
            storm%t_prev = t0
            storm%u_prev = 0
            storm%v_prev = 0
            ! Ambient pressure may not be properly set in parallel
            !storm%ambient_pressure = 101.3d3 ! 101300 Pascals
            storm%p_prev = storm%ambient_pressure
            storm%eye_prev = storm%eye_next
        else
            ! Read in the second storm data snapshot as 'next',
            !  update 'prev' with old 'next' data,
            !  and increment storm%last_storm_index to 2
            call read_wrf_storm(storm,t0)
        endif

        ! Initialize current storm module data
        storm%t = t0
        ! Interpolate wind & pressure fields using prev and next snapshots
        call storm_interpolate(storm)

    end subroutine set_wrf_storm


    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function date_to_seconds(year,months,days,hours,minutes, &
                                               seconds) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: year, months, days, hours, minutes
        real(kind=8), intent(in) :: seconds

        ! Local storage
        integer :: total_days

        ! Count number of days
        total_days = days

        ! Add days for months that have already passed
        if (months > 1) total_days = total_days + 31
        if (months > 2) then
            if (int(year / 4) * 4 == year) then
                total_days = total_days + 29
            else
                total_days = total_days + 28
            endif
        endif
        if (months > 3)  total_days = total_days + 30
        if (months > 4)  total_days = total_days + 31
        if (months > 5)  total_days = total_days + 30
        if (months > 6)  total_days = total_days + 31
        if (months > 7)  total_days = total_days + 30
        if (months > 8)  total_days = total_days + 31
        if (months > 9)  total_days = total_days + 30
        if (months > 10) total_days = total_days + 31
        if (months > 11) total_days = total_days + 30

        ! Convert everything to seconds since the beginning of the year
        time = real((total_days - 1) * 86400 + hours * 3600 + minutes * 60,kind=8)
        time = time + seconds

    end function date_to_seconds

    ! ==========================================================================
    !  real(kind=8) pure ymdh_to_seconds(ymdh)
    !    Convert time from year, month, day, hour (YEARMODAHR)
    !    to seconds since the beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function ymdh_to_seconds(ymdh) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: ymdh

        ! Local storage
        integer :: year, month, day, hour

        ! collect datetime (input format is YEARMODAHR)
        hour = mod(ymdh,100)
        day = mod(ymdh/100,100)
        month = mod(ymdh/10000,100)
        year = ymdh / 1000000

        ! process as seconds from start of year
        time = date_to_seconds(year,month,day,hour,0,0.d0)

    end function ymdh_to_seconds

    ! ==========================================================================
    !  wrf_storm_location(t,storm)
    !    Returns location of hurricane at the current time
    ! ==========================================================================
    function wrf_storm_location(t,storm) result(location)

        use amr_module, only: rinfinity

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(wrf_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Local storage
        real(kind=8) :: alpha
        integer :: eye(2)

        ! Estimate location based on two nearest snapshots 

        ! If no data at this time, return infinity
        if ((t < storm%t_prev - eps) .OR. (t > storm%t_next + eps)) then
            location = [rinfinity,rinfinity]
        else if ((storm%eye_prev(1) == 0) .AND. (storm%eye_next(1) == 0)) then
            location = [rinfinity,rinfinity]
        else
            ! Otherwise check if there is a low pressure system
            !  and if so interpolate eye location from snapshots
            if (storm%eye_prev(1) == 0) then
                eye = storm%eye_next
            else if (storm%eye_next(1) == 0) then
                eye = storm%eye_prev
            else
                ! Determine the linear interpolation parameter (in time)
                if (storm%t_next-storm%t_prev < eps) then
                    print *, "t_next = ", storm%t_next,"t_prev = ", storm%t_prev
                    print *, "t = ", t, "storm%t = ", storm%t
                    alpha = 0
                else
                    alpha = (t-storm%t_prev) / (storm%t_next-storm%t_prev)
                endif
                ! Estimate location index of storm center at time t
                eye = storm%eye_prev + NINT((storm%eye_next - storm%eye_prev) * alpha)
            endif
            ! Convert to lat-lon
            location(1) = storm%lon(eye(1))
            location(2) = storm%lat(eye(2))
        endif

    end function wrf_storm_location

    ! ==========================================================================
    !  read_wrf_storm_data_file()
    !    Opens storm data file and reads next storm entry
    !    Currently only for ASCII file
    ! ==========================================================================
    subroutine read_wrf_storm_file(data_path,storm_array,num_lats,last_storm_index,timestamp)

        implicit none

        ! Subroutine I/O
        real(kind=8), intent(in out) :: storm_array(:,:)
        character(len=*), intent(in) :: data_path
        integer, intent(in) :: num_lats, last_storm_index
        integer, intent(inout) :: timestamp

        ! Local storage
        integer :: j, k, iostatus
        integer, parameter :: data_file = 701

        ! Open the input file
        !
        open(unit=data_file,file=data_path,status='old', &
                action='read',iostat=iostatus)
        if (iostatus /= 0) then
            print *, "Error opening data file: ",trim(data_path)
            print *, "Status = ", iostatus
            stop 
        endif            
        ! Advance to the next time step to be read in
        ! Skip entries based on total number previously read
        do k = 1, last_storm_index
            do j = 1, num_lats
                read(data_file, *, iostat=iostatus)
                ! Exit loop if we ran into an error or we reached the end of the file
                if (iostatus /= 0) then
                    print *, "Unexpected end-of-file reading ",trim(data_path)
                    print *, "Status = ", iostatus
                    if (DEBUG) print *, "k, laststormindex = ", k, last_storm_index
                    if (DEBUG) print *, "j, num_lats = ", j, num_lats
                    timestamp = -1
                    close(data_file) 
                    return
                endif
            enddo
        enddo
        ! Read in next time snapshot 
        do j = 1, num_lats
            read(data_file, *, iostat=iostatus) timestamp, storm_array(:,j) 
            ! Exit loop if we ran into an error or we reached the end of the file
            if (iostatus /= 0) then
                print *, "Unexpected end-of-file reading ",trim(data_path)
                print *, "Status = ", iostatus
                if (DEBUG) print *, "j, num_lats = ", j, num_lats
                !if (DEBUG) print *, "storm_array(:,",j,") = ", storm_array(:,j) 
                timestamp = -1
                close(data_file) 
                return
            endif
        enddo

        timestamp = ymdh_to_seconds(timestamp)
        close(data_file) 

    end subroutine read_wrf_storm_file

    ! ==========================================================================
    !  read_wrf_storm_data()
    !    Reads storm fields for next time snapshot
    !    Currently only for ASCII files
    ! ==========================================================================

    subroutine read_wrf_storm(storm,t)

        implicit none

        ! Subroutine I/O
        type(wrf_storm_type), intent(in out) :: storm
        real(kind=8), intent(in) :: t

        ! Local storage
        character(len=4096) :: data_path
        real(kind=8) :: lowest_p

        ! Reading buffer variables
        integer :: timestamp

        ! Overwrite older storm states with newer storm states
        storm%t_prev = storm%t_next
        storm%u_prev = storm%u_next 
        storm%v_prev = storm%v_next 
        storm%p_prev = storm%p_next 
        storm%eye_prev = storm%eye_next
        
        ! Current time t currently unused in favor of storm%last_storm_index.
        ! This should probably be changed in the future.

        ! Read the u-velocity file
        data_path = trim(storm%data_path_root) // "u10.dat"
        call read_wrf_storm_file(data_path,storm%u_next,storm%num_lats,storm%last_storm_index,timestamp)
        ! Error handling: set to clear skies if file ended
        if (timestamp == -1) then
            storm%u_next = 0
            storm%t_next = storm%t_next + 365*24*60*60
        else
            ! Save timestamp (sec) of next snapshot
            storm%t_next = timestamp
        endif

        ! Read v-velocity file
        data_path = trim(storm%data_path_root) // "v10.dat"
        call read_wrf_storm_file(data_path,storm%v_next,storm%num_lats,storm%last_storm_index,timestamp)
        ! Error handling: set to clear skies if file ended
        if (timestamp == -1) then
            storm%v_next = 0
        endif

        ! Read pressure file
        data_path = trim(storm%data_path_root) // "pmsl.dat"
        call read_wrf_storm_file(data_path,storm%p_next,storm%num_lats,storm%last_storm_index,timestamp)
        ! Error handling: set to clear skies if file ended
        if (timestamp == -1) then
            storm%p_next = storm%ambient_pressure ! causes SIGSEGV - module init not threadsafe?
            !storm%p_next = 101.3d3 ! workaround 
            !if (DEBUG) print *, "ambient pressure = ", storm%ambient_pressure
            !if (DEBUG) print *, "shape(p_next) = ", shape(storm%p_next) 
            !if (DEBUG) print *, "shape(ambient_pressure) = ", shape(storm%ambient_pressure) 
            !if (DEBUG) print *, "shape(timestamp) = ", shape(timestamp) 
            !if (DEBUG) print *, "p_next(:,:) = ", storm%p_next(:,:) 
            storm%eye_next = [0,0]
        else
            ! Convert pressure units: mbar to Pa
            storm%p_next = storm%p_next * 1.0e2
            ! Estimate storm center location based on lowest pressure
            ! (only the array index is saved)
            storm%eye_next = MINLOC(storm%p_next)
            ! If no obvious low pressure area, set storm center to 0 instead
            lowest_p = storm%p_next(storm%eye_next(1),storm%eye_next(2))
            if (lowest_p > storm%ambient_pressure*0.99) then
                storm%eye_next = [0,0]
            endif
        endif

        ! Update number of storm snapshots read in
        storm%last_storm_index = storm%last_storm_index + 1
        if (DEBUG) print *, "last_storm_index=", storm%last_storm_index

    end subroutine read_wrf_storm

    ! ==========================================================================
    !  integer pure get_lat_index(lat)
    !    Returns index of latitude array of the storm data
    !    corresponding to input lat.
    ! ==========================================================================
    pure integer function get_lat_index(lat,storm) result(i)
      
        implicit none

        ! Input
        real(kind=8), intent(in) :: lat
        type(wrf_storm_type), intent(in) :: storm

        ! Local storage
        real(kind=8) :: dy

        ! Out-of-bound conditions:
        if (lat < storm%lat(1)) then
            i = 1
        else if (lat > storm%lat(storm%num_lats)) then
            i = storm%num_lats
        else
            ! Find spacing between latitude values
            dy = (storm%lat(storm%num_lats) - storm%lat(1)) / (storm%num_lats-1)
            ! Determine index based on spacing
            i = 1 + (lat - storm%lat(1)) / dy
        endif

    end function get_lat_index

    ! ==========================================================================
    !  integer pure get_lon_index(lon)
    !    Returns index of longitude array of the storm data
    !    corresponding to input lon.
    ! ==========================================================================
    pure integer function get_lon_index(lon,storm) result(i)
      
        implicit none

        ! Input
        real(kind=8), intent(in) :: lon
        type(wrf_storm_type), intent(in) :: storm

        ! Local storage
        real(kind=8) :: dx

        ! Out-of-bound conditions:
        if (lon < storm%lon(1)) then
            i = 1
        else if (lon > storm%lon(storm%num_lons)) then
            i = storm%num_lons
        else
            ! Find spacing between longitude values
            dx = (storm%lon(storm%num_lons) - storm%lon(1)) / (storm%num_lons-1)
            ! Determine index based on spacing
            i = 1 + (lon - storm%lon(1)) / dx
        endif

    end function get_lon_index

    ! ==========================================================================
    !  storm_interpolate()
    !  Determines intermediate storm values
    !   for time t between t_prev and t_next.
    !  If distinct storms are present at both times,
    !   the storm centers are shifted to an intermediate point
    !   and a weighted average is taken.
    !  Otherwise, no shift occurs and the
    !   weighted average is performed in place.
    ! ==========================================================================
    subroutine storm_interpolate(storm)

        implicit none

        ! Storm description, need "in out" here since will update the storm
        ! values at time t
        type(wrf_storm_type), intent(in out) :: storm

        ! Check if there are distinct storm "eyes"
        ! If not, interpolate wind & pressure fields in place.
        ! If so, spatially shift storm snapshots then interpolate. 
        if (storm%eye_prev(1) == 0 .or. storm%eye_next(1) == 0) then
            call storm_inplace_interpolate(storm)
        else
            call storm_shift_interpolate(storm)
            ! Optional: no weighted average, only shift prev snapshot in space
            !call storm_shift_only(storm)
        endif

    end subroutine storm_interpolate

    ! ==========================================================================
    !  storm_inplace_interpolate()
    !  Determines intermediate storm values
    !   for time t between t_prev and t_next
    !   based on simple weighted average.
    ! ==========================================================================
    subroutine storm_inplace_interpolate(storm)

        implicit none

        ! Storm description, need "in out" here since will update the storm
        ! values at time t
        type(wrf_storm_type), intent(in out) :: storm

        ! Local storage
        real(kind=8) :: alpha

        ! This is just a simple weighted average.
        ! Note that this a poor approach for a tropical cyclone:
        !  intensity is smoothed out between intervals
        !  so intermediate values may appear less intense
        ! For a more realistic storm field, use storm_shift_interpolate()

        ! Determine the linear interpolation parameter (in time)
        alpha = (storm%t-storm%t_prev) / (storm%t_next-storm%t_prev)

        ! Take weighted average of two storm fields
        storm%u = storm%u_prev + &
                (storm%u_next - storm%u_prev) * alpha
        storm%v = storm%v_prev + &
                (storm%v_next - storm%v_prev) * alpha
        storm%p = storm%p_prev + &
                (storm%p_next - storm%p_prev) * alpha

    end subroutine storm_inplace_interpolate

    ! ==========================================================================
    !  storm_shift_interpolate()
    !  Determines intermediate storm values
    !   for time t between t_prev and t_next
    !   both in time (linearly) and in space (approximate) 
    ! ==========================================================================
    subroutine storm_shift_interpolate(storm)

        implicit none

        ! Storm description, need "in out" here since will update the storm
        ! values at time t
        type(wrf_storm_type), intent(in out) :: storm

        ! Local storage
        real(kind=8) :: alpha
        integer :: i,j
        integer :: pi,pj,ni,nj
        integer :: prev_shift(2), next_shift(2)

        ! Determine the linear interpolation parameter (in time)
        alpha = (storm%t-storm%t_prev) / (storm%t_next-storm%t_prev)

        ! Estimate relative location of storm center at time t
        ! Note: The spatial grid is constant in time
        !  so we don't translate to lat-long
        prev_shift = NINT((storm%eye_next - storm%eye_prev) * alpha)
        next_shift = NINT((storm%eye_next - storm%eye_prev) * (alpha - 1))

        ! Now shift the two storm fields onto the intermediate
        !  storm center and use time-weighted average of their values.
        do j = 1,storm%num_lats
            ! If index would be out of bounds, use edge value
            pj = MIN(MAX(1,j-prev_shift(2)),storm%num_lats)
            nj = MIN(MAX(1,j-next_shift(2)),storm%num_lats)
            do i = 1,storm%num_lons
                ! If index would be out of bounds, use edge value
                pi = MIN(MAX(1,i-prev_shift(1)),storm%num_lons)
                ni = MIN(MAX(1,i-next_shift(1)),storm%num_lons)
                ! Perform shift & interpolate
                storm%u(i,j) = storm%u_prev(pi,pj) + &
                    (storm%u_next(ni,nj)-storm%u_prev(pi,pj)) * alpha
                storm%v(i,j) = storm%v_prev(pi,pj) + &
                    (storm%v_next(ni,nj)-storm%v_prev(pi,pj)) * alpha
                storm%p(i,j) = storm%p_prev(pi,pj) + &
                    (storm%p_next(ni,nj)-storm%p_prev(pi,pj)) * alpha
            enddo
        enddo
                

    end subroutine storm_shift_interpolate

    ! ==========================================================================
    !  storm_shift_only()
    !  Determines intermediate storm values
    !   for time t between t_prev and t_next
    !   by shifting storm data towards next position
    !  By not taking averages, this preserves large values
    ! ==========================================================================
    subroutine storm_shift_only(storm)

        implicit none

        ! Storm description, need "in out" here since will update the storm
        ! values at time t
        type(wrf_storm_type), intent(in out) :: storm

        ! Local storage
        real(kind=8) :: alpha
        integer :: i,j
        integer :: pi,pj
        integer :: prev_shift(2)

        ! Determine the linear interpolation parameter (in time)
        alpha = (storm%t-storm%t_prev) / (storm%t_next-storm%t_prev)

        ! Estimate relative location of storm center at time t
        ! Note: The spatial grid is constant in time
        !  so we don't translate to lat-long
        prev_shift = NINT((storm%eye_next - storm%eye_prev) * alpha)

        ! Now shift the earlier storm field
        !  onto the intermediate storm center
        do j = 1,storm%num_lats
            ! If index would be out of bounds, use edge value
            pj = MIN(MAX(1,j-prev_shift(2)),storm%num_lats)
            do i = 1,storm%num_lons
                ! If index would be out of bounds, use edge value
                pi = MIN(MAX(1,i-prev_shift(1)),storm%num_lons)
                ! Perform shift
                storm%u(i,j) = storm%u_prev(pi,pj)
                storm%v(i,j) = storm%v_prev(pi,pj)
                storm%p(i,j) = storm%p_prev(pi,pj)
            enddo
        enddo
                
    end subroutine storm_shift_only

    ! ==========================================================================
    !  set_wrf_storm_fields()
    ! ==========================================================================
    subroutine set_wrf_storm_fields(maux,mbc,mx,my,xlower, &
                                    ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, storm)

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need "in out" here since we may update the storm
        ! if at next time point
        type(wrf_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y
        integer :: i,j,k,l

        if (t < storm%t_prev - eps) then
            print *, "Simulation time precedes storm data in memory. &
                        Race condition?"
            print *, "t=",t,"< t_prev=",storm%t_prev,"t_next=",storm%t_next
        endif

        if (t > storm%t_next + eps) then
            ! Load two snapshots into memory, at times t_next and t_prev
            !$OMP CRITICAL (READ_STORM)
            do while (t > storm%t_next + eps)
            ! update all storm data, including value of t_next
                if (DEBUG) print *,"loading new storm snapshot ",&
                                    "t=",t,"old t_next=",storm%t_next
                call read_wrf_storm(storm,t)
                if (DEBUG) print *,"new t_next=",storm%t_next
                ! If storm data ends, the final storm state is used.
            enddo
            !$OMP END CRITICAL (READ_STORM)
        endif
        
        ! Interpolate storm data in time
        ! t_prev <= t <= t_next
        if (t > storm%t + eps) then
            !$OMP CRITICAL (INTERP_STORM)
            if (t > storm%t + eps) then
                ! Update storm data by interpolation
                call storm_interpolate(storm)
                ! Update current time in storm module (race condition?)
                storm%t = t
            endif
            !$OMP END CRITICAL (INTERP_STORM)
        endif

        ! Set fields
        ! Determine lat/long of each cell in layer,
        !  determine corresponding storm cell indices
        !  (or nearest cell index if out-of-bound)
        !  then get value of corresponding storm data cell.
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            k = get_lat_index(y,storm) ! storm index of latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude
                l = get_lon_index(x,storm) ! storm index of longitude
                ! Set pressure field
                aux(pressure_index,i,j) = storm%p(l,k)
                ! Set velocity components of storm 
                !aux(wind_index,i,j)   = storm%u(l,k) * 1.38 ! TEMPORARY!
                !aux(wind_index+1,i,j) = storm%v(l,k) * 1.38 ! TEMPORARY!
                aux(wind_index,i,j)   = storm%u(l,k)
                aux(wind_index+1,i,j) = storm%v(l,k)
            enddo
        enddo

    end subroutine set_wrf_storm_fields

end module wrf_storm_module

