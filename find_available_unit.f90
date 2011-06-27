integer function find_available_unit()

  integer, parameter :: min_unit = 7
  integer, parameter :: max_unit = 99

  logical :: exists, already_used

  continue

  do find_available_unit = min_unit, max_unit
    inquire(unit = find_available_unit, exist = exists, opened = already_used)
    if (exists .and. .not. already_used) return
  end do

end function find_available_unit
