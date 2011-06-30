! Will hold new q1d_nozzle code

  call read_nml

  call read_grid

  call allocate_soln(cells)

  call initialize_soln
    if ( restart ) call read_restart inside initialize_soln

  if ( explicit ) then
    call explicit_solve
  else
    call implicit_solve
  end if


  subroutine explicit_solve()

    implicit none

    continue

    do n = 1, iterations
      call set_time_step( x_f, prim_cc, dt )

      do rk = 1, rk_order
        call create_residual( cells, faces, n, prim_cc,                        &
                              area_f, area_cc, dadx_cc, residual )

        if ( rk == 1 ) cons_cc_0 = cons_cc

        do cell = 2, cells+1
          cons_cc(:,cell) = cons_cc_0(:,cell)                                  &
                          + dt*residual(:,cell)/real(1+rk_order-rk,dp)
          prim_cc(:,cell) = conserved_to_primitive_1D(cons_cc(:,cell))
          prim_cc(:,cell) = floor_primitive_vars_1D(prim_cc(:,cell))
        end do

        call set_inflow(prim_cc(:,1), prim_cc(:,2), prim_cc(:,3))
        call set_outflow(prim_cc(:,cells), prim_cc(:,cells+1),                 &
                         prim_cc(:,cells+2))

        do cell = 1, cells+2
          cons_cc(:,cell) = primitive_to_conserved_1D(prim_cc(:,cell))
        end do
      end do

      call check_convergence(residual, convergence_flag)

      if ( convergence_flag ) break

    end do

  end subroutine explicit_solve
