PROGRAM  Main
  use omp_lib
  ! $OMP PARALLEL

  print *,OMP_GET_NUM_PROCS()

  ! $OMP END PARALLEL

END PROGRAM Main
