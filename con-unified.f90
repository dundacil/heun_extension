program con_unified
  use mpi
  use precision_kind_mod, only: dp
  use random_numbers_mod, only: ran2, rnor, init_random_numbers
  implicit none

  ! Parametri della griglia temporale e del campionamento
  integer, parameter :: ndh = 13, npoichk = 16
  integer, parameter :: nnoise = npoichk * 2**ndh

  ! Variabili MPI e controllo I/O
  integer :: ierr, rank, size, root, ios
  integer :: idum, iave, idh, ii, inoise, ipoints, inner, iw, k
  integer :: ninner, nskip, ibase

  ! Parametri configurabili via namelist (con valori di default)
  character(len=16) :: scheme = 'euler'
  integer :: nave = 100000
  integer :: ndh0 = 1
  real(dp) :: d0 = 5.d-2

  ! Variabili fisiche e numeriche
  real(dp) :: dd, dh0, dhm, dh, d1, s3, pi, aa, ratio, xs
  real(dp) :: x0, x1, x2, x3, time, wtot
  real(dp) :: z1, z2, z3, zd2, zd3
  real(dp) :: f0x, f1x, f2x, g0x, g1x, g2x, g3x
  logical :: got, use_z2, use_z3, is_linear_exact

  ! Array per rumore, traiettorie e accumulo residui
  real(dp) :: noise(nnoise), noise2(nnoise), noise3(nnoise)
  real(dp) :: z1a(nnoise), z2a(nnoise), z3a(nnoise)
  real(dp) :: xsave(npoichk, 0:ndh), xexact(npoichk)
  real(dp) :: res_loc(npoichk, ndh), res_glo(npoichk, ndh)

  ! Namelist per la scelta dello schema e dei parametri di simulazione
  namelist /sim_params/ scheme, nave, ndh0, d0

  ! --- Inizializzazione MPI e RNG ---
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  call init_random_numbers(rank)

  root = 0
  idum = -17 - 1000 * rank
  pi = 4.d0 * datan(1.d0)

  ! Lettura della namelist dal file 'input.nml' sul rank 0 e broadcast
  if (rank == root) then
    open(unit=99, file='con_unified.nml', status='old', iostat=ios)
    if (ios == 0) then
      read(99, nml=sim_params)
      close(99)
    else
      print *, "Attenzione: 'con_unified.nml' non trovato, uso valori di default."
    end if
    print *, "Schema selezionato: ", trim(scheme), " | nave =", nave, " | d0 =", d0
  end if

  call MPI_Bcast(scheme, 16, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nave, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndh0, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(d0, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

  ! --- Configurazione specifica per lo schema scelto ---
  use_z2 = .false.
  use_z3 = .false.
  is_linear_exact = .false.
  dh0 = 1.d0 / npoichk

  select case (trim(adjustl(scheme)))
  case ('euler', 'stra', 'heun', 'hest', 'pccc', 'hpcm', 'miln', 'milm', 'rk')
    ! Schemi che richiedono solo gli incrementi di Wiener Z1
  case ('hone', 't3d2', 't3d2d')
    ! Schemi che richiedono gli integrali doppi (Z1, Z2)
    use_z2 = .true.
  case ('hon2')
    ! Honeycutt su SDE lineare con soluzione analitica esatta e intervallo [0, 4]
    use_z2 = .true.
    is_linear_exact = .true.
    dh0 = 4.d0 / npoichk
  case ('cup1', 'cup2', 's2')
    ! Schemi di ordine forte 2.0 che richiedono (Z1, Z2, Z3)
    use_z2 = .true.
    use_z3 = .true.
  case default
    if (rank == root) print *, "ERRORE: scheme non riconosciuto -> ", trim(scheme)
    call MPI_Finalize(ierr)
    stop
  end select

  dd = dsqrt(2.d0 * d0)
  dhm = dh0 / 2**ndh
  d1  = dsqrt(dhm)
  zd2 = dsqrt(dhm) / dsqrt(3.d0)
  zd3 = dhm * dhm / 6.d0

  res_glo = 0.d0
  res_loc = 0.d0

  ! --- Ciclo principale sulle realizzazioni Monte Carlo ---
  do iave = 1, nave

    ! 1. Generazione del rumore alla scala piu' fine (ndh) per tutti gli schemi
    do inoise = 1, nnoise
      noise(inoise) = rnor(idum)
      z1a(inoise)   = d1 * noise(inoise)

      if (use_z2) then
        noise2(inoise) = rnor(idum)
        z2a(inoise)    = dhm * 0.5d0 * (z1a(inoise) + noise2(inoise) * zd2)
      end if

      if (use_z3) then
        noise3(inoise) = rnor(idum)
        select case (trim(adjustl(scheme)))
        case ('cup1')
          z3a(inoise) = dhm * (z1a(inoise)**2 - dhm) / 6.d0 + noise3(inoise) * zd3
        case ('cup2')
          z3a(inoise) = z1a(inoise) * z2a(inoise) - 0.5d0 * dhm * dhm - zd3 * noise3(inoise)
        case ('s2')
          z3a(inoise) = dhm * dhm * ((noise(inoise)**2 - 1.d0) / 6.d0 - &
                        (noise2(inoise)**2 - 1.d0) / 10.d0 + &
                        noise3(inoise) * dsqrt(7.d0) / 30.d0)
        end select
      end if
    end do

    ! 2. Estrazione del punto iniziale dalla distribuzione stazionaria di Stratonovich
    aa = 1.d0 / (2.d0 * d0)
    got = .false.
    do while (.not. got)
      xs = dtan(pi * (ran2(idum) - 0.5d0))
      ratio = 1.d0 / (1.d0 + xs**2)**aa
      if (ratio > ran2(idum)) got = .true.
    end do

    xsave = 0.d0

    ! 3. Integrazione ai diversi passi temporali dh = dh0 / 2**idh
    do idh = ndh0, ndh
      dh = dh0 / 2**idh
      ninner = 2**idh
      nskip = 2**(ndh - idh)
      time = 0.d0
      wtot = 0.d0
      x0 = xs
      s3 = dsqrt(3.d0 * dh)

      do ipoints = 1, npoichk
        do inner = 1, ninner
          ibase = (ipoints - 1) * ninner * nskip + (inner - 1) * nskip

          ! Aggregazione del rumore sulla scala dh (unico blocco per tutti gli schemi)
          z1 = z1a(ibase + 1)
          if (use_z2) z2 = z2a(ibase + 1)
          if (use_z3) z3 = z3a(ibase + 1)

          do iw = 2, nskip
            if (use_z3) then
              z3 = z3 + z3a(ibase + iw) + z2 * z1a(ibase + iw) &
                   + z1 * (dhm * z1a(ibase + iw) - z2a(ibase + iw))
            end if
            if (use_z2) then
              z2 = z2 + z2a(ibase + iw) + dhm * z1
            end if
            z1 = z1 + z1a(ibase + iw)
          end do

          ! Avanzamento di un passo dh secondo lo schema selezionato
          select case (trim(adjustl(scheme)))

          case ('euler')
            x0 = x0 + dh * f0(x0) + g0(x0) * z1

          case ('stra')
            x0 = x0 + dh * f0(x0) + g0(x0) * z1 * (1.d0 + 0.5d0 * z1 * g1(x0))

          case ('heun')
            x1 = x0 + dh * f0(x0) + g0(x0) * z1
            x2 = x0 + dh * f0(x1) + g0(x1) * z1
            x0 = 0.5d0 * (x1 + x2)

          case ('hest')
            x1 = x0 + dh * f0(x0) + g0(x0) * z1 * (1.d0 + 0.5d0 * z1 * g1(x0))
            x2 = x0 + dh * f0(x1) + g0(x1) * z1 * (1.d0 - 0.5d0 * z1 * g1(x1))
            x0 = 0.5d0 * (x1 + x2)

          case ('pccc')
            x1 = x0 + dh * f0(x0) + g0(x0) * z1
            x2 = x0 + dh * f0(x1) + g0(x1) * z1
            do k = 1, 3
              x3 = 0.5d0 * (x1 + x2)
              x2 = x0 + dh * f0(x3) + g0(x3) * z1
            end do
            x0 = 0.5d0 * (x1 + x2)

          case ('hpcm')
            x1 = x0 + dh * f0(x0) + g0(x0) * z1 + 0.5d0 * z1 * z1 * g0(x0) * g1(x0)
            x2 = x0 + dh * f0(x1) + g0(x1) * z1 - 0.5d0 * z1 * z1 * g0(x1) * g1(x1)
            do k = 1, 3
              x3 = 0.5d0 * (x1 + x2)
              x2 = x0 + dh * f0(x3) + g0(x3) * z1 - 0.5d0 * z1 * z1 * g0(x3) * g1(x3)
            end do
            x0 = 0.5d0 * (x1 + x2)

          case ('miln')
            x1 = x0 + dh * (f0(x0) - 0.5d0 * g0(x0) * g1(x0)) + g0(x0) * z1 + 0.5d0 * z1 * z1 * g0(x0) * g1(x0)
            x2 = x0 + dh * (f0(x1) - 0.5d0 * g0(x1) * g1(x1)) + g0(x1) * z1 + 0.5d0 * z1 * z1 * g0(x1) * g1(x1)
            x0 = 0.5d0 * (x1 + x2)

          case ('milm')
            x1 = x0 + dh * (f0(x0) - 0.5d0 * g0(x0) * g1(x0)) + g0(x0) * z1 + 0.5d0 * z1 * z1 * g0(x0) * g1(x0)
            x2 = x0 + dh * (f0(x1) + 0.5d0 * g0(x1) * g1(x1)) + g0(x1) * z1 - 0.5d0 * z1 * z1 * g0(x1) * g1(x1)
            x0 = 0.5d0 * (x1 + x2)

          case ('rk')
            f1x = frk(x0)
            g1x = g0(x0)
            x1  = x0 + dh * f1x + z1 * g1x
            f2x = frk(x1)
            g2x = g0(x0 - (2.d0 / 3.d0) * g1x * (z1 + s3))
            g3x = g0(x0 + (2.d0 / 9.d0) * g1x * (3.d0 * z1 + s3))
            g4x = g0(x0 - (20.d0 / 27.d0) * f1x * dh + (10.d0 / 27.d0) * (g2x - g1x) * z1 - (10.d0 / 27.d0) * g2x * s3)
            x0  = x0 + 0.5d0 * dh * (f1x + f2x) + (z1 / 40.d0) * (37.d0 * g1x + 30.d0 * g3x - 27.d0 * g4x) &
                  + (s3 / 16.d0) * (8.d0 * g1x + g2x - 9.d0 * g3x)

          case ('hone')
            x0 = x0 + dh * frk(x0) + g0(x0) * z1 + 0.5d0 * (z1 * z1 - dh) * g0(x0) * g1(x0) &
                 + frk1(x0) * g0(x0) * z2 + 0.5d0 * (frk(x0) * frk1(x0) + 0.5d0 * g0(x0)**2 * frk2(x0)) * dh * dh &
                 + (z1 * dh - z2) * (frk(x0) * g1(x0) + 0.5d0 * g2(x0) * g0(x0)**2) &
                 + 0.5d0 * g0(x0) * (g0(x0) * g2(x0) + g1(x0)**2) * (z1 * z1 / 3.d0 - dh) * z1

          case ('hon2')
            ! Caso test lineare: f(x) = x, g(x) = dd*x
            g0x = dd * x0
            g1x = dd
            x0 = x0 + dh * x0 + g0x * z1 + 0.5d0 * (z1 * z1 - dh) * g0x * g1x &
                 + 1.d0 * g0x * z2 + 0.5d0 * x0 * dh * dh &
                 + (z1 * dh - z2) * (x0 * g1x) &
                 + 0.5d0 * g0x * (g1x**2) * (z1 * z1 / 3.d0 - dh) * z1

          case ('t3d2')
            g0x = g0(x0); g1x = g1(x0); g2x = g2(x0)
            f0x = f0(x0); f1x = f1(x0)
            x0 = x0 + dh * f0x + g0x * z1 + 0.5d0 * z1 * z1 * g0x * g1x &
                 + z2 * (g0x * f1x - f0x * g1x) + dh * z1 * f0x * g1x &
                 + z1 * z1 * z1 * (g1x * g1x + g2x * g0x) * g0x / 6.d0

          case ('t3d2d')
            g0x = g0(x0); g1x = g1(x0); g2x = g2(x0)
            f0x = f0(x0); f1x = f1(x0)
            x0 = x0 + dh * f0x + g0x * z1 + 0.5d0 * z1 * z1 * g0x * g1x &
                 + z2 * (g0x * f1x - f0x * g1x) + dh * z1 * f0x * g1x &
                 + z1 * z1 * z1 * (g1x * g1x + g2x * g0x) * g0x / 6.d0 &
                 + 0.5d0 * f0x * f1x * dh * dh

          case ('cup1', 'cup2', 's2')
            g0x = g0(x0); g1x = g1(x0); g2x = g2(x0); g3x = g3(x0)
            f0x = f0(x0); f1x = f1(x0); f2x = f2(x0)
            x0 = x0 + dh * f0x + g0x * z1 + 0.5d0 * z1 * z1 * g0x * g1x &
                 + z2 * (g0x * f1x - f0x * g1x) + dh * z1 * f0x * g1x &
                 + z1 * z1 * z1 * (g1x * g1x + g2x * g0x) * g0x / 6.d0 &
                 + 0.5d0 * f0x * f1x * dh * dh &
                 + 0.5d0 * g0x * (g1x * f1x + g0x * f2x) * (z1 * z2 - z3) &
                 + g0x * (f1x * g1x + f0x * g2x) * z3 &
                 + f0x * (g1x * g1x + g0x * g2x) * (0.5d0 * dh * z1 * z1 - 0.5d0 * (z1 * z2 + z3)) &
                 + g0x * (g1x**3 + 4.d0 * g0x * g1x * g2x + g0x * g0x * g3x) * (z1**4) / 24.d0

          end select

          time = time + dh
          if (is_linear_exact) wtot = wtot + z1 / d1
        end do

        if (idh == ndh .and. is_linear_exact) then
          xexact(ipoints) = xs * dexp((1.d0 - 0.5d0 * dd * dd) * time + d1 * dd * wtot)
        end if
        xsave(ipoints, idh) = x0
      end do
    end do

    ! 4. Scrittura traiettorie campionate (prima e ultima realizzazione sul rank 0)
    if (rank == root .and. iave == 1) then
      do ipoints = 1, npoichk
        if (is_linear_exact) then
          write(100, 9) ipoints * dh * 2**ndh, (xsave(ipoints, ii), ii = 1, ndh), xexact(ipoints)
        else
          write(100, 8) ipoints * dh * 2**ndh, (xsave(ipoints, ii), ii = 1, ndh)
        end if
      end do
      close(100)
    end if

    if (rank == root .and. iave == nave) then
      do ipoints = 1, npoichk
        if (is_linear_exact) then
          write(101, 9) ipoints * dh * 2**ndh, (xsave(ipoints, ii), ii = 1, ndh), xexact(ipoints)
        else
          write(101, 8) ipoints * dh * 2**ndh, (xsave(ipoints, ii), ii = 1, ndh)
        end if
      end do
      close(101)
    end if

    ! 5. Accumulo dell'errore forte rispetto al passo fine (o esatto per hon2)
    do ipoints = 1, npoichk
      do idh = ndh0, ndh
        if (is_linear_exact) then
          res_loc(ipoints, idh) = res_loc(ipoints, idh) + dabs(xexact(ipoints) - xsave(ipoints, idh))
        else
          res_loc(ipoints, idh) = res_loc(ipoints, idh) + dabs(xsave(ipoints, ndh) - xsave(ipoints, idh))
        end if
      end do
    end do
  end do

  ! 6. Riduzione MPI e scrittura dei risultati su fort.11 .. fort.26
  call MPI_Reduce(res_loc, res_glo, npoichk * ndh, MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

  if (rank == root) then
    res_glo = res_glo / (size * nave)
    do idh = ndh0, ndh
      dh = dh0 / 2**idh
      do ipoints = 1, npoichk
        write(10 + ipoints, *) dh, res_glo(ipoints, idh)
      end do
    end do
  end if

8 format(15(1x, e10.4))
9 format(16(1x, e10.4))

  call MPI_Finalize(ierr)

contains

  ! --- Definizione delle funzioni di drift, diffusione e derivate ---
  pure real(dp) function f0(x)
    real(dp), intent(in) :: x
    f0 = -x * (1.d0 + x * x)
  end function f0

  pure real(dp) function f1(x)
    real(dp), intent(in) :: x
    f1 = -1.d0 - 3.d0 * x * x
  end function f1

  pure real(dp) function f2(x)
    real(dp), intent(in) :: x
    f2 = -6.d0 * x
  end function f2

  pure real(dp) function g0(x)
    real(dp), intent(in) :: x
    g0 = dd * (x * x + 1.d0)
  end function g0

  pure real(dp) function g1(x)
    real(dp), intent(in) :: x
    g1 = dd * (2.d0 * x)
  end function g1

  pure real(dp) function g2(x)
    real(dp), intent(in) :: x
    g2 = dd * 2.d0
  end function g2

  pure real(dp) function g3(x)
    real(dp), intent(in) :: x
    g3 = 0.d0
  end function g3

  pure real(dp) function frk(x)
    real(dp), intent(in) :: x
    frk = f0(x) + 0.5d0 * g0(x) * g1(x)
  end function frk

  pure real(dp) function frk1(x)
    real(dp), intent(in) :: x
    frk1 = f1(x) + 0.5d0 * (g0(x) * g2(x) + g1(x) * g1(x))
  end function frk1

  pure real(dp) function frk2(x)
    real(dp), intent(in) :: x
    frk2 = f2(x) + 0.5d0 * (g1(x) * g2(x) + g0(x) * g3(x) + 2.d0 * g2(x) * g1(x))
  end function frk2

end program con_unified
