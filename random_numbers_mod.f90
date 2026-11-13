module pcg64_module
  implicit none
  private
  public :: pcg64_state, pcg64_init, pcg64_next_u64, pcg64_next_f64

  ! 128-bit integer, 64-bit integer, and 64-bit real kinds
  integer, parameter :: i128 = selected_int_kind(38)
  integer, parameter :: i64  = selected_int_kind(18)
  integer, parameter :: r64  = selected_real_kind(15, 307)

  ! PCG64 LCG Multiplier (0x2360ED051FC65DA44385DF649FCCF645)
  integer(i128), parameter :: PCG_MULT = 47054244927992906478985114704771449669_i128
  ! DXSM Permutation Multiplier (0xDA942042E4DD58B5)
  integer(i128), parameter :: DXSM_MULT = 15750247268752932021_i128

  ! Thread-local state object
  type :: pcg64_state
     integer(i128) :: state = 0_i128
     integer(i128) :: inc   = 1_i128
  end type pcg64_state

contains

  ! Initialize generator state and stream ID for thread-safe parallel usage
  subroutine pcg64_init(rng, seed, seq)
    type(pcg64_state), intent(out) :: rng
    integer(i64), intent(in)       :: seed
    integer(i64), intent(in)       :: seq

    ! Increment must be odd: inc = (seq << 1) | 1
    rng%inc = ishft(int(seq, i128), 1) + 1_i128
    rng%state = 0_i128

    ! Step LCG to load initial state
    rng%state = rng%state * PCG_MULT + rng%inc
    rng%state = rng%state + int(seed, i128)
    rng%state = rng%state * PCG_MULT + rng%inc
  end subroutine pcg64_init

  ! Generate next 64-bit pseudo-random unsigned integer
  function pcg64_next_u64(rng) result(res)
    type(pcg64_state), intent(inout) :: rng
    integer(i64) :: res
    integer(i128) :: old_state, hi, xorshifted

    old_state = rng%state
    ! Advance LCG 128-bit state
    rng%state = old_state * PCG_MULT + rng%inc

    ! PCG64-DXSM output permutation
    hi = ishft(old_state, -64)
    xorshifted = ieor(hi, ishft(hi, -32))
    xorshifted = xorshifted * DXSM_MULT
    res = int(ieor(xorshifted, ishft(rng%state, -64)), i64)
  end function pcg64_next_u64

  ! Generate double precision float in uniform [0.0, 1.0) interval
  function pcg64_next_f64(rng) result(res)
    type(pcg64_state), intent(inout) :: rng
    real(r64) :: res
    integer(i64) :: raw

    raw = pcg64_next_u64(rng)
    ! Extract upper 53 bits using logical right shift and scale by 2^(-53)
    res = real(ishft(raw, -11), r64) * (1.0_r64 / 9007199254740992.0_r64)
  end function pcg64_next_f64

end module pcg64_module

MODULE random_numbers_mod
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32, INT64, REAL64
  USE pcg64_module, ONLY: pcg64_state, pcg64_init, pcg64_next_f64
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_random_numbers, pcg_random, get_host_id
  PUBLIC :: ran2, rnor, power, levy, expdev

  ! Precision kind matching precision_kind_mod
  INTEGER, PARAMETER :: dp = REAL64

  ! Rank-local PCG64 state instance
  TYPE(pcg64_state), SAVE :: global_rng

CONTAINS

  ! ===========================================================================
  ! Generate positive 32-bit integer hash from hostname / username
  ! ===========================================================================
  INTEGER FUNCTION get_host_id()
    CHARACTER(LEN=128) :: hostname
    INTEGER :: status, i
    INTEGER(INT32) :: host_hash

    CALL get_environment_variable('HOSTNAME', hostname, status=status)
    IF (status /= 0) THEN
      CALL get_environment_variable('USER', hostname, status=status)
    END IF

    host_hash = 0_INT32
    DO i = 1, LEN_TRIM(hostname)
      host_hash = host_hash + IACHAR(hostname(i:i))
    END DO

    get_host_id = ABS(host_hash) + 1
  END FUNCTION get_host_id

  ! ===========================================================================
  ! Initialize PCG64 using Date/Time, Host ID, and MPI Rank
  ! ===========================================================================
  SUBROUTINE init_random_numbers(rank, custom_seed)
    INTEGER, INTENT(IN) :: rank
    INTEGER(INT64), INTENT(IN), OPTIONAL :: custom_seed

    INTEGER :: values(8)
    INTEGER(INT64) :: base_seed, stream_seq
    INTEGER :: host_id

    host_id = get_host_id()

    IF (PRESENT(custom_seed)) THEN
      base_seed = custom_seed
    ELSE
      ! Retrieve wall-clock time array:
      ! values(3)=day, values(5)=hour, values(6)=min, values(7)=sec, values(8)=ms
      CALL DATE_AND_TIME(VALUES=values)

      ! Unique 64-bit seed combining time and host ID
      base_seed = INT(values(8), INT64) &
                + INT(values(7), INT64) * 1000_INT64 &
                + INT(values(6), INT64) * 60000_INT64 &
                + INT(values(5), INT64) * 3600000_INT64 &
                + INT(values(3), INT64) * 86400000_INT64 &
                + INT(host_id, INT64)   * 1000000000_INT64
    END IF

    ! Unique stream sequence ID per process across different machines
    stream_seq = INT(host_id, INT64) * 1000000_INT64 + INT(rank + 1, INT64)

    CALL pcg64_init(global_rng, base_seed, stream_seq)
  END SUBROUTINE init_random_numbers

  ! ===========================================================================
  ! Uniform double-precision float in [0, 1)
  ! ===========================================================================
  FUNCTION pcg_random() RESULT(random_val)
    REAL(dp) :: random_val
    random_val = REAL(pcg64_next_f64(global_rng), dp)
  END FUNCTION pcg_random

  ! ===========================================================================
  ! Legacy ran2 interface wrapper
  ! ===========================================================================
  FUNCTION ran2(idum) RESULT(random_val)
    INTEGER, INTENT(INOUT), OPTIONAL :: idum
    REAL(dp) :: random_val

    IF (PRESENT(idum)) THEN
      IF (idum < 0) THEN
        CALL pcg64_init(global_rng, INT(ABS(idum), INT64), 1_INT64)
        idum = 1
      END IF
    END IF

    random_val = pcg_random()
  END FUNCTION ran2

  ! ===========================================================================
  ! Normal (Gaussian) distribution
  ! ===========================================================================
  FUNCTION rnor(idum) RESULT(normal_deviate)
    INTEGER, INTENT(INOUT), OPTIONAL :: idum
    REAL(dp) :: normal_deviate

    REAL(dp), PARAMETER :: TWopi = 8.0_dp * ATAN(1.0_dp)
    INTEGER, SAVE :: iv_saved = 0
    REAL(dp), SAVE :: r_saved, t_saved
    REAL(dp) :: u1, u2
    
    IF (PRESENT(idum)) CONTINUE ! Evita il warning di "unused dummy argument"
    
    IF (iv_saved == 0) THEN
      u1 = pcg_random()
      DO WHILE (u1 == 0.0_dp)
        u1 = pcg_random()
      END DO
      u2 = pcg_random()

      r_saved = SQRT(-2.0_dp * LOG(u1))
      t_saved = TWopi * u2

      normal_deviate = r_saved * COS(t_saved)
      iv_saved = 1
    ELSE
      normal_deviate = r_saved * SIN(t_saved)
      iv_saved = 0
    END IF
  END FUNCTION rnor

  ! ===========================================================================
  ! Power-law distribution
  ! ===========================================================================
  FUNCTION power(idum, alpha, t0) RESULT(power_val)
    REAL(dp), INTENT(IN) :: alpha, t0
    INTEGER, INTENT(INOUT), OPTIONAL :: idum
    REAL(dp) :: power_val

    IF (PRESENT(idum)) CONTINUE ! Evita il warning di "unused dummy argument"
    
    power_val = t0 * (pcg_random()**alpha - 1.0_dp)
  END FUNCTION power

  ! ===========================================================================
  ! Stable Lévy distribution
  ! ===========================================================================
  FUNCTION levy(idum, alpha, sigma, beta, mu) RESULT(levy_val)
    REAL(dp), INTENT(IN) :: alpha, sigma, beta, mu
    INTEGER, INTENT(INOUT), OPTIONAL :: idum
    REAL(dp) :: levy_val

    REAL(dp), PARAMETER :: pihalf = 2.0_dp * ATAN(1.0_dp)
    LOGICAL, SAVE :: is_initialized = .FALSE.
    REAL(dp), SAVE :: balbe, salbe, scale, alphainv
    REAL(dp) :: v, w, u1, u2

    IF (.NOT. is_initialized) THEN
      IF (alpha /= 1.0_dp) THEN
        balbe = ATAN(beta * TAN(pihalf * alpha)) / alpha
        salbe = (1.0_dp + (beta * TAN(pihalf * alpha))**2)**(0.5_dp / alpha)
      END IF
      scale = sigma * beta * LOG(sigma) / pihalf + mu
      alphainv = 1.0_dp / alpha
      is_initialized = .TRUE.
    END IF

    IF (PRESENT(idum)) CONTINUE ! Evita il warning di "unused dummy argument"
    
    u1 = pcg_random()
    u2 = pcg_random()
    DO WHILE (u2 == 0.0_dp)
      u2 = pcg_random()
    END DO

    v = pihalf - 2.0_dp * pihalf * u1
    w = -LOG(u2)

    IF (alpha == 1.0_dp) THEN
      levy_val = ((pihalf + beta * v) * TAN(v) - &
                 beta * LOG(w * COS(v) / (pihalf + beta * v))) / pihalf
      levy_val = sigma * levy_val + scale
    ELSE
      levy_val = salbe * SIN(alpha * (v + balbe)) / (COS(v))**alphainv * &
                (COS(v - alpha * (v + balbe)) / w)**((1.0_dp - alpha) * alphainv)
      levy_val = sigma * levy_val + mu
    END IF
  END FUNCTION levy

  ! ===========================================================================
  ! Exponential distribution
  ! ===========================================================================
  FUNCTION expdev(idum) RESULT(expdev_val)
    INTEGER, INTENT(INOUT), OPTIONAL :: idum
    REAL(dp) :: expdev_val
    REAL(dp) :: y

    IF (PRESENT(idum)) CONTINUE ! Evita il warning di "unused dummy argument"
    
    DO
      y = pcg_random()
      IF (y > 0.0_dp) EXIT
    END DO

    expdev_val = -LOG(y)
  END FUNCTION expdev

END MODULE random_numbers_mod


