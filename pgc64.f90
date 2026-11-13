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
