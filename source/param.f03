module param

implicit none

public


! ================================================================
! particle pusher
! ================================================================
integer, parameter :: p_push2_std = 0, p_push2_robust = 1, p_push2_robust_subcyc = 2, &
                      p_push2_clamp = 3, p_push2_std_pgc = 4, p_push2_robust_pgc = 5


end module param
