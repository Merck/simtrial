# 2024-02-22: Converted `example("sim_gs_n")` to tests from commit 306de0d
# https://github.com/Merck/simtrial/tree/306de0dbe380fdb1e906a59f34bf3871d3ee5312

skip_if_not_installed("gsDesign2")

# parameters for enrollment
enroll_rampup_duration <- 4 # duration for enrollment ramp up
enroll_duration <- 16 # total enrollment duration
enroll_rate <- gsDesign2::define_enroll_rate(
  duration = c(
    enroll_rampup_duration,
    enroll_duration - enroll_rampup_duration
    ),
  rate = c(10, 30)
)

# parameters for treatment effect
delay_effect_duration <- 3 # delay treatment effect in months
median_ctrl <- 9 # survival median of the control arm
median_exp <- c(9, 14) # survival median of the experimental arm
dropout_rate <- 0.001
fail_rate <- gsDesign2::define_fail_rate(
    duration = c(delay_effect_duration, 100),
    fail_rate = log(2) / median_ctrl,
    hr = median_ctrl / median_exp,
    dropout_rate = dropout_rate
)

# other related parameters
alpha <- 0.025 # type I error
beta <- 0.1 # type II error
ratio <- 1 # randomization ratio (exp:ctrl)
# Define cuttings of 2 IAs and 1 FA
# IA1
# The 1st interim analysis will occur at the later of the following 3 conditions:
# - At least 20 months have passed since the start of the study
# - At least 100 events have occurred
# - At least 20 months have elapsed after enrolling 200/400 subjects, with a
#   minimum of 20 months follow-up
# However, if events accumulation is slow, we will wait for a maximum of 24 months.
ia1_cut <- create_cut(
  planned_calendar_time = 20,
  target_event_overall = 100,
  max_extension_for_target_event = 24,
  min_n_overall = 200,
  min_followup = 20
)

# IA2
# The 2nd interim analysis will occur at the later of the following 3 conditions:
# - At least 32 months have passed since the start of the study
# - At least 250 events have occurred
# - At least 10 months after IA1
# However, if events accumulation is slow, we will wait for a maximum of 34 months.
ia2_cut <- create_cut(
  planned_calendar_time = 32,
  target_event_overall = 200,
  max_extension_for_target_event = 34,
  min_time_after_previous_analysis = 10
)

# FA
# The final analysis will occur at the later of the following 2 conditions:
# - At least 45 months have passed since the start of the study
# - At least 300 events have occurred
fa_cut <- create_cut(
  planned_calendar_time = 45,
  target_event_overall = 350
)

cut <- list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)

test_that("regular logrank test", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = wlr,
    cut = cut,
    weight = fh(rho = 0, gamma = 0)
  ) |>
    dplyr::select(-c(info, info0))
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep("FH(rho=0, gamma=0)", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -14.8967761757316, -22.6791676993923, -29.4104630799085,
      -13.3530329578349, -17.9272135845997, -25.9789692783839,
      -15.7016927028295, -22.7155477802184, -34.9030500614426
    ),
    se = c(
      7.79986198971421, 8.70892718893558, 9.29168172400568,
      7.65165409688827, 8.77164253357172, 9.41875140383468,
      7.44263362582829, 8.42150520256931, 9.20559144909002
    ),
    z = -c(
      -1.90987689210094, -2.60412875287388, -3.16524650257064,
      -1.74511717188905, -2.04376928448542, -2.75821795952773,
      -2.10969577332675, -2.6973263370173,  -3.79150544041283
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  expect_equal(observed, expected, ignore_attr = TRUE)
})

test_that("regular logrank test parallel", {
  set.seed(2024)
  plan("multisession", workers = 2)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = wlr,
    cut = cut,
    weight = fh(rho = 0, gamma = 0)
  ) |>
    dplyr::select(-c(info, info0))
  plan("sequential")
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep("FH(rho=0, gamma=0)", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -14.8967761757316, -22.6791676993923, -29.4104630799085,
      -13.3530329578349, -17.9272135845997, -25.9789692783839,
      -15.7016927028295, -22.7155477802184, -34.9030500614426
    ),
    se = c(
      7.79986198971421, 8.70892718893558, 9.29168172400568,
      7.65165409688827, 8.77164253357172, 9.41875140383468,
      7.44263362582829, 8.42150520256931, 9.20559144909002
    ),
    z = -c(
      -1.90987689210094, -2.60412875287388, -3.16524650257064,
      -1.74511717188905, -2.04376928448542, -2.75821795952773,
      -2.10969577332675, -2.6973263370173,  -3.79150544041283
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  expect_equal(observed, expected, ignore_attr = TRUE)
})

test_that("weighted logrank test by FH(0, 0.5)", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = wlr,
    cut = cut,
    weight = fh(rho = 0, gamma = 0.5)
  ) |>
    dplyr::select(-c(info, info0))
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep("FH(rho=0, gamma=0.5)", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -11.5775033174745, -18.1616545384682, -24.0266481696649,
      -8.44736378922731, -12.0875497249867, -19.3796342804342,
      -9.95948396485636, -15.7876875029804, -26.5907291815672
    ),
    se = c(
      4.34496240294236, 5.37217800719008, 6.12650200890773,
      4.16899365283586, 5.47959154932447, 6.30461298220442,
      3.99084589541075, 5.03739766754931, 6.02201456006756
    ),
    z = -c(
      -2.66458078201881, -3.3806874072603,  -3.92175635211266,
      -2.02623570402444, -2.20592166700397, -3.0738816696815,
      -2.49558219632314, -3.13409592510117, -4.41558699606811
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  expect_equal(observed, expected, ignore_attr = TRUE)
})

test_that("weighted logrank test by MB(3)", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = wlr,
    cut = cut,
    weight = mb(delay = 3)
  ) |>
    dplyr::select(-c(info, info0))
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep("MB(delay = 3, max_weight = Inf)", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -19.5592474225953, -29.5687542054386, -38.2263688427585,
      -16.3881266934583, -22.1124941540687, -32.188879862274,
      -19.8774758099929, -28.6176512547267, -43.8048638497507
    ),
    se = c(
      9.61093979157535, 10.8257248464033, 11.5995352667424,
      9.23731341784668, 10.6834095617449, 11.5139237244706,
      8.93809056261512, 10.198207959654,  11.2011384564222
    ),
    z = -c(
      -2.03510248183433, -2.73134174616145, -3.29550865303709,
      -1.7741226211721,  -2.0697974767577,  -2.79564817629134,
      -2.22390628856832, -2.80614509607408, -3.91075103840314
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  expect_equal(observed, expected, ignore_attr = TRUE)
})

test_that("weighted logrank test by early zero (6)", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = wlr,
    cut = cut,
    weight = early_zero(6)
  ) |>
    dplyr::select(-c(info, info0))
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep("Xu 2017 with first 6 months of 0 weights", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -13.1571999820816, -20.9395915057423, -27.6708868862585,
      -12.9176255501546, -17.4918061769194, -25.5435618707036,
      -5.11317406520125, -12.1270291425901, -24.3145314238143
    ),
    se = c(
      4.90920032621625, 6.25362403463101, 7.04256700674935, 5.10417050246453,
      6.66681774436398, 7.49784129646925, 4.65652183734504, 6.10017624419681,
      7.14376051256763
    ),
    z = -c(
      -2.68011063060906, -3.34839309011608, -3.92909103452473,
      -2.5307982058823,  -2.62371146889484, -3.4067888156998,
      -1.09806723640677, -1.98798012666056, -3.40360393955524
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  expect_equal(observed, expected, ignore_attr = TRUE)
})

test_that("RMST", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = rmst,
    cut = cut,
    tau = 20
  )
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("RMST", 9L),
    parameter = rep(20, 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      1.29434319527997, 1.03579407229687, 1.05740233262724, 1.27673666400642,
      1.34542683554887, 1.355651597816,   1.52706388710948, 1.33389386467127,
      1.34645188089274
    ),
    se = c(
      0.772544846741024, 0.738409622827227, 0.738497623667433,
      0.771657256733873, 0.738351348672298, 0.738159118543527,
      0.782892168148695, 0.739072978624375, 0.736694772438538
    ),
    z = c(
      1.67542790653533, 1.40273642200249, 1.4318290252257,  1.65453853101876,
      1.82220407393879, 1.83653031407491, 1.95054178498237, 1.80482023189919,
      1.82769300294591
    )
  )
  expect_equal(observed, expected)
})

test_that("Milestone", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = milestone,
    cut = cut,
    ms_time = 10,
    test_type = "naive"
  )
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("milestone", 9L),
    parameter = rep(10, 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      0.0570882489260, 0.0550000000000, 0.0550000000000,
      0.0982363112857, 0.0995332208091, 0.0995332208091,
      0.0314033376240, 0.0475912427493, 0.0475912427493
    ),
    se = c(
      0.0503182801363, 0.0499186838769, 0.0499186838769,
      0.0500025291675, 0.0498435166464, 0.0498435166464,
      0.0512870375551, 0.0498652399931, 0.0498652399931
    ),
    z = c(
      1.1345429289596, 1.1017918688652, 1.1017918688652,
      1.9646268483073, 1.9969140924625, 1.9969140924625,
      0.6123055477761, 0.9543971463075, 0.9543971463075
    )
  )
  expect_equal(observed, expected)
})

test_that("WLR with fh(0, 0.5) test at IA1, WLR with mb(6, Inf) at IA2, and milestone test at FA", {
  ia1_test <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
  ia2_test <- create_test(wlr, weight = mb(delay = 6, w_max = Inf))
  fa_test <- create_test(milestone, ms_time = 10, test_type = "naive")

  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = list(ia1 = ia1_test, ia2 = ia2_test, fa = fa_test),
    cut = cut
  )

  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep(c("WLR", "WLR", "milestone"), 3),
    parameter = rep(c("FH(rho=0, gamma=0.5)", "MB(delay = 6, max_weight = Inf)", "10"), 3),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -11.5775033174745,  -36.9093856541259,  0.05500000,
      -8.44736378922731,  -25.8424460996795,  0.09953322,
      -9.95948396485636,  -32.6844032640339,  0.04759124
    ),
    se = c(
      4.34496240294236,   12.4451486506265,   0.04991868,
      4.16899365283586,   12.0591010341806,   0.04984352,
      3.99084589541075,   11.6181044782549,   0.04986524
    ),
    z = c(
      2.66458078201881, 2.96576494908061, 1.1017919,
      2.02623570402444, 2.14298279999738, 1.9969141,
      2.49558219632314, 2.81323027566226, 0.9543971
    ),
    info = c(
      18.234728155604408, 154.73996276625465, NA, 17.323744020480003,
      146.078846032117, NA, 15.524151518925096, 134.8967385215799, NA
    ),
    info0 = c(
      18.97777882184863, 157.11581663819095, NA, 17.488976711530597,
      146.64084968485673, NA, 15.969731304813632, 136.27941618722392, NA
    )
  )
  expect_equal(observed, expected)
})

test_that("MaxCombo (WLR-FH(0,0) + WLR-FH(0, 0.5))", {
  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = maxcombo,
    cut = cut,
    rho = c(0, 0),
    gamma = c(0, 0.5)
  )
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("MaxCombo", 9L),
    parameter = rep("FH(0, 0) + FH(0, 0.5)", 9L),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    z = I(list(
      c(-1.90987689210094, -2.66458078201881),
      c(-2.60412875287388, -3.3806874072603),
      c(-3.16524650257064, -3.92175635211266),
      c(-1.74511717188905, -2.02623570402444),
      c(-2.04376928448542, -2.20592166700397),
      c(-2.75821795952773, -3.0738816696815),
      c(-2.10969577332675, -2.49558219632314),
      c(-2.6973263370173,  -3.13409592510117),
      c(-3.79150544041283, -4.41558699606811)
    )),
    p_value = c(
      0.00540201317191025,  0.000532495483968387, 6.72595178197177e-05,
      0.0283259325313511,   0.0184175606190634,   0.00152059320052655,
      0.00872501551218574,  0.00124820253291602,  7.95668913800007e-06
    )
  )
  expect_equal(observed, expected)
})

test_that("sim_gs_n() accepts different tests per cutting", {
  wlr_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  wlr_cut2 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
  wlr_cut3 <- create_test(wlr, weight = fh(rho = 0.5, gamma = 0))

  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = list(wlr_cut1, wlr_cut2, wlr_cut3),
    cut = cut
  )

  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep("WLR", 9L),
    parameter = rep(c("FH(rho=0, gamma=0)", "FH(rho=0, gamma=0.5)", "FH(rho=0.5, gamma=0)"), 3),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -14.8967761757316, -18.1616545384682, -16.9584344921688,
      -13.3530329578349, -12.0875497249867, -16.4421705167356,
      -15.7016927028295, -15.7876875029804, -21.4586952735303
    ),
    se = c(
      7.79986198971421, 5.37217800719008, 6.98579432813984,
      7.65165409688827, 5.47959154932447, 6.99748048599332,
      7.44263362582829, 5.03739766754931, 6.96263273237168
    ),
    z = -c(
      -1.90987689210094, -3.3806874072603,  -2.42755994459466,
      -1.74511717188905, -2.20592166700397, -2.34972724106162,
      -2.10969577332675, -3.13409592510117, -3.08198006391483
    ),
    info = c(
      60.30737704918032, 28.86098470336026, 49.40288707640943,
      58.36595744680852, 30.32561220505373, 49.14163369066358,
      54.73873873873875, 25.34221836733899, 48.86012654643182
    ),
    info0 = c(
      61.00000000000000, 29.50363707013600, 49.54649841796000,
      58.75000000000000, 30.38923168829244, 49.34042361519703,
      55.50000000000000, 25.73690356687072, 49.15868877945109
    )
  )
  class(expected) <- c("simtrial_gs_wlr", class(expected))
  attr(expected, "method") <- unique(observed$parameter[observed$sim_id == 1])
  expect_equal(observed, expected)
})

test_that("sim_gs_n() requires a test for each cutting", {
  skip_if_not_installed("gsDesign2")

  wlr_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  wlr_cut2 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))

  set.seed(2024)
  expect_error(
    sim_gs_n(
      n_sim = 3,
      sample_size = 400,
      enroll_rate = enroll_rate,
      fail_rate = fail_rate,
      test = list(wlr_cut1, wlr_cut2),
      cut = cut
    ),
    "If you want to run different tests at each cutting"
  )
})

test_that("sim_gs_n() can combine wlr(), rmst(), and milestone() tests", {
  test_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  test_cut2 <- create_test(rmst, tau = 20)
  test_cut3 <- create_test(milestone, ms_time = 10, test_type = "naive")

  set.seed(2024)
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    test = list(test_cut1, test_cut2, test_cut3),
    cut = cut
  )
  expected <- data.frame(
    sim_id = rep(1:3, each = 3L),
    method = rep(c("WLR", "RMST", "milestone"), 3),
    parameter = rep(c("FH(rho=0, gamma=0)", "20", "10"), 3),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 45, 24, 32, 45),
    n = rep(400L, 9L),
    event = c(244, 307, 362, 235, 310, 361, 222, 286, 351),
    estimate = c(
      -14.8967761757316,  1.03579407229687,   0.05500000,
      -13.3530329578349,  1.34542683554887,   0.09953322,
      -15.7016927028295,  1.33389386467127,   0.04759124
    ),
    se = c(
      7.79986198971421,   0.738409622827227,  0.04991868,
      7.65165409688827,   0.738351348672298,  0.04984352,
      7.44263362582829,   0.739072978624375,  0.04986524
    ),
    z = c(
      1.90987689210094, 1.40273642200249,  1.1017919,
      1.74511717188905, 1.82220407393879,  1.9969141,
      2.10969577332675, 1.80482023189919,  0.9543971
    ),
    info = c(
      60.307377049180324, NA, NA, 58.365957446808515, NA,
      NA, 54.738738738738746, NA, NA
    ),
    info0 = c(61, NA, NA, 58.75, NA, NA, 55.5, NA, NA)
  )
  expect_equal(observed, expected)
})

test_that("convert_list_to_df_w_list_cols() is robust to diverse input", {
  x <- list(
    num_single = 0.5,
    num_multi = seq(0, 1, by = 0.1),
    chr_single = "a",
    chr_multi = letters,
    int_1 = 1L,
    int_multi = 1L:10L
  )
  observed <- convert_list_to_df_w_list_cols(x)
  expected <- data.frame(
    num_single = 0.5,
    num_multi = I(list(
      c(
        0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6000000000000001,
        0.7000000000000001, 0.8, 0.9, 1
      )
    )),
    chr_single = "a",
    chr_multi = I(list(
      c(
        "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o",
        "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
      )
    )),
    int_1 = 1L,
    int_multi = I(list(1:10))
  )
  expect_equal(observed, expected)
})

test_that("create_cut() can accept variables as arguments", {
  # https://github.com/Merck/simtrial/issues/260
  ratio <- 1

  enroll_rate <- gsDesign2::define_enroll_rate(duration = c(2, 2, 8),
                                               rate = c(1, 2, 3))

  fail_rate <- gsDesign2::define_fail_rate(duration = c(4, Inf),
                                           fail_rate = log(2) / 12,
                                           hr = c(1, .6),
                                           dropout_rate = .001)

  alpha <- 0.025
  beta <- 0.1

  upper <- gsDesign2::gs_spending_bound
  upar <- list(sf = gsDesign::sfLDOF, total_spend = alpha)
  test_upper <- rep(TRUE, 2)

  lower <- gsDesign2::gs_spending_bound
  lpar <- list(sf = gsDesign::sfLDOF, total_spend = beta)
  test_lower <- c(TRUE, FALSE)
  binding <- FALSE

  info_frac = NULL
  analysis_time = c(24, 36)

  x <- gsDesign2::gs_design_ahr(enroll_rate = enroll_rate, fail_rate = fail_rate,
                                alpha = alpha, beta = beta, ratio = ratio,
                                info_frac = info_frac, analysis_time = analysis_time,
                                upper = upper, upar = upar, test_upper = test_upper,
                                lower = lower, lpar = lpar, test_lower = test_lower,
                                binding = binding) |> gsDesign2::to_integer()

  ia_cut <- simtrial::create_cut(planned_calendar_time = x$analysis$time[1])
  fa_cut <- simtrial::create_cut(planned_calendar_time = x$analysis$time[2])

  # Must run the parallel version first for 2 reasons:
  #
  # 1. In order to reproduce the bug that inspired this test, the cutting
  # functions must not have been evaluated prior to running them in parallel
  #
  # 2. In order to avoid an R CMD check NOTE about detritus in the temp
  # directory, need to run `future::plan("sequential")` to shut down parallel
  # cluster
  # https://future.futureverse.org/articles/future-7-for-package-developers.html#making-sure-to-stop-parallel-workers
  future::plan("multisession", workers = 2)
  set.seed(1)
  results_parallel <- simtrial::sim_gs_n(
    n_sim = 1e2,
    sample_size = x$analysis$n[2],
    enroll_rate = x$enroll_rate,
    fail_rate = x$fail_rate,
    test = simtrial::wlr,
    cut = list(ia = ia_cut, fa = fa_cut),
    weight = simtrial::fh(rho = 0, gamma = 0))

  future::plan("sequential")
  set.seed(1)
  results_sequential <- simtrial::sim_gs_n(
    n_sim = 1e2,
    sample_size = x$analysis$n[2],
    enroll_rate = x$enroll_rate,
    fail_rate = x$fail_rate,
    test = simtrial::wlr,
    cut = list(ia = ia_cut, fa = fa_cut),
    weight = simtrial::fh(rho = 0, gamma = 0))

  expect_equal(results_parallel, results_sequential)
})

test_that("Updating bounds changes the simulation results", {
  x <- gsDesign2::gs_design_ahr(analysis_time = 1:3*12) |>
    gsDesign2::to_integer()

  # No boundary updates
  set.seed(1)
  run1 <- sim_gs_n(
    n_sim = 1,
    sample_size = max(x$analysis$n),
    enroll_rate = x$enroll_rate,
    fail_rate = x$fail_rate,
    test = wlr,
    cut = list(ia1 = create_cut(planned_calendar_time = x$analysis$time[1]),
               ia2 = create_cut(planned_calendar_time = x$analysis$time[2]),
               fa = create_cut(planned_calendar_time = x$analysis$time[3])),
    weight = fh(rho = 0, gamma = 0)
  )

  # With boundary updates
  set.seed(1)
  run2 <- sim_gs_n(
    n_sim = 1,
    sample_size = max(x$analysis$n),
    enroll_rate = x$enroll_rate,
    fail_rate = x$fail_rate,
    test = wlr,
    cut = list(ia1 = create_cut(planned_calendar_time = x$analysis$time[1]),
               ia2 = create_cut(planned_calendar_time = x$analysis$time[2]),
               fa = create_cut(planned_calendar_time = x$analysis$time[3])),
    weight = fh(rho = 0, gamma = 0),
    original_design = x
  )

  expect_equal(run1, run2[, colnames(run1)], ignore_attr = TRUE)

  expected <- data.frame(
    planed_upper_bound = c(3.870248012128966, 2.3566552618098884, 2.009757742407378),
    planed_lower_bound = c(-1.705270817327003, 0.9601286375623664, 2.004752252887608),
    updated_upper_bound = c(3.870248012128966, 2.3867954048423474, 2.0074221828251764),
    updated_lower_bound = c(-1.6671962217546439, 0.9631736579151768, 2.1126105535696467)
  )
  observed <- run2[, c("planned_upper_bound", "planned_lower_bound",
                       "updated_upper_bound", "updated_lower_bound")]
  expect_equal(observed, expected, ignore_attr = TRUE)
})
