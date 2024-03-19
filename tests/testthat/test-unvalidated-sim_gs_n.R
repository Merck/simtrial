# 2024-02-22: Converted `example("sim_gs_n")` to tests from commit 306de0d
# https://github.com/Merck/simtrial/tree/306de0dbe380fdb1e906a59f34bf3871d3ee5312

# See helper-sim_gs_n.R for helper functions

test_that("Test 1: regular logrank test", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = wlr,
    cuts = test_cutting(),
    seed = 2024,
    weight = fh(rho = 0, gamma = 0)
  )
  expected <- data.frame(
    rho = numeric(9),
    gamma = numeric(9),
    z = c(
      -3.7486049782713247, -4.53034007934394, -4.316452743033609,
      -3.4771440155825752, -3.8631501353780324, -3.2777779731288317,
      -3.075862925191481, -3.619345457605645, -4.2225917786532925
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 2: weighted logrank test by FH(0, 0.5)", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = wlr,
    cuts = test_cutting(),
    seed = 2024,
    weight = fh(rho = 0, gamma = 0.5)
  )
  expected <- data.frame(
    rho = numeric(9),
    gamma = rep(0.5, 9L),
    z = c(
      -4.149161171743935, -4.778107819550277, -4.2607297587160256,
      -3.605092910242299, -3.945081123231263, -2.919179640988388,
      -3.1432278107909206, -3.640458610667732, -4.243289152457
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 3: weighted logrank test by MB(3)", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = wlr,
    cuts = test_cutting(),
    seed = 2024,
    weight = mb(delay = 3)
  )
  expected <- data.frame(
    z = c(
      -3.797133894694147, -4.581330588107247, -4.3496437937060906,
      -3.5011312494121394, -3.886541892591609, -3.2792862684447983,
      -3.114079263266195, -3.6587146250230145, -4.2632793831797855
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 4: weighted logrank test by early zero (6)", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = wlr,
    cuts = test_cutting(),
    seed = 2024,
    weight = early_zero(6)
  )
  expected <- data.frame(
    z = c(
      -4.552617167258777, -5.188572984743822, -4.686073828268738,
      -3.185533497487861, -3.5975030245947046, -2.786930008687834,
      -2.3673440974318556, -3.0630537456426414, -3.7816194091003705
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 5: RMST", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = rmst,
    cuts = test_cutting(),
    seed = 2024,
    tau = 20
  )
  expected <- data.frame(
    rmst_arm1 = c(
      12.466259284156251, 12.444204897288326, 12.425100778728808,
      12.392111715564337, 12.496963791557544, 12.479119007501355, 12.62769367846186,
      12.737915554271744, 12.740241766667666
    ),
    rmst_arm0 = c(
      9.585107633112955, 9.591073977478539, 9.590592780789704, 9.824721964671674,
      10.097271436421035, 10.110783864663125, 10.340195893022198,
      10.289798076615766, 10.261299533752227
    ),
    rmst_diff = c(
      2.8811516510432966, 2.8531309198097876, 2.834507997939104, 2.567389750892662,
      2.3996923551365086, 2.36833514283823, 2.287497785439662, 2.4481174776559786,
      2.478942232915438
    ),
    z = c(
      3.7899815357169184, 3.991862864282945, 3.980100861311682, 3.474868814723485,
      3.2950209410683957, 3.2541151987300845, 2.9805344295194454,
      3.3009521580248022, 3.3504301652133
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 6: Milestone", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = milestone,
    cuts = test_cutting(),
    seed = 2024,
    ms_time = 10
  )
  expected <- data.frame(
    method = rep("milestone", 9L),
    z = c(
      9.252619142383594, 12.078380683791904, 12.078380683791904, 5.565741269919053,
      5.457930240636103, 5.457930240636103, 9.051772787302813, 9.054982526543846,
      9.054982526543846
    ),
    ms_time = rep(10, 9L),
    surv_ctrl = c(
      0.40800409626773176, 0.40972689075630214, 0.40972689075630214,
      0.4718268722892688, 0.46670065754089335, 0.46670065754089335,
      0.46149611243704863, 0.46199999999999974, 0.46199999999999974
    ),
    surv_exp = c(
      0.568975019886668, 0.5849999999999997, 0.5849999999999997, 0.5922853919588814,
      0.5840900715499292, 0.5840900715499292, 0.6150543366195163,
      0.6139773404060171, 0.6139773404060171
    ),
    surv_diff = c(
      0.16097092361893622, 0.1752731092436976, 0.1752731092436976,
      0.12045851966961263, 0.11738941400903585, 0.11738941400903585,
      0.15355822418246762, 0.1519773404060174, 0.1519773404060174
    ),
    std_err_ctrl = c(
      0.03693587681297664, 0.034952703615152854, 0.034952703615152854,
      0.03614098127448581, 0.035432630739150366, 0.035432630739150366,
      0.035815727559287504, 0.03540131462139614, 0.03540131462139614
    ),
    std_err_exp = c(
      0.03662189834863626, 0.03484070894801079, 0.03484070894801079,
      0.035312669921649095, 0.034912158581439694, 0.034912158581439694,
      0.03505127094114008, 0.034738243333119145, 0.034738243333119145
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("Test 7: MaxCombo (WLR-FH(0,0) + WLR-FH(0, 0.5))", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = maxcombo,
    cuts = test_cutting(),
    seed = 2024,
    rho = c(0, 0),
    gamma = c(0, 0.5)
  )
  expected <- data.frame(
    p_value = c(
      2.6155386454673746e-05, 1.4330486162172917e-06, 1.247801863046849e-05,
      0.0002358380298724816, 6.130077643518028e-05, 0.0007667834024346343,
      0.001216230102102256, 0.00020471863687732128, 1.7249355113824194e-05
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("sim_gs_n() accepts different tests per cutting", {
  wlr_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  wlr_cut2 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
  wlr_cut3 <- create_test(wlr, weight = fh(rho = 0.5, gamma = 0))

  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    tests = list(wlr_cut1, wlr_cut2, wlr_cut3),
    cuts = test_cutting(),
    seed = 2024
  )
  expected <- data.frame(
    rho = rep(c(0, 0, 0.5), 3),
    gamma = rep(c(0, 0.5, 0), 3),
    z = c(
      -3.7486049782713247, -4.778107819550277, -4.189693884801371,
      -3.4771440155825752, -3.945081123231263, -3.438138809871842,
      -3.075862925191481, -3.640458610667732, -3.9489173860678495
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  expect_equal(observed, expected)
})

test_that("sim_gs_n() requires a test for each cutting", {
  wlr_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  wlr_cut2 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))

  expect_error(
    sim_gs_n(
      n_sim = 3,
      sample_size = 400,
      enroll_rate = test_enroll_rate(),
      fail_rate = test_fail_rate(),
      tests = list(wlr_cut1, wlr_cut2),
      cuts = test_cutting(),
      seed = 2024
    ),
    "If you want to run different tests at each cutting"
  )
})
