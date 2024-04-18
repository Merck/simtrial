# 2024-02-22: Converted `example("sim_gs_n")` to tests from commit 306de0d
# https://github.com/Merck/simtrial/tree/306de0dbe380fdb1e906a59f34bf3871d3ee5312

# See helper-sim_gs_n.R for helper functions

test_that("Test 1: regular logrank test", {
  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    test = wlr,
    cut = test_cutting(),
    seed = 2024,
    weight = fh(rho = 0, gamma = 0)
  )
  expected <- data.frame(
    method = rep("WLR", 9L),
    parameter = rep("FH(rho=0, gamma=0)", 9L),
    estimation = c(
      -28.194825408790173, -38.32580538077858, -39.49229553865729,
      -26.84871111584948, -32.548237296118835, -30.06631062297029,
      -23.063020152157016, -30.16329862679027, -38.75506042018556
    ),
    se = c(
      7.521418120132856, 8.459807588292295, 9.149247748025164, 7.7214837796562,
      8.425309955739992, 9.17277218574699, 7.498065002594769, 8.333909813280053,
      9.178026778744327
    ),
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
    test = wlr,
    cut = test_cutting(),
    seed = 2024,
    weight = fh(rho = 0, gamma = 0.5)
  )
  expected <- data.frame(
    method = rep("WLR", 9L),
    parameter = rep("FH(rho=0, gamma=0.5)", 9L),
    estimation = c(
      -16.934217242190208, -24.448179085866208, -25.51076208491462,
      -15.500239367897708, -19.967690764549445, -17.5390556887186,
      -12.664624037110103, -18.05051250570102, -25.59217169575864
    ),
    se = c(
      4.081359229309616, 5.116707284384241, 5.987416130471112, 4.299539499761723,
      5.061414490811455, 6.008213897648367, 4.02917790229253, 4.958307300296486,
      6.031210878226377
    ),
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
    test = wlr,
    cut = test_cutting(),
    seed = 2024,
    weight = mb(delay = 3)
  )
  expected <- data.frame(
    method = rep("WLR", 9L),
    parameter = rep("MB(delay = 3, max_weight = Inf)", 9L),
    estimate = c(
      -34.1924345680359, -46.745479781695614, -48.190848712798775,
      -32.192766832733724, -39.186116293163025, -36.14077883676622,
      -27.895635278073794, -36.66945854377384, -47.28630966948101
    ),
    se = c(
      9.004800861990685, 10.203472306286518, 11.079263268070516, 9.194961439431635,
      10.082514836095871, 11.020928299103936, 8.957907914269184, 10.022497598741575,
      11.091534337637594
    ),
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
    test = wlr,
    cut = test_cutting(),
    seed = 2024,
    weight = early_zero(6)
  )
  expected <- data.frame(
    method = rep("WLR", 9L),
    parameter = rep("Xu 2017 with first 6 months of 0 weights", 9L),
    estimate = c(
      -21.998993527998245, -32.129973499986654, -33.29646365786535,
      -17.199406900467533, -22.89893308073689, -20.417006407588342,
      -11.776058510868394, -18.876336985501645, -27.468098778896934
    ),
    se = c(
      4.8321641639910355, 6.192448982496682, 7.105407400328277, 5.399223368403168,
      6.36522969520413, 7.325984629661098, 4.974375513742725, 6.162587585135993,
      7.2635809708390155
    ),
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
    test = rmst,
    cut = test_cutting(),
    seed = 2024,
    tau = 20
  )
  expected <- data.frame(
    method = rep("RMST", 9L),
    parameter = rep(20, 9L),
    estimation = c(
      2.8811516510432966, 2.8531309198097876, 2.834507997939104, 2.567389750892662,
      2.3996923551365086, 2.36833514283823, 2.287497785439662, 2.4481174776559786,
      2.478942232915438
    ),
    se = c(
      0.7602020283980866, 0.7147367073498636, 0.7121698913441514,
      0.7388450867596175, 0.7282783320819864, 0.7277969580678861,
      0.7674790677752639, 0.7416397937499536, 0.7398877489385366
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
    test = milestone,
    cut = test_cutting(),
    seed = 2024,
    ms_time = 10
  )
  expected <- data.frame(
    method = rep("milestone", 9L),
    parameter = rep(10, 9L),
    estimation = c(
      0.16097092361893622, 0.1752731092436976, 0.1752731092436976,
      0.12045851966961263, 0.11738941400903585, 0.11738941400903585,
      0.15355822418246762, 0.1519773404060174, 0.1519773404060174
    ),
    se = I(list(
      c(0.03693587681297664, 0.03662189834863626),
      c(0.034952703615152854, 0.03484070894801079),
      c(0.034952703615152854, 0.03484070894801079),
      c(0.03614098127448581, 0.035312669921649095),
      c(0.035432630739150366, 0.034912158581439694),
      c(0.035432630739150366, 0.034912158581439694),
      c(0.035815727559287504, 0.03505127094114008),
      c(0.03540131462139614, 0.034738243333119145),
      c(0.03540131462139614, 0.034738243333119145)
    )),
    z = c(
      9.252619142383594, 12.078380683791904, 12.078380683791904, 5.565741269919053,
      5.457930240636103, 5.457930240636103, 9.051772787302813, 9.054982526543846,
      9.054982526543846
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
    test = maxcombo,
    cut = test_cutting(),
    seed = 2024,
    rho = c(0, 0),
    gamma = c(0, 0.5)
  )
  expected <- data.frame(
    method = rep("MaxCombo", 9L),
    parameter = rep("FH(0, 0) + FH(0, 0.5)", 9L),
    z = I(list(
      c(-3.7486049782713247, -4.149161171743935),
      c(-4.53034007934394, -4.778107819550277),
      c(-4.316452743033609, -4.2607297587160256),
      c(-3.4771440155825752, -3.605092910242299),
      c(-3.8631501353780324, -3.945081123231263),
      c(-3.2777779731288317, -2.919179640988388),
      c(-3.075862925191481, -3.1432278107909206),
      c(-3.619345457605645, -3.640458610667732),
      c(-4.2225917786532925, -4.243289152457)
    )),
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
    test = list(wlr_cut1, wlr_cut2, wlr_cut3),
    cut = test_cutting(),
    seed = 2024
  )
  expected <- data.frame(
    method = rep("WLR", 9L),
    parameter = rep(c("FH(rho=0, gamma=0)", "FH(rho=0, gamma=0.5)", "FH(rho=0.5, gamma=0)"), 3),
    estimation = c(
      -28.194825408790173, -24.448179085866208, -28.98456223760244,
      -26.84871111584948, -19.967690764549445, -23.830324019953483,
      -23.063020152157016, -18.05051250570102, -27.319166131937404
    ),
    se = c(
      7.521418120132856, 5.116707284384241, 6.918062043326721, 7.7214837796562,
      5.061414490811455, 6.931169838614448, 7.498065002594769, 4.958307300296486,
      6.9181407107482125
    ),
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
  skip_if_not_installed("gsDesign2")

  wlr_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  wlr_cut2 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))

  expect_error(
    sim_gs_n(
      n_sim = 3,
      sample_size = 400,
      enroll_rate = test_enroll_rate(),
      fail_rate = test_fail_rate(),
      test = list(wlr_cut1, wlr_cut2),
      cut = test_cutting(),
      seed = 2024
    ),
    "If you want to run different tests at each cutting"
  )
})

test_that("sim_gs_n() can combine wlr(), rmst(), and milestone() tests", {
  test_cut1 <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
  test_cut2 <- create_test(rmst, tau = 20)
  test_cut3 <- create_test(milestone, ms_time = 10)

  observed <- sim_gs_n(
    n_sim = 3,
    sample_size = 400,
    enroll_rate = test_enroll_rate(),
    fail_rate = test_fail_rate(),
    test = list(test_cut1, test_cut2, test_cut3),
    cut = test_cutting(),
    seed = 2024
  )
  expected <- data.frame(
    method = rep(c("WLR", "RMST", "milestone"), 3),
    parameter = rep(c("FH(rho=0, gamma=0)", "20", "10"), 3),
    estimation = c(
      -28.194825408790173, 2.8531309198097876, 0.1752731092436976,
      -26.84871111584948, 2.3996923551365086, 0.11738941400903585,
      -23.063020152157016, 2.4481174776559786, 0.1519773404060174
    ),
    se = I(list(
      7.521418120132856, 0.7147367073498636,
      c(0.034952703615152854, 0.03484070894801079), 7.7214837796562,
      0.7282783320819864, c(0.035432630739150366, 0.034912158581439694),
      7.498065002594769, 0.7416397937499536,
      c(0.03540131462139614, 0.034738243333119145)
    )),
    z = c(
      -3.7486049782713247, 3.991862864282945, 12.078380683791904,
      -3.4771440155825752, 3.2950209410683957, 5.457930240636103,
      -3.075862925191481, 3.3009521580248022, 9.054982526543846
    ),
    analysis = rep(1:3, 3),
    cut_date = c(24, 32, 45, 24, 32, 46.219327415802894, 24, 32, 50.86585486314699),
    sim_id = rep(1:3, each = 3L),
    n = rep(400L, 9L),
    event = c(229, 295, 355, 241, 290, 350, 226, 282, 350)
  )
  class(expected$se) <- NULL # I() adds the class "AsIs"
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
