Tests and Coverage
================
23 August, 2020 13:01:17

  - [Coverage](#coverage)
  - [Unit Tests](#unit-tests)

This output is created by
[covrpage](https://github.com/metrumresearchgroup/covrpage).

## Coverage

Coverage summary is created using the
[covr](https://github.com/r-lib/covr) package.

| Object                                              | Coverage (%) |
| :-------------------------------------------------- | :----------: |
| simtrial                                            |    20.07     |
| [R/cutData.r](../R/cutData.r)                       |     0.00     |
| [R/cutDataAtCount.R](../R/cutDataAtCount.R)         |     0.00     |
| [R/getCutDateForCount.r](../R/getCutDateForCount.r) |     0.00     |
| [R/pMaxCombo.r](../R/pMaxCombo.r)                   |     0.00     |
| [R/pwexpfit.r](../R/pwexpfit.r)                     |     0.00     |
| [R/simfix.R](../R/simfix.R)                         |     0.00     |
| [R/simfix2simPWSurv.r](../R/simfix2simPWSurv.r)     |     0.00     |
| [R/simPWSurv.r](../R/simPWSurv.r)                   |     0.00     |
| [R/tenFH.R](../R/tenFH.R)                           |     0.00     |
| [R/tenFHCorr.r](../R/tenFHCorr.r)                   |     0.00     |
| [R/wMB.R](../R/wMB.R)                               |     0.00     |
| [R/rpwenroll.r](../R/rpwenroll.r)                   |    60.00     |
| [R/rpwexp.r](../R/rpwexp.r)                         |    88.24     |
| [R/tensurv.r](../R/tensurv.r)                       |    88.46     |
| [R/fixedBlockRand.r](../R/fixedBlockRand.r)         |    100.00    |

<br>

## Unit Tests

Unit Test summary is created using the
[testthat](https://github.com/r-lib/testthat) package.

| file                                                  | n | time | error | failed | skipped | warning | icon |
| :---------------------------------------------------- | -: | ---: | ----: | -----: | ------: | ------: | :--- |
| [testfixedBlockRand.r](testthat/testfixedBlockRand.r) | 1 | 0.01 |     0 |      0 |       0 |       0 |      |
| [testrpwenroll.r](testthat/testrpwenroll.r)           | 1 | 0.08 |     0 |      0 |       0 |       0 |      |
| [testrpwexp.r](testthat/testrpwexp.r)                 | 2 | 0.00 |     0 |      0 |       0 |       0 |      |
| [testtensurv.r](testthat/testtensurv.r)               | 6 | 0.27 |     0 |      0 |       0 |       2 | \-   |

<details open>

<summary> Show Detailed Test Results </summary>

| file                                                     | context            | test                                                     | status  | n | time | icon |
| :------------------------------------------------------- | :----------------- | :------------------------------------------------------- | :------ | -: | ---: | :--- |
| [testfixedBlockRand.r](testthat/testfixedBlockRand.r#L3) | testfixedBlockRand | fixedBlockRand returns an appropriate object             | PASS    | 1 | 0.01 |      |
| [testrpwenroll.r](testthat/testrpwenroll.r#L5)           | testrpwenroll      | rpwenroll handles 0 enrollment rate properly             | PASS    | 1 | 0.08 |      |
| [testrpwexp.r](testthat/testrpwexp.r#L6)                 | testrpwexp         | rpwexp handles 0 fail rate properly for one period       | PASS    | 1 | 0.00 |      |
| [testrpwexp.r](testthat/testrpwexp.r#L15)                | testrpwexp         | rpwexp handles 0 fail rate properly for multiple periods | PASS    | 1 | 0.00 |      |
| [testtensurv.r](testthat/testtensurv.r#L53)              | testtensurv        | Counting Process Format without ties                     | WARNING | 3 | 0.14 | \-   |
| [testtensurv.r](testthat/testtensurv.r#L70)              | testtensurv        | Counting Process Format with ties                        | WARNING | 3 | 0.13 | \-   |

| Failed | Warning | Skipped |
| :----- | :------ | :------ |
| \!     | \-      | \+      |

</details>

<details>

<summary> Session Info </summary>

| Field    | Value                          |
| :------- | :----------------------------- |
| Version  | R version 3.6.3 (2020-02-29)   |
| Platform | i386-w64-mingw32/i386 (32-bit) |
| Running  | Windows 10 x64 (build 17763)   |
| Language | English\_United States         |
| Timezone | America/New\_York              |

| Package  | Version |
| :------- | :------ |
| testthat | 2.3.2   |
| covr     | 3.5.0   |
| covrpage | 0.0.70  |

</details>

<!--- Final Status : skipped/warning --->
