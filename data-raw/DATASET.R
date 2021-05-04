## code to prepare `DATASET` dataset goes here
library(tidyr)
set.seed(6671)
ds <- simPWSurv(n = 200,
                enrollRates = tibble(rate = 200 / 12, duration = 12),
                failRates = tribble(
                  ~Stratum, ~Period, ~Treatment,     ~duration, ~rate,
                  "All",        1,   "Control",      42,        log(2) / 15,
                  "All",        1,   "Experimental", 6,         log(2) / 15,
                  "All",        2,   "Experimental", 36,        log(2) / 15 * 0.6),
                dropoutRates = tribble(
                  ~Stratum, ~Period, ~Treatment,     ~duration, ~rate,
                  "All",        1,   "Control",      42,        0,
                  "All",        1,   "Experimental", 42,        0)
)
# cut data at 24 months after final enrollment
MBdelayed <- ds %>% cutData(max(ds$enrollTime) + 24)

usethis::use_data("MBdelayed")
