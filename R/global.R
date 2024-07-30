#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the simtrial program.
#
#  simtrial is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# These global variables are declared to eliminate associated R CMD check warnings.
# There is no other identified functional impact of these global declarations.
# However, only variables defined in this package should be added to the list.
# If R CMD check warns about undefined global functions or variables that are
# from another package, they need to be added to this package's NAMESPACE.

utils::globalVariables(
  c(".",
    "analysis",
    "asy_n",
    "bound",
    "ex1_delayed_effect",
    "cte",
    "cross_lower",
    "cross_upper",
    "dropout_rate",
    "dropout_time",
    "duration",
    "ends_with",
    "enroll_time",
    "event",
    "eventCount",
    "events",
    "fail",
    "fail_time",
    "finish",
    "hr",
    "i",
    "lambda",
    "lower",
    "lower_bound",
    "max_weight",
    "mtte",
    "N",
    "n",
    "n_cross_lower",
    "n_cross_upper",
    "n_event_tol",
    "n_risk_tol",
    "n_risk_trt",
    "o_minus_e",
    "one",
    "origin",
    "period",
    "probability",
    "rate",
    "s",
    "sim_id",
    "sim_lower_prob",
    "sim_n",
    "sim_upper_prob",
    "status",
    "stratum",
    "time",
    "treatment",
    "tte",
    "upper",
    "upper_bound",
    "v",
    "var_o_minus_e",
    "weight",
    "z"
  )
)

# Workaround to remove `R CMD check` NOTE "All declared Imports should be used."
# https://r-pkgs.org/dependencies-in-practice.html#how-to-not-use-a-package-in-imports
ignore_unused_imports <- function() {
  utils::globalVariables
}
