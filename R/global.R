#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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

utils::globalVariables(
  c(
    "atrisk",
    "Count",
    "cte",
    "dropoutRate",
    "dropoutTime",
    "enrollTime",
    "eventCount",
    "events",
    "Ex1delayedEffect",
    "fail",
    "fail_rate",
    "fail_time",
    "hr",
    "duration",
    "enrollTime",
    "event",
    "finish",
    "lambda",
    "mtte",
    "N",
    "OminusE",
    "one",
    "origin",
    "period",
    "rate",
    "s",
    "S",
    "status",
    "stratum",
    "time",
    "treatment",
    "tte",
    "txevents",
    "Var",
    "w",
    "wOminusE",
    "wVar",
    "txatrisk"
  )
)
