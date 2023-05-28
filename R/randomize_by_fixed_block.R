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

#' Permuted fixed block randomization
#'
#' Fixed block randomization. The `block` input should repeat each
#' treatment code the number of times it is to be included within each block.
#' The final block will be a partial block if `n` is not an exact multiple
#' of the block length.
#'
#' @param n Sample size to be randomized.
#' @param block Vector of treatments to be included in each block.
#'
#' @return A treatment group sequence (vector) of length `n` with
#'   treatments from `block` permuted within each block having
#'   block size equal to the length of `block`.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Example 1
#' # 2:1 randomization with block size 3, treatments "A" and "B"
#' tibble(x = 1:10) %>% mutate(Treatment = randomize_by_fixed_block(block = c("A", "B", "B")))
#'
#' # Example 2
#' # Stratified randomization
#' tibble(stratum = c(rep("A", 10), rep("B", 10))) %>%
#'   group_by(stratum) %>%
#'   mutate(Treatment = randomize_by_fixed_block())
randomize_by_fixed_block <- function(n = 10, block = c(0, 0, 1, 1)) {
  length_block <- length(block)
  n_block <- ceiling(n / length_block)
  n_total <- n_block * length_block

  block_to_sample <- rep(block, each = n_block)
  sample_order <- rep(1:n_block, length_block)
  sample_order <- sample_order + stats::runif(n_total)

  (block_to_sample[order(sample_order)])[1:n]
}
