library(ggplot2)
library(showtext)

# Pattern derived from <https://github.com/Ijeamakaanyene/patterns>

set.seed(42)

num_lines <- 3000
r_vals <- seq(0.04, 2, by = 0.1)
circles_list <- list()

for (i in seq_along(r_vals)) {
  r <- r_vals[i]
  circle <- data.frame(len = seq(0, 2 * pi, length.out = num_lines))
  circle$x <- r * sin(circle$len)
  circle$y <- r * cos(circle$len)

  circles_list[[i]] <- circle[sample(1:num_lines, num_lines / i, replace = FALSE), ]
}

final_circle <- data.table::rbindlist(circles_list)

font_add("Invention Bold", "inst/logo/Invention_Bd.ttf")
showtext_auto()

dots <- ggplot() +
  geom_point(
    data = final_circle,
    aes(x = x, y = y),
    color = "#f7f7f7",
    size = 0.1,
    alpha = 0.6
  ) +
  coord_fixed() +
  theme(
    plot.background = element_rect(fill = "#424242"),
    panel.background = element_rect(fill = "#424242"),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  annotate(
    "text",
    x = -0.5825, y = -0.2,
    label = "simtrial",
    family = "Invention Bold",
    size = 82,
    color = "#ffffff"
  )

ggsave(
  file.path(tempdir(), "logo-main.png"),
  device = ragg::agg_png,
  scale = 4, width = 256 * 2, height = 256 * 2,
  units = "px", dpi = 300
)

showtext_auto(FALSE)
