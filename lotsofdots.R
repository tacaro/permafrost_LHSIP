library(tidyverse)

df <- crossing(
  x = seq(1, 100, by = 1),
  y = seq(1, 100, by = 1)
) |> 
  mutate(
    z = case_when(
      x == 70 & y == 23 ~ FALSE,
      TRUE ~ TRUE)
  )

df |> 
  ggplot(
    aes(x, y, color = z)
  ) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("red", "white")) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "None",
    panel.background = element_rect(fill = "black")
  )