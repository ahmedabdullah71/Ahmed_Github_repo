# Make the example reproducible each time the script is run.
set.seed(42)

random_numbers <- rnorm(50)

svg("test_plot.svg", width = 9, height = 6)
plot(
  random_numbers,
  type = "o",
  col = "steelblue",
  pch = 19,
  xlab = "Observation",
  ylab = "Random value",
  main = "Test Plot of Random Numbers"
)
abline(h = 0, col = "gray50", lty = 2)
dev.off()

cat("Plotted", length(random_numbers), "random numbers to test_plot.svg\n")
