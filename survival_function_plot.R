# ── Parameters ──────────────────────────────────────────────
c0  <- 100
lam <- 0.08
mu  <- 0.0025

# ── Survival function (piecewise) ───────────────────────────
# 0 < t < c0:  S(t) = 1 - mu*t/(lam+mu*t) * exp(-lam*(c0/t - 1))
# t >= c0:     S(t) = lam * exp(mu*c0 - mu*t) / (lam + mu*t)
S <- function(t) {
  result <- numeric(length(t))
  m1 <- t > 0 & t < c0
  m2 <- t >= c0
  t1 <- t[m1]
  result[m1] <- 1 - (mu * t1) / (lam + mu * t1) * exp(-lam * (c0 / t1 - 1))
  t2 <- t[m2]
  result[m2] <- lam * exp(mu * c0 - mu * t2) / (lam + mu * t2)
  result
}

# ── x range ─────────────────────────────────────────────────
tau_scan    <- seq(0.1, 6000, length.out = 600000)
x_max       <- max(tau_scan[S(tau_scan) >= 1e-6]) * 1.10
lam_over_mu <- lam / mu     # 32
one_over_mu <- 1.0 / mu     # 400

rb        <- c(0, lam_over_mu, c0, one_over_mu, one_over_mu * 1.3, x_max)
bg_colors <- c("#fce8e8", "#fdefd8", "#fdf7de", "#e8f0fb", "#e8f5e8")

t_plot <- seq(0.3, x_max, length.out = 8000)
s_plot <- S(t_plot)
keep   <- s_plot >= 1e-6
t_plot <- t_plot[keep]
s_plot <- s_plot[keep]

y_min <- 1e-6

# ── Open PNG ─────────────────────────────────────────────────
png("survival_function_plot.png",
    width = 2000, height = 1100, res = 200, bg = "white")

par(mar = c(5.5, 5.2, 6.5, 2))

# Empty plot on log-y scale
plot(NA, NA,
     xlim = c(0, x_max), ylim = c(y_min, 1.5),
     log  = "y",
     xlab = "", ylab = "",
     axes = FALSE, frame.plot = FALSE)

# ── Background shading ───────────────────────────────────────
for (i in seq_along(bg_colors)) {
  rect(rb[i], 1e-9, rb[i + 1], 1e2,
       col = bg_colors[i], border = NA)
}

# ── Boundary dashed lines ────────────────────────────────────
abline(v = rb[2:(length(rb) - 1)],
       lty = 2, col = "#55555599", lwd = 1.2)

# ── Main survival curve ──────────────────────────────────────
lines(t_plot, s_plot, col = "#1a3f7a", lwd = 2.3)

# ── Axes ─────────────────────────────────────────────────────
x_ticks <- seq(0, floor(x_max / 400) * 400, by = 400)
axis(1, at = x_ticks, cex.axis = 0.85, lwd = 1.2)

y_ticks  <- 10^seq(-6, 0)
y_labels <- c(expression(10^-6), expression(10^-5), expression(10^-4),
              expression(10^-3), expression(10^-2), expression(10^-1),
              expression(10^0))
axis(2, at = y_ticks, labels = y_labels, las = 1, cex.axis = 0.85, lwd = 1.2)
box(bty = "l", lwd = 1.2)

mtext("Time in antibiotic,  t", side = 1, line = 3.5, cex = 1.0)
mtext(expression("Survival fraction  " * italic(S(t))),
      side = 2, line = 3.8, cex = 1.0)

# ── Boundary labels below x-axis ────────────────────────────
mtext(expression(lambda / mu), side = 1, at = lam_over_mu, line = 4.6, cex = 0.92, col = "#444444")
mtext(expression(c[0]),        side = 1, at = c0,           line = 4.6, cex = 0.92, col = "#444444")
mtext(expression(1 / mu),      side = 1, at = one_over_mu,  line = 4.6, cex = 0.92, col = "#444444")

# ── Top region labels above panel ───────────────────────────
usr <- par("usr")  # c(x1, x2, log10(y1), log10(y2)) when log="y"
y_top_line1 <- 10^(usr[4] + 0.52)
y_top_line2 <- 10^(usr[4] + 0.13)

cx3 <- (c0 + one_over_mu) / 2
cx5 <- (one_over_mu * 1.3 + x_max) / 2

text(cx3, y_top_line1, "III. Apparent power-law-like",
     col = "#9e6000", cex = 0.9, font = 2, xpd = TRUE, adj = c(0.5, 0))
text(cx3, y_top_line2, expression(c[0] ~ "<" ~ t ~ "<" ~ 1/mu),
     col = "#9e6000", cex = 0.85, xpd = TRUE, adj = c(0.5, 0))

text(cx5, y_top_line1, "V. True exponential tail",
     col = "#1a7a1a", cex = 0.9, font = 2, xpd = TRUE, adj = c(0.5, 0))
text(cx5, y_top_line2, expression(t ~ ">>" ~ 1/mu),
     col = "#1a7a1a", cex = 0.85, xpd = TRUE, adj = c(0.5, 0))

# ── Title ────────────────────────────────────────────────────
title(main = expression(bold(
  "Predicted survival function " * italic(S(t)) * "  -  linear x, log y"
)), cex.main = 1.1, line = 5.0)

# ── Watermark ────────────────────────────────────────────────
text(x_max * 0.99, y_min * 1.2, "x-axis linear; y-axis log",
     adj = c(1, 0), col = "#aaaaaa", cex = 0.65, xpd = TRUE)

dev.off()
message("Saved survival_function_plot.png")