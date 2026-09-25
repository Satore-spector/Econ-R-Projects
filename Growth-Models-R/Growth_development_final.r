# Solow & RCK growth models - Republic of Korea
# data: Penn World Tables 11.0
#
# models: Solow growth model, RCK model, Barro convergence regressions
#
# PWT variables used:
#   rgdpe  - expenditure-side real GDP (2017 USD PPP)
#   cn     - capital stock at current PPPs
#   emp    - employment (used as labour force proxy)
#   hc     - human capital index
#   csh_i  - investment share of GDP (savings rate proxy)
#   irr    - real internal rate of return
#   delta  - average depreciation rate
#   csh_c  - consumption share of GDP
#   rtfpna - TFP at constant national prices (used for g)
#   rconna - real consumption (national accounts)
#   labsh  - labour share of output (used to derive capital share alpha)


# libraries

library(tidyverse)   # data wrangling and ggplot2
library(readxl)      # reading PWT xlsx file
library(estimatr)    # lm_robust for HC-robust regressions
library(zoo)         # rollmean for conditional Barro regression
library(writexl)     # optional: exporting to Excel
library(ggthemes)    # WSJ theme and palette


# WSJ colour palette, defined once so I can retheme every plot from here
# navy = primary series, red = secondary/reference, black = theoretical loci
# green = third series, grey = zero lines/neutral, fill = CI ribbon

WSJ_NAVY  <- "#0E4D78"
WSJ_RED   <- "#C1272D"
WSJ_BLACK <- "#000000"
WSJ_GREEN <- "#336633"
WSJ_GREY  <- "#666666"
WSJ_FILL  <- "#DCCFCC"

# legend position, applied after theme_wsj() on every plot below
LEGEND_POS <- "bottom"


# 1. load and filter data

pwt <- read_excel(
  "C:\\Users\\teert\\OneDrive\\Documents\\Data\\Growth_development_data\\solow\\pwt110.xlsx",
  sheet = "Data"
)

# filter to korea, grab the variables I need, rename so they make sense
pwt_korea <- pwt |>
  filter(country == "Republic of Korea") |>
  select(
    year, rgdpe, cn, pop, hc, csh_i, irr, delta,
    emp, csh_c, rtfpna, rconna, labsh
  ) |>
  rename(
    real_gdp                  = rgdpe,
    capital_stock             = cn,
    population                = pop,
    human_capital_index       = hc,
    s                         = csh_i,   # investment share (savings rate proxy)
    real_interest_rate        = irr,
    depreciation_rate         = delta,
    labour_force              = emp,
    consumption_share         = csh_c,
    total_factor_productivity = rtfpna,
    consumption               = rconna,
    labour_share_output       = labsh
  ) |>
  slice(-1:-9)  # dropping pre-1960, data gets patchy that far back


# 2. build per-worker terms and growth rates

pwt_korea <- pwt_korea |>
  mutate(
    y             = real_gdp / labour_force,             # output per worker
    k             = capital_stock / labour_force,        # capital per worker
    n             = (labour_force - lag(labour_force)) / lag(labour_force), # labour force growth
    growth_rate   = log(real_gdp) - log(lag(real_gdp)), # log GDP growth rate
    g             = (total_factor_productivity - lag(total_factor_productivity)) /
                      lag(total_factor_productivity),    # TFP growth (national prices)
    capital_share = 1 - labour_share_output              # α = 1 - labour share
  )


# 3. calibrated parameters
# g = long-run mean TFP growth (rtfpna). this mixes catch-up growth with
#     actual frontier growth, so it's a simplification, but it's the
#     standard approach given what PWT gives us
# alpha = mean of (1 - labour share) over the sample
# theta = 2 (inverse IES), pretty standard value in the RCK literature
# rho and A get backed out below from the Euler equation / production function

x     <- mean(pwt_korea$g, na.rm = TRUE)                 # mean TFP growth (g ≈ 0.0138)
alpha <- mean(pwt_korea$capital_share, na.rm = TRUE)      # capital share (α ≈ 0.35)
theta <- 2                                                 # inverse elasticity of substitution

# back out A from y = A * k^alpha  ->  A = y / k^alpha
pwt_korea <- pwt_korea |>
  mutate(A_estimate = y / (k ^ alpha))

A <- mean(pwt_korea$A_estimate, na.rm = TRUE)


# 4. solow model
# steady state: s*f(k) = (n + g + delta)*k
# delta_k = s*y - (n + g + delta)*k, positive means we're below steady state

pwt_korea <- pwt_korea |>
  mutate(
    actual_investment     = s * y,                             # s * f(k)
    break_even_investment = (n + x + depreciation_rate) * k,  # (n + g + δ)*k
    delta_k               = actual_investment - break_even_investment  # predicted Δk
  )


# 5. RCK model
# euler equation: g_c = (1/theta)*(r - rho)  ->  rho = r - theta*g_c
# k-dot = 0 locus: c = f(k) - (n + g + delta)*k
# c-dot = 0 locus: vertical line at k* where f'(k*) = rho + theta*g
# cobb-douglas gives f'(k) = alpha*A*k^(alpha-1), solve for k*:
# k* = (alpha*A / (rho + theta*g)) ^ (1/(1-alpha))

pwt_korea <- pwt_korea |>
  mutate(
    c                        = consumption / labour_force,                        # consumption per worker
    g_c_dot                  = (c - lag(c)) / lag(c),                            # consumption growth rate
    estimated_rho            = real_interest_rate - (g_c_dot * theta),           # implied discount rate
    rck_capital_accumulation = y - c - ((n + x + depreciation_rate) * k),        # k-dot
    c_k_dot_zero             = y - ((n + x + depreciation_rate) * k)             # k-dot = 0 locus
  )

# back out mean rho from the euler equation
rho <- mean(pwt_korea$estimated_rho, na.rm = TRUE)

# k* where c-dot = 0
k_star <- (alpha * A / (rho + theta * x)) ^ (1 / (1 - alpha))


# 6. barro convergence regressions
# growth = b0 + b1*lag(y) + b2*lag(controls) + e
# note: this is one country over time, not a cross-country panel, so it's
# really testing mean reversion within Korea, not "convergence" in the
# textbook cross-country Barro sense. same idea though

# 6a. unconditional
barro_uc_reg <- lm_robust(
  growth_rate ~ lag(y),
  data    = pwt_korea,
  se_type = "HC1"
)
summary(barro_uc_reg)

# grab fitted values so I can plot them
pwt_korea_barro <- pwt_korea |>
  mutate(g_cap = NA_real_) |>
  slice(-1) |>   # drop first row (NA from lag)
  mutate(g_cap = fitted(barro_uc_reg))

# plot 6a - unconditional convergence scatter
ggplot(data = pwt_korea_barro, aes(x = lag(y), y = growth_rate)) +
  geom_point(aes(color = "Observations"), size = 2, alpha = 0.7) +
  geom_smooth(
    aes(color = "OLS Trend (95% CI)"),
    method = "lm", fill = WSJ_FILL, se = TRUE
  ) +
  scale_color_manual(
    name   = "",
    values = c("Observations" = WSJ_NAVY, "OLS Trend (95% CI)" = WSJ_RED)
  ) +
  labs(
    title    = "Unconditional Convergence: Republic of Korea",
    subtitle = "Barro Regression — Growth Rate vs. Lagged Output per Worker",
    x        = "Lagged Output per Worker",
    y        = "Growth Rate",
    caption  = "Shaded band: 95% confidence interval"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 6b - actual vs predicted growth (unconditional)
ggplot(data = pwt_korea_barro, aes(x = year)) +
  geom_line(aes(y = growth_rate, color = "Actual Growth Rate"),    linewidth = 1) +
  geom_line(aes(y = g_cap,       color = "Predicted Growth Rate"), linewidth = 1) +
  scale_color_manual(
    name   = "",
    values = c("Actual Growth Rate" = WSJ_NAVY, "Predicted Growth Rate" = WSJ_RED)
  ) +
  labs(
    title    = "Actual vs. Predicted Growth Rate (Unconditional)",
    subtitle = "Regression-predicted growth rate against realised growth",
    x        = "Year",
    y        = "Growth Rate"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 6b. conditional, controls are a 5yr rolling mean of the real interest rate,
# human capital, and depreciation. depreciation barely moves in PWT so it's
# probably not doing much here but leaving it in
barro_c_reg <- lm_robust(
  growth_rate ~ lag(y) +
    lag(rollmean(real_interest_rate, 5, fill = NA)) +
    lag(human_capital_index) +
    lag(depreciation_rate),
  data    = pwt_korea,
  se_type = "HC1"
)
summary(barro_c_reg)

# Attach fitted values for plotting
pwt_korea_barro_c <- pwt_korea_barro |>
  mutate(g_cap_c = NA_real_) |>
  slice(-1:-3) |>  # drop rows with NAs from rolling mean lag
  mutate(g_cap_c = fitted(barro_c_reg))

# plot 6c - actual vs predicted growth, conditional model fit
ggplot(data = pwt_korea_barro_c, aes(x = g_cap_c, y = growth_rate)) +
  geom_point(aes(color = "Observations"), size = 2, alpha = 0.7) +
  geom_smooth(aes(color = "OLS Trend"), method = "lm", se = FALSE) +
  scale_color_manual(
    name   = "",
    values = c("Observations" = WSJ_NAVY, "OLS Trend" = WSJ_RED)
  ) +
  labs(
    title   = "Conditional Model Fit: Actual vs. Predicted Growth",
    x       = "Predicted Growth Rate (Conditional on Controls)",
    y       = "Actual Growth Rate"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 6d - human capital vs predicted growth
ggplot(data = pwt_korea_barro_c, aes(x = g_cap_c, y = human_capital_index)) +
  geom_point(aes(color = "Observations"), size = 2, alpha = 0.7) +
  geom_smooth(aes(color = "OLS Trend"), method = "lm", se = FALSE) +
  scale_color_manual(
    name   = "",
    values = c("Observations" = WSJ_NAVY, "OLS Trend" = WSJ_RED)
  ) +
  labs(
    title   = "Human Capital vs. Predicted Growth (Conditional)",
    x       = "Predicted Growth Rate (Conditional on Controls)",
    y       = "Human Capital Index"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 6e - actual vs predicted growth (conditional)
ggplot(data = pwt_korea_barro_c, aes(x = year)) +
  geom_line(aes(y = growth_rate, color = "Actual Growth Rate"),    linewidth = 1) +
  geom_line(aes(y = g_cap_c,     color = "Predicted Growth Rate"), linewidth = 1) +
  scale_color_manual(
    name   = "",
    values = c("Actual Growth Rate" = WSJ_NAVY, "Predicted Growth Rate" = WSJ_RED)
  ) +
  labs(
    title    = "Actual vs. Predicted Growth Rate (Conditional)",
    subtitle = "Regression-predicted growth rate against realised growth",
    x        = "Year",
    y        = "Growth Rate"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 7. solow plots
# f(k) = A * K^alpha * L^(1-alpha), steady state where delta_k = 0

# plot 7a - change in capital per worker, smoothed
ggplot(data = pwt_korea, aes(x = year, y = delta_k)) +
  geom_point(aes(color = "Annual Δk"),    size = 2, alpha = 0.7) +
  geom_smooth(aes(color = "Smoothed Trend"), linewidth = 1, alpha = 0.7, se = FALSE) +
  scale_color_manual(
    name   = "",
    values = c("Annual Δk" = WSJ_NAVY, "Smoothed Trend" = WSJ_RED)
  ) +
  labs(
    title   = "Solow Model: Change in Capital per Worker (Δk)",
    x       = "Year",
    y       = "Δk (Capital Accumulation per Worker)",
    caption = "Δk > 0 indicates economy is below steady state."
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 7b - same thing but as a connected line instead of smoothed
ggplot(data = pwt_korea, aes(x = year, y = delta_k)) +
  geom_point(aes(color = "Annual Observations"), size = 2, alpha = 0.7) +
  geom_line(aes(color  = "Year-on-Year Path"),   linewidth = 1, alpha = 0.7) +
  scale_color_manual(
    name   = "",
    values = c("Annual Observations" = WSJ_RED, "Year-on-Year Path" = WSJ_NAVY)
  ) +
  labs(
    title   = "Solow Model: Change in Capital per Worker (Δk)",
    x       = "Year",
    y       = "Δk (Capital Accumulation per Worker)",
    caption = "Δk > 0 indicates economy is below steady state."
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 7c - actual vs break-even investment, smoothed
ggplot(data = pwt_korea, aes(x = k)) +
  geom_point(aes(y = actual_investment,    color = "Actual Investment s·f(k)"),   size = 2, alpha = 0.7) +
  geom_smooth(aes(y = actual_investment,   color = "Actual Investment s·f(k)"),   linewidth = 1, se = FALSE) +
  geom_point(aes(y = break_even_investment, color = "Break-Even (n+g+δ)·k"),      size = 2) +
  geom_smooth(aes(y = break_even_investment, color = "Break-Even (n+g+δ)·k"),     linewidth = 1, se = FALSE) +
  scale_color_manual(
    name   = "",
    values = c("Actual Investment s·f(k)" = WSJ_NAVY, "Break-Even (n+g+δ)·k" = WSJ_RED)
  ) +
  labs(
    title   = "Solow Model: Actual vs. Break-Even Investment",
    x       = "Capital per Worker (k)",
    y       = "Investment per Worker",
    caption = "Steady state where s·y = (n + g + δ)·k"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 7d - same, connected line
ggplot(data = pwt_korea, aes(x = k)) +
  geom_point(aes(y = actual_investment,     color = "Actual Investment s·f(k)"),  size = 2, alpha = 0.7) +
  geom_line(aes(y = actual_investment,      color = "Actual Investment s·f(k)"),  linewidth = 1) +
  geom_point(aes(y = break_even_investment,  color = "Break-Even (n+g+δ)·k"),     size = 2) +
  geom_line(aes(y = break_even_investment,   color = "Break-Even (n+g+δ)·k"),     linewidth = 1) +
  scale_color_manual(
    name   = "",
    values = c("Actual Investment s·f(k)" = WSJ_NAVY, "Break-Even (n+g+δ)·k" = WSJ_RED)
  ) +
  labs(
    title   = "Solow Model: Actual vs. Break-Even Investment",
    x       = "Capital per Worker (k)",
    y       = "Investment per Worker",
    caption = "Steady state where s·y = (n + g + δ)·k"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 7e - output per worker vs steady state
y_star <- A * k_star ^ alpha

ggplot(data = pwt_korea, aes(x = year)) +
  geom_line(aes(y = y, color = "Output per Worker (y)"), linewidth = 1) +
  geom_hline(
    aes(yintercept = y_star, color = "y* Steady State"),
    linetype  = "dashed",
    linewidth = 0.8
  ) +
  scale_color_manual(
    name   = "",
    values = c("Output per Worker (y)" = WSJ_NAVY, "y* Steady State" = WSJ_RED)
  ) +
  annotate(
    "text",
    x     = min(pwt_korea$year) + 1,
    y     = y_star * 1.04,
    label = "y* (steady state)",
    color = WSJ_RED,
    hjust = 0,
    size  = 3.5
  ) +
  labs(
    title   = "Output per Worker vs. Steady-State Level: Republic of Korea",
    x       = "Year",
    y       = "Output per Worker (y)",
    caption = "Red dashed line: implied steady-state output y* = A·(k*)^α"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 7f - capital per worker converging to steady state
ggplot(data = pwt_korea, aes(x = year)) +
  geom_line(aes(y = k, color = "Capital per Worker (k)"), linewidth = 1) +
  geom_hline(
    aes(yintercept = k_star, color = "k* Steady State"),
    linetype  = "dashed",
    linewidth = 0.8
  ) +
  scale_color_manual(
    name   = "",
    values = c("Capital per Worker (k)" = WSJ_NAVY, "k* Steady State" = WSJ_RED)
  ) +
  annotate(
    "text",
    x     = min(pwt_korea$year) + 1,
    y     = k_star * 1.03,
    label = "k* (steady state)",
    color = WSJ_RED,
    hjust = 0,
    size  = 3.5
  ) +
  labs(
    title   = "Capital per Worker Convergence: Republic of Korea",
    x       = "Year",
    y       = "Capital per Worker (k)",
    caption = "Red dashed line: steady-state capital k*"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 8. RCK plots

# plot 8a - k-dot over time
ggplot(data = pwt_korea, aes(x = year, y = rck_capital_accumulation)) +
  geom_point(aes(color  = "Annual Observations"), size = 1.5) +
  geom_line(aes(color   = "Year-on-Year Path"),   linewidth = 0.8) +
  geom_smooth(aes(color = "Smoothed Trend"),       linewidth = 1, se = FALSE) +
  scale_color_manual(
    name   = "",
    values = c(
      "Annual Observations" = WSJ_RED,
      "Year-on-Year Path"   = WSJ_NAVY,
      "Smoothed Trend"      = WSJ_BLACK
    )
  ) +
  labs(
    title = "RCK Model: Capital Accumulation (k-dot) over Time",
    x     = "Year",
    y     = "k-dot  (y − c − (n + g + δ)·k)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 8b - empirical path, consumption vs capital, with arrows showing direction
# points are coloured by year, grey path is just the chronological trajectory
ggplot(data = pwt_korea, aes(x = k, y = c)) +
  geom_smooth(color = WSJ_GREY, se = FALSE, linewidth = 0.8) +
  geom_path(
    color = WSJ_GREY,
    arrow = arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  geom_point(aes(color = year), size = 2.5) +
  scale_color_gradient(
    low  = WSJ_RED,
    high = WSJ_NAVY,
    name = "Year"
  ) +
  labs(
    title    = "Empirical RCK Path: Republic of Korea",
    subtitle = "Coloured points: economy by year | Grey path: chronological trajectory | Grey line: smoothed trend",
    x        = "Capital per Worker (k)",
    y        = "Consumption per Worker (c)",
    caption  = "Arrows indicate chronological direction."
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 8c - phase diagram, k-dot = 0 locus plus the actual trajectory
# black triangles/line = k-dot = 0 locus, grey path = saddle path
ggplot(data = pwt_korea, aes(x = k)) +

  # k-dot = 0 locus (inverted-U curve)
  geom_point(aes(y = c_k_dot_zero), color = WSJ_BLACK, size = 2.5, shape = "triangle") +
  geom_line(aes(y = c_k_dot_zero),  color = WSJ_BLACK, linewidth = 0.8) +

  # Actual saddle path with directional arrows
  geom_path(
    aes(y = c),
    color = WSJ_GREY,
    arrow = arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  geom_point(aes(y = c, color = year), size = 2.5) +

  scale_color_gradient(
    low  = WSJ_RED,
    high = WSJ_NAVY,
    name = "Year"
  ) +
  labs(
    title    = "Empirical Phase Diagram: Republic of Korea",
    subtitle = "Black triangles/line: k-dot = 0 locus | Coloured points: actual economy path | Grey: trajectory",
    x        = "Capital per Worker (k)",
    y        = "Consumption per Worker (c)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 8d - same phase diagram but smoothed, easier to read
ggplot(data = pwt_korea, aes(x = k)) +
  geom_smooth(
    aes(y = c_k_dot_zero, color = "k-dot = 0 Locus"),
    se        = FALSE,
    linetype  = "dashed",
    linewidth = 1
  ) +
  geom_smooth(
    aes(y = c, color = "Consumption Path"),
    se        = FALSE,
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "",
    values = c("k-dot = 0 Locus" = WSJ_BLACK, "Consumption Path" = WSJ_NAVY)
  ) +
  labs(
    title    = "Empirical Phase Diagram: Republic of Korea",
    subtitle = "Smoothed loci — dashed: k-dot = 0 locus | solid: consumption path",
    x        = "Capital per Worker (k)",
    y        = "Consumption per Worker (c)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 8e - full phase diagram: both loci plus the consumption path
ggplot(data = pwt_korea, aes(x = k)) +

  # k-dot = 0 locus
  geom_smooth(
    aes(y = c_k_dot_zero, color = "k-dot = 0 Locus"),
    se        = FALSE,
    linetype  = "dashed",
    linewidth = 1
  ) +

  # consumption path
  geom_smooth(
    aes(y = c, color = "Consumption Path"),
    se        = FALSE,
    linewidth = 1
  ) +

  # c-dot = 0 locus, vertical line at k*
  geom_vline(
    aes(xintercept = k_star, color = "c-dot = 0 Locus (k*)"),
    linetype  = "dashed",
    linewidth = 0.8
  ) +

  scale_color_manual(
    name   = "",
    values = c(
      "k-dot = 0 Locus"      = WSJ_BLACK,
      "Consumption Path"     = WSJ_NAVY,
      "c-dot = 0 Locus (k*)" = WSJ_RED
    )
  ) +

  # label k* on the plot
  annotate(
    "text",
    x     = k_star * 1.03,
    y     = max(pwt_korea$c, na.rm = TRUE) * 0.5,
    label = "ċ = 0\n(k*)",
    color = WSJ_RED,
    hjust = 0,
    size  = 3.5
  ) +

  labs(
    title    = "Empirical Phase Diagram: Republic of Korea",
    subtitle = "Full RCK phase diagram with both nullclines and saddle path",
    x        = "Capital per Worker (k)",
    y        = "Consumption per Worker (c)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 9. diagnostic plots

# plot 9a - cobb-douglas fit vs actual GDP, just a sanity check
ggplot(data = pwt_korea, aes(x = year)) +
  geom_line(
    aes(y = real_gdp, color = "Actual Real GDP"),
    linewidth = 1.5
  ) +
  geom_line(
    aes(
      y     = A * (capital_stock ^ alpha) * (labour_force ^ (1 - alpha)),
      color = "Cobb-Douglas f(k)"
    ),
    linewidth = 1.5
  ) +
  scale_color_manual(
    name   = "",
    values = c("Actual Real GDP" = WSJ_RED, "Cobb-Douglas f(k)" = WSJ_NAVY)
  ) +
  labs(
    title    = "Actual vs. Cobb-Douglas Production Function",
    subtitle = "Plotting the production function against real GDP",
    x        = "Year",
    y        = "Output"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)

# plot 9b - raw delta k vs solow delta k vs RCK k-dot, all together
ggplot(data = pwt_korea, aes(x = year)) +
  geom_line(aes(y = (k - lag(k)),            color = "Raw Δk/L"),  linewidth = 1) +
  geom_line(aes(y = delta_k,                  color = "Solow Δk"), linewidth = 1) +
  geom_line(aes(y = rck_capital_accumulation, color = "RCK k-dot"), linewidth = 1) +
  geom_hline(
    aes(yintercept = 0, color = "Zero Reference"),
    linetype  = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(
    name   = "",
    values = c(
      "Raw Δk/L"      = WSJ_NAVY,
      "Solow Δk"      = WSJ_GREEN,
      "RCK k-dot"     = WSJ_RED,
      "Zero Reference" = WSJ_GREY
    )
  ) +
  labs(
    title    = "Three Measures of Capital Accumulation per Worker",
    subtitle = "Comparison of raw data, Solow model, and RCK model capital accumulation",
    x        = "Year",
    y        = "Change in Capital per Worker (USD PPP)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 10. OLG model
# law of motion: k(t+1) = ((1-delta)*k(t) + i(t)) / (1+n(t))
# steady state is where the transition curve crosses the 45 degree line
# dynamic inefficiency (over-saving) shows up when r < n, per Aaron 1966

# 10a. build the OLG variables

pwt_korea <- pwt_korea |>
  mutate(
    i = actual_investment, # investment per worker (already computed)
    k_next = ((1 - depreciation_rate) * k + i) / (1 + n), # OLG law of motion
    ln_k = log(k), # log capital per worker
    ln_k_next = log(k_next), # log next-period capital per worker
    fortyfive = ln_k # 45-degree reference line
  )

# 10b. transition function, log-linear regression
# ln(k_t+1) = b0 + b1*ln(k_t), b1 < 1 means it converges to a steady state

olg_reg <- lm(ln_k_next ~ ln_k, data = pwt_korea)
summary(olg_reg)

pwt_korea <- pwt_korea |>
  mutate(ln_k_next_hat = predict(olg_reg, newdata = pwt_korea))

# 10c. find years close to steady state, |ln(k_t+1) - ln(k_t)| < 0.01

olg_ss <- pwt_korea |>
  filter(!is.na(ln_k_next_hat)) |>
  mutate(diff = ln_k_next_hat - ln_k) |>
  filter(abs(diff) < 0.01) |>
  select(year, ln_k, diff)

cat("OLG Steady-State Years (|ln k_next_hat - ln k| < 0.01):\n")
print(olg_ss)

# 10d. dynamic inefficiency check
# flag years where r < n + g -> economy is over-saving
pwt_korea <- pwt_korea |>
  mutate(
    total_growth_rate = n + x, # x is your calibrated mean g
    dyn_ineff = real_interest_rate < total_growth_rate,
    dyn_ineff_label = if_else(
      dyn_ineff,
      "r < n+g  (Dynamically Inefficient)",
      "r ≥ n+g  (Dynamically Efficient)"
    )
  )

# plot 10a - OLG phase diagram, raw trajectory vs 45 degree line
# where they cross is the steady state

ggplot(
  data = pwt_korea |> filter(!is.na(ln_k), !is.na(ln_k_next)) |> arrange(ln_k),
  aes(x = ln_k)
) +
  geom_line(aes(y = ln_k_next, color = "k(t+1) — OLG Transition"), linewidth = 1) +
  geom_line(aes(y = fortyfive, color = "45° Line  (k(t+1) = k(t))"),
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "k(t+1) — OLG Transition" = WSJ_NAVY,
      "45° Line  (k(t+1) = k(t))" = WSJ_RED
    )
  ) +
  labs(
    title    = "OLG Phase Diagram: Republic of Korea",
    subtitle = "Steady state where the transition curve crosses the 45° line",
    x        = "ln(k_t)  — Log Capital per Worker",
    y        = "ln(k_t+1)  — Log Next-Period Capital per Worker",
    caption  = "Law of motion: k(t+1) = [(1 - δ)·k + i] / (1 + n)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# plot 10b - fitted transition curve instead of raw, easier to read

ggplot(
  data = pwt_korea |> filter(!is.na(ln_k), !is.na(ln_k_next_hat)) |> arrange(ln_k),
  aes(x = ln_k)
) +
  geom_line(aes(y = ln_k_next_hat, color = "Fitted Transition  ln(k_t+1) = β0 + β1·ln(k_t)"),
    linewidth = 1.2
  ) +
  geom_line(aes(y = fortyfive, color = "45° Line  (k(t+1) = k(t))"),
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "Fitted Transition  ln(k_t+1) = β0 + β1·ln(k_t)" = WSJ_NAVY,
      "45° Line  (k(t+1) = k(t))" = WSJ_RED
    )
  ) +
  labs(
    title = "OLG Fitted Transition Function: Republic of Korea",
    subtitle = paste0(
      "β1 = ", round(coef(olg_reg)["ln_k"], 3),
      " — slope < 1 confirms convergence to steady state"
    ),
    x = "ln(k_t)  — Log Capital per Worker",
    y = "ln(k_t+1)  — Log Next-Period Capital per Worker",
    caption = "Fitted from OLS: ln(k_t+1) = β0 + β1·ln(k_t)"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# plot 10c - r vs n+g over time, shaded where r < n+g (inefficient)

ggplot(
  data = pwt_korea |> filter(!is.na(real_interest_rate), !is.na(n)),
  aes(x = year)
) +
  geom_ribbon(
    aes(
      ymin = pmin(real_interest_rate - depreciation_rate, total_growth_rate),
      ymax = pmax(real_interest_rate - depreciation_rate, total_growth_rate),
      fill = (real_interest_rate -depreciation_rate) < total_growth_rate
    ),
    alpha = 0.2
  ) +
  geom_line(aes(y = real_interest_rate - depreciation_rate, color = "Real Interest Rate (r)"), linewidth = 1) +
  geom_line(aes(y = total_growth_rate, color = "Total Growth Rate (n + g)"), linewidth = 1) +
  scale_color_manual(
    name = "",
    values = c(
      "Real Interest Rate (r)"    = WSJ_NAVY,
      "Total Growth Rate (n + g)" = WSJ_RED
    )
  ) +
  scale_fill_manual(
    name = "",
    values = c("FALSE" = WSJ_NAVY, "TRUE" = WSJ_RED),
    labels = c(
      "FALSE" = "r ≥ n+g  (Efficient)",
      "TRUE"  = "r < n+g  (Inefficient)"
    )
  ) +
  labs(
    title    = "Dynamic Inefficiency Check: Republic of Korea",
    subtitle = "Shaded red where r < n+g — economy over-accumulates capital",
    x        = "Year",
    y        = "Rate",
    caption  = "Golden rule condition: r = n + g  |  Inefficiency threshold includes technology growth g"
  ) +
  theme_wsj() +
  theme(legend.position = LEGEND_POS)


# 11. export to excel, commented out, uncomment if you want the file

# write_xlsx(
#   pwt_korea,
#   "C:/Users/teert/OneDrive/Documents/Data/Growth_development_data/solow/pwt_kr_solow_final.xlsx"
# )

#rm(list = ls())
