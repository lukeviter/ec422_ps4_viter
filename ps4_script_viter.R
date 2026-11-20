library(tidyverse)
library(fpp3)
library(BVAR)

#2.1a
liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv",
  na = ".",
  col_types = cols(
    observation_date = col_date(format = "%Y-%m-%d"),
    MRTSSM4453USN = col_double()
  )
)

liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month <= yearmonth("2019 Dec"))

has_gaps(liquor)

#2.1b
liquor |>
  autoplot(sales) +
  labs(
    title = "Monthly US liquor sales",
    x = "Year",
    y = "Millions of dollars"
  )

liquor <- liquor |>
  mutate(log_sales = log(sales))

liquor |>
  gg_season(log_sales) +
  labs(
    title = "Seasonal plot: Monthly US liquor sales",
    x = "Month",
    y = "Log of millions of dollars"
  )

# There is a definite upward trend to the data, with the size of the monthly swing increasing steadily as time progresses. This suggests that logs could be a good tool to stabilize the data. In the seasonal plot, the summer months as well as December seem to be peaks, while winter and fall are generally lower.

#2.2a
liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)),
        ets   = ETS(log(sales)))

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

# The ARIMA model selected is an ARIMA(1, 0, 2)(2,1, 1)[12] model, meaning the series is differenced once at the seasonal interval to try and achieve stationarity, but not at the regular interval.
# The ETS model selected is an ETS(A, A, A) model, meaning it has additive error, additive trend, and additive seasonality.

#2.2b
liquor_fit |> select(arima) |> gg_tsresiduals()

liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

# The residuals do not look like white noise, and we reject the null hypothesis of no autocorrelation in the Ljung-Box test.

#2.3a
liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(sales),
        ets    = ETS(log(sales)),
        arima  = ARIMA(log(sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")

# 97 forecast origins

#2.3b
liquor_cv_fc <- liquor_cv_fit |>
  forecast(h = 1)

liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

liquor_cv |>
  as_tibble() |>
  group_by(.id) |>
  summarise(origin = max(Month)) |>
  summarise(
    first_origin = min(origin),
    last_origin = max(origin)
  )

# The first forecast is December of 2001, and the last forecast is November of 2019. According to the accuracy report, the ARIMA model forecasts more accurately.

#2.4
liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs$arima, liquor_errs$snaive, h = 1, power = 2)

# ARIMA vs ETS, DM = -2.9062, p-value = 0.004042
# ARIMA vs snaive, DM = -10.755, p-value < 2.2e-16
# We can conclude that the ARIMA model has the lowest forecast errors of the three models by a good margin, especially relative to the snaive model. If I were the liquor store, I would use the ARIMA model.

#3a
liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)

liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

# Combo ranks above ETS and snaive but below ARIMA.

#3b
cor(liquor_errs$arima, liquor_errs$ets)

# If combining two models is most efficient when their errors are structurally different, then it makes sense that the combination did not beat the best model in this context. ARIMA and ETS have a correlation of 0.89, which is quite high, meaning they are structurally similar and so combining them does not yield a better model. Moreover, we know that ARIMA has less estimated loss on average than ETS, so combining ETS with ARIMA is more likely to add loss than be efficient.

#4.1a
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2))

#4.1b
growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),     # real GDP
            cons   = growth(PCECC96),   # real consumption
            inv    = growth(GPDIC1),    # real investment
            emp    = growth(PAYEMS),    # payroll employment
            infl   = growth(GDPCTPI),   # GDP price index
            unrate = difference(UNRATE),     # unemployment rate
            ffr    = difference(FEDFUNDS),   # federal funds rate
            gs10   = difference(GS10)) |>    # 10-year Treasury yield
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()

nrow(y)
head(y)

#4.1c
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

# 240 - 120 = 120 forecast origins of each variable
# The first forecast will be 1990 Q1, and the last forecast will be 2019 Q4

#4.2
var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)            # columns: y_t, then its lags 1 to p
  M <- ncol(z)                    # number of variables
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))   # constant, latest p lags
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)

# The VAR(1) model has a lower RMSE for all of the variables except the inflation rate. It is possible that the VAR(1) forecasts better because the additional paramters of the VAR(4) model may lead to data overfitting and worse forecasts when taken out-of-sample.

#4.3a
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))         # constant and 4 lags of each variable
    fc <- rbind(fc, x %*% coef(fit_T))   # forecast of y_{T+1}
  }
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

#4.3b
set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)

# The BVAR(4) model with lambda = 0.3 beats the VAR(4) model for all variables. This difference results from the different criterion by which each model is estimated. VAR(4) simply minimizes OLS, whereas BVAR(4) introduces a prior that structurally alters the effect of lags, other variables' lags, and distant lags, resulting in a different RMSE value.
# Compared to the VAR(1) model, BVAR(4) with a lambda 0.3 performs worse in every variable except for the Inflation Rate and the GS10.
# When lambda = 0.05, we know that there will be much stronger shrinkage and the model will place heavier weight on its own prior assumptions. Moreover, choosing a model based off of this sample is not a guarantee that it will perform well once it is out-of-sample, as selecting the model based on a low in-sample RMSE does not necessarily mean that it will have that same RMSE out-of-sample, as that could have just been a fluke.

#4.3c
e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

e_var4 <- actual[, "infl"] - fc_var4[, 1]
e_bvar <- actual[, "infl"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

# Based on these tests, we see that the BVAR(4) model performs better than the VAR(4) model for GDP and inflation by relatively significant amounts, despite the RMSE values of both of these models being relatively close for these variables. This suggests that RMSE is not the sole indicator of a model's forecasting ability, and that a wide range of components should be considered. It is more nuanced than 4.3b would suggest on its own.

#5.1a
md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,                      # level
         "2" = d(x),                   # first difference
         "3" = d(d(x)),                # second difference
         "4" = log(x),                 # log
         "5" = 100 * d(log(x)),        # growth rate
         "6" = 100 * d(d(log(x))),     # change in the growth rate
         "7" = d(x / dplyr::lag(x) - 1))
}

#5.1b
fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

fm <- fm[, colSums(is.na(fm)) == 0]

dim(fm)

#5.1c
tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")]

# INDPRO and PAYEMS receive a growth rate transformations, CPIAUCSL receeives a change in growth rate transformation, and FEDFUNDS receives a first difference transformation. The transformations of INDPRO and FEDFUNDS are likely different because they measure different economic components and therefor undergo different structural changes over time, requiring distinct statistical tools to stabilize them.

#5.2a
Z <- fm |> select(-date) |> as.matrix()

pc <- prcomp(scale(Z))

summary(pc)$importance[, 1:5]

# PC1 explains 0.1498 of the variance, PC2 explains 0.07407 of the variance, and PC3 explains 0.0706 of the variance.

#5.2b
F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))

tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

#5.2c
sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]

# Based on the correlations, I would assume that the first factor measuers overall economic activity, such as employment or industrial production. The strongest correlations related to the second factor, however, make me believe it is more related to prices and inflation, based on the FRED-MD index.





