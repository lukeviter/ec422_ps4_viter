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




