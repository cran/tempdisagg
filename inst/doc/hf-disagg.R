## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width  = 7.1,
  fig.height = 3.5,
  comment = "#>"
)


## -----------------------------------------------------------------------------
library(tempdisagg)
data(tempdisagg)
head(gdp.q)

## -----------------------------------------------------------------------------
library(tsbox)
ts_plot(gdp.q, title = "Swiss GDP", subtitle = "real, not seasonally adjusted")

## -----------------------------------------------------------------------------
m.d.noind <- td(gdp.q ~ 1, to = "daily", method = "fast")
summary(m.d.noind)

## -----------------------------------------------------------------------------
gdp.d.noind <- predict(m.d.noind)
ts_plot(
  ts_scale(
    ts_c(gdp.d.noind, gdp.q)
  ),
  title = "Daily disaggregated GDP",
  subtitle = "no indicator"
)

## -----------------------------------------------------------------------------
all.equal(ts_frequency(gdp.d.noind, "quarter", aggregate = "sum"), gdp.q)

## -----------------------------------------------------------------------------
ts_plot(spi.d, title = "Swiss Performance Index", subtitle = "daily values, interpolated")

## -----------------------------------------------------------------------------
m.d.stocks <- td(gdp.q ~ spi.d, method = "chow-lin-fixed", fixed.rho = 0.9)
summary(m.d.stocks)

## -----------------------------------------------------------------------------
gdp.d.stocks <- predict(m.d.stocks)
ts_plot(
  ts_scale(
    ts_c(gdp.d.stocks, gdp.q)
  ),
  title = "Daily disaggregated GDP",
  subtitle = "one indicator"
)

