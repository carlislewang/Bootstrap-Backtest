setwd("~/Desktop")
library(xts)
library(readr)
library(HighFreq)
library(nortest)
library(tseries)
library(forecast)
library(boot)
#read data
XLK <- read_csv("~/Desktop/XLK.csv", col_types = cols(Date = col_date(format = "%Y-%m-%d")))
ohlcva <- xts::xts(XLK, order.by = XLK$Date)
open_ <- as.numeric(ohlcva[,2])
high_ <- as.numeric(ohlcva[,3])
low_ <- as.numeric(ohlcva[,4])
close_ <- as.numeric(ohlcva[,7])
oh_lc <- xts(open_, order.by = index(ohlcva))
oh_lc <- cbind(oh_lc, high_, low_, close_)
colnames(oh_lc) <- c("open", "high", "low", "close")
head(oh_lc)
cl_ose <- Cl(oh_lc)

#model the time series###########
log_price <- log(cl_ose)
lag_log_price <- rutils::lag_xts(log_price)
head(log_price)
head(lag_log_price)
returns <- (log_price - lag_log_price)
returns <- returns[-1]
#test if returns are stationary
adf.test(returns)
#p-value < 0.01, which means the returns are stationary
#just choose a model for the TS. Arima or simple ar(1) is no big difference
model <- ar(returns, aic=FALSE, order.max = 1)
model
beta <- model$ar[1]
model_res <- residuals(model)
model_res <- model_res[-1]

# simulate single EWMA model using historical oh_lc data
simu_ewma <- function(cl_ose, lamb_da=0.05, win_dow=50) {
  #for multi-processes
  #stopifnot("package:quantmod" %in% search() || require("quantmod", quietly=TRUE))
  
  # calculate EWMA prices
  weight_s <- exp(-lamb_da*1:win_dow)
  weight_s <- weight_s/sum(weight_s)
  ew_ma <- filter(as.numeric(cl_ose), filter=weight_s, sides=1)
  ew_ma[1:(win_dow-1)] <- cl_ose[1:(win_dow-1)]
  # determine dates right after EWMA has crossed prices
  in_dic <- xts(sign(as.numeric(cl_ose) - ew_ma), order.by=index(oh_lc))
  trade_dates <- (rutils::diff_xts(in_dic) != 0)
  trade_dates <- which(trade_dates) + 1
  trade_dates <- trade_dates[trade_dates<NROW(oh_lc)]
  # calculate positions, either: -1, 0, or 1
  po_sitions <- rep(NA_integer_, NROW(cl_ose))
  po_sitions[1] <- 0
  po_sitions[trade_dates] <- rutils::lag_xts(in_dic)[trade_dates]
  po_sitions <- xts(na.locf(po_sitions), order.by=index(oh_lc))
  op_en <- rutils::lag_xts(cl_ose)
  prices_lag <- rutils::lag_xts(cl_ose)
  position_lagged <- rutils::lag_xts(po_sitions)
  # calculate daily profits and losses
  re_turns <- po_sitions*(cl_ose - prices_lag)
  re_turns[trade_dates] <- po_sitions[trade_dates] * (cl_ose[trade_dates] - op_en[trade_dates])
  
  out_put <- cbind(po_sitions, re_turns)
  colnames(out_put) <- c("po_sitions", "re_turns")
  out_put
}  # end simu_ewma

agg_regate <- function(cl_ose, lamb_das, trend=1){
  sharpe_ratios <- sapply(lamb_das, function(lamb_da) {
    re_turns <- trend*(simu_ewma(cl_ose = cl_ose,
                          lamb_da=lamb_da)[, 2])
    # calculate annualized Sharpe ratio of strategy returns
    sqrt(260)*sum(re_turns)/sd(re_turns)/NROW(re_turns)
    
  })  # end sapply
  sharpe_ratios
} #end of agg_regate

lamb_das <- seq(0.01, 0.8, 0.01)
sharpe_ratios <- agg_regate(cl_ose, lamb_das)
plot(lamb_das, sharpe_ratios, type = "l")
ewma_trending <- simu_ewma(cl_ose, lamb_da = lamb_das[which.max(sharpe_ratios)])
sharpe_ratios <- agg_regate(cl_ose, lamb_das, trend = -1)
plot(lamb_das, sharpe_ratios, type = "l")
ewma_revert <- -simu_ewma(cl_ose, lamb_da = lamb_das[which.max(sharpe_ratios)])

#trending
po_sitions <- ewma_trending[,1]
pn_l <- cumsum(ewma_trending[, 2])
#combine
po_sitions <- (ewma_trending[,1]+ewma_revert[,1])/2
returns <- (ewma_trending[,2]+ewma_revert[,2])/2
pn_l <- cumsum(returns)
#mean-reverting
po_sitions <- ewma_revert[,1]
daily_return <- ewma_revert[,2]
pn_l <- cumsum(daily_return)
pn_l <- cbind(close_ - as.numeric(close_[1]),
              pn_l)
colnames(pn_l) <- c("XLK", "EWMA PnL")
plot_theme <- chart_theme()
plot_theme$col$line.col <- c("orange", "blue")
chart_Series(pn_l, theme=plot_theme,
             name="Performance of EWMA Strategy")
add_TA(po_sitions ==1, on=-1,
       col="lightgreen", border="lightgreen")
legend("topleft", legend=colnames(pn_l),
       inset=0.05, bg="white", lty=c(1, 1), lwd=c(2,2), col=plot_theme$col$line.col, bty="n")

#draw pnl distribution
plot(density(daily_return))
returns <- pn_l[,2]
returns <- returns[diff_xts(po_sitions) != 0]
returns_lag <- lag_xts(returns)
returns <- returns - returns_lag
trade_close <- close_[diff_xts(po_sitions) != 0]
trade_close <- xts(trade_close, order.by = index(returns))
close_lag <- lag_xts(trade_close)
returns_round_trade <- as.numeric(returns/close_lag)
hist(returns_round_trade, xlab = "returns percentage")
plot(density(returns_round_trade))

#bootstrapping##############
bt_statistics <- function(beta, data, indices, cl_ose, model_mean, lamb_das){
  tmp_res <- data[indices]
  price <- array(as.numeric(cl_ose[1]))
  price <- append(price, as.numeric((cl_ose[2])))
  re_turn <- as.numeric(log(cl_ose[2])) - as.numeric(log(cl_ose[1]))
  for(i in tmp_res){
    re_turn <- beta*(re_turn - model_mean) + i + model_mean
    price <- append(price, price[length(price)] * exp(re_turn))
  }
  price <- xts::xts(price, order.by = index(cl_ose))
  #trading stategy is choosing the biggest sharpe ratio
  sharpe_ratios <- agg_regate(price, lamb_das, trend = -1)
  ewma_revert <- -simu_ewma(price, lamb_da = lamb_das[which.max(sharpe_ratios)])
  #get results of trading
  po_sitions <- ewma_revert[,1]
  daily_return <- ewma_revert[,2]
  pn_l <- cumsum(daily_return)
  #criterias####
  #trade_num
  trade_num <- sum(sign(abs(diff_xts(po_sitions))))
  #total return
  total_return <- as.numeric(pn_l[nrow(pn_l)])/as.numeric(close_[1])
  #winning rate
  returns <- pn_l[diff_xts(po_sitions) != 0]
  returns_lag <- lag_xts(returns)
  winning <- sum((returns - returns_lag) > 0)
  winning_rate <- winning / trade_num
  #max drawdwon
  drawdown <- (pn_l - cummax(pn_l))/as.numeric(cl_ose[match(cummax(pn_l), pn_l)])
  max_drawdown <- min(drawdown, na.rm = TRUE)
  #sharpe ratio
  sharpe_ratio <- sqrt(260)*sum(daily_return)/sd(daily_return)/NROW(daily_return)
  
  return(cbind(tradeNum = trade_num, total_return = total_return, winning_rate = winning_rate, 
               max_drawdown = max_drawdown, sharpe_ratio = sharpe_ratio))
}

result <- boot(data = model_res, statistic = bt_statistics, R=1000, beta=beta, cl_ose=cl_ose, model_mean=model$x.mean, lamb_das = lamb_das)
result$t0
data <- result$t
data <- as.data.frame(data)
colnames(data) <- colnames(result$t0)
apply(data, 2, mean)
data <- rbind(result$t0, data, apply(data, 2, mean))
write_csv(data, path = "bootstrap_Result.csv")
