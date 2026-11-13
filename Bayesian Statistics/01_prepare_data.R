# =====================================================================
# 01_prepare_data.R
# Builds the Stan data list for Egidi, Pauli & Torelli (2018).
# Mirrors "1 Data.R" + "2 Bookmakers probabilities.R" of the original code.
#
#   source("01_prepare_data.R")
#   prep <- prepare_league("PL")      # "PL", "BundesLiga", "Liga", "serieA"
# =====================================================================
library(readr)
library(skellam)
# ---- Skellam probabilities: y1 - y2 ~ PD(l1, l2) ----------------------
# probabilities (home win, draw, away win)
skellam_3way <- function(l1, l2) {
  c(win  = 1 - pskellam(0, l1, l2),     # P(y1 >  y2)
    draw = dskellam(0, l1, l2),         # P(y1 == y2)
    loss = pskellam(-1, l1, l2))        # P(y1 <  y2)
}
# <<<

# ---- odds -> probabilities -------------------------------------------
# (A) basic normalisation, eq. (2.1)
basic_prob <- function(odds) {            # odds: array [games, bookies, 3]
  inv <- 1 / odds
  sweep(inv, c(1, 2), apply(inv, c(1, 2), sum), "/")
}
# (B) Shin's procedure, eq. (2.2); z found by minimising (sum pi(z) - 1)^2
shin_one <- function(o) {                 # o: inverse odds, length 3
  b  <- sum(o)
  pz <- function(z) (sqrt(z^2 + 4 * (1 - z) * o^2 / b) - z) / (2 * (1 - z))
  z  <- optimize(function(z) (sum(pz(z)) - 1)^2, c(0, 0.5))$minimum
  pz(z)
}
shin_prob <- function(odds) {
  inv <- 1 / odds
  out <- array(NA_real_, dim(odds))
  for (n in seq_len(dim(odds)[1])) for (g in seq_len(dim(odds)[2]))
    out[n, g, ] <- shin_one(inv[n, g, ])
  out
}

# ---- implicit scoring rates, eq. (3.2) --------------------------------
# Exactly the least-squares "dist" + optim(c(1,1)) of the original code
implicit_rates <- function(prob) {        # prob: array [games, bookies, 3]
  dist <- function(par, pr) {
    l <- exp(par)
    sum((skellam_3way(l[1], l[2]) - pr)^2)
  }
  G <- dim(prob)[1]; S <- dim(prob)[2]
  th1 <- th2 <- matrix(NA_real_, G, S)
  for (n in 1:G) for (s in 1:S) {
    fit <- optim(c(1, 1), dist, pr = prob[n, s, ])
    th1[n, s] <- exp(fit$par[1]); th2[n, s] <- exp(fit$par[2])
  }
  list(theta1 = th1, theta2 = th2)
}

# ---- reading the football-data.co.uk files ----------------------------
seasons_all <- c("07-08","08-09","09-10","10-11","11-12",
                 "12-13","13-14","14-15","15-16","16-17")
# >>> ÄNDRAT (artikeln, avsnitt 4.1: Sportingbet är den femte bookmakern, inte Pinnacle)
# bookies <- c("B365", "BW", "IW", "LB", "SB_or_PS", "WH", "VC")
bookies <- c("B365", "BW", "IW", "LB", "SB", "WH", "VC")
# <<<

read_season <- function(league, s, data_dir = ".") {
  d <- suppressMessages(read_csv(file.path(data_dir, paste0(league, s, ".csv")),
                                 show_col_types = FALSE))
  d <- d[!is.na(d$HomeTeam) & !is.na(d$FTHG), ]          # empty trailing rows
  # >>> ÄNDRAT (artikeln, avsnitt 4.1: alltid Sportingbet)
  # # 5th bookmaker: Sportingbet until it disappears, then Pinnacle (as in 1 Data.R)
  # fifth <- if ("SBH" %in% names(d)) "SB" else "PS"
  # bk <- c("B365", "BW", "IW", "LB", fifth, "WH", "VC")
  bk <- bookies
  # <<<
  odds <- array(NA_real_, c(nrow(d), length(bk), 3))
  for (g in seq_along(bk)) for (k in 1:3) {
    col <- paste0(bk[g], c("H", "D", "A")[k])
    if (col %in% names(d)) odds[, g, k] <- suppressWarnings(as.numeric(d[[col]]))
  }
  list(df = data.frame(home = d$HomeTeam, away = d$AwayTeam,
                       y1 = d$FTHG, y2 = d$FTAG, ftr = d$FTR,
                       stringsAsFactors = FALSE),
       odds = odds)
}

# prepare_league <- function(league = "PL", data_dir = ".",
#                            seasons = seasons_all, jags_param = 1,
#                            compute_ppc = 1, verbose = TRUE) {
prepare_league <- function(league = "PL", data_dir = ".",
                           seasons = seasons_all, verbose = TRUE) {
  # <<<
  raw <- lapply(seasons, function(s) read_season(league, s, data_dir))
  df  <- do.call(rbind, lapply(seq_along(raw), function(i)
    cbind(raw[[i]]$df, season = i)))
  odds <- do.call(abind_1, lapply(raw, `[[`, "odds"))
  
  # Serie A: the original code drops four matches with unusable odds
  if (league == "serieA") {
    drop <- c(90, 745, 1121, 1935)
    df <- df[-drop, ]; odds <- odds[-drop, , , drop = FALSE]
  }
  
  # missing odds for a (match, bookmaker): replace by the mean of the
  # other bookmakers for that match (as in 1 Data.R)
  for (n in seq_len(dim(odds)[1])) for (k in 1:3) {
    miss <- is.na(odds[n, , k])
    if (any(miss)) odds[n, miss, k] <- mean(odds[n, !miss, k])
  }
  
  # team index: order of first appearance as home team (as in 1 Data.R)
  teams <- unique(df$home)
  T_    <- length(seasons)
  test  <- df$season == T_
  train <- !test
  
  if (verbose) message("Converting odds to probabilities ...")
  p_basic <- basic_prob(odds)
  p_shin  <- shin_prob(odds)
  if (verbose) message("Solving for implicit Skellam rates (", dim(odds)[1],
                       " x ", dim(odds)[2], " optimisations) ...")
  rates <- implicit_rates(p_basic)        # the paper's code uses basic probabilities
  
  stan_data <- list(
    M = sum(train), M_test = sum(test), nteams = length(teams), T = T_,
    S = dim(odds)[2],
    team1 = match(df$home[train], teams), team2 = match(df$away[train], teams),
    season = df$season[train],
    y1 = as.integer(df$y1[train]), y2 = as.integer(df$y2[train]),
    team1_test = match(df$home[test], teams), team2_test = match(df$away[test], teams),
    season_test = df$season[test],
    theta_bm1 = rates$theta1[train, ], theta_bm2 = rates$theta2[train, ],
    theta_bm1_test = rates$theta1[test, ], theta_bm2_test = rates$theta2[test, ],
    alpha1 = rowMeans(rates$theta1[train, ]), alpha2 = rowMeans(rates$theta2[train, ]),
    alpha1_test = rowMeans(rates$theta1[test, ]), alpha2_test = rowMeans(rates$theta2[test, ])
  )
  list(stan_data = stan_data, df = df, teams = teams, odds = odds,
       p_basic = p_basic, p_shin = p_shin, rates = rates,
       train = train, test = test, league = league)
}

abind_1 <- function(...) {                # bind 3-d arrays along dim 1
  a <- list(...)
  out <- array(NA_real_, c(sum(sapply(a, function(x) dim(x)[1])), dim(a[[1]])[2:3]))
  i <- 0
  for (x in a) { out[i + seq_len(dim(x)[1]), , ] <- x; i <- i + dim(x)[1] }
  out
}
