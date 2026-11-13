# =====================================================================
# 03_compare_with_paper.R
# Figures 2, 3, 4(a), 6, 7, 8 and Tables 2-3 (EPL) of
# Egidi, Pauli & Torelli (2018), computed from our Stan fit.
#
# Every quantity is defined the same way as in the authors' original
# scripts ("4 Results.R", "5 Prediction Betting Power.R"), so the numbers
# can be compared directly with the paper. Where the original code makes
# a choice that is not described in the paper, this is noted with
# "AUTHORS' DEFINITION".
# =====================================================================
library(rstan)
library(ggplot2)
library(skellam)
source("01_prepare_data.R")

league <- "PL"
prep   <- readRDS(paste0("prep_", league, ".rds"))
fit    <- readRDS(paste0("fit_", league, ".rds"))
d      <- prep$stan_data
teams  <- prep$teams
n_tr   <- length(d$y1)            # training matches
n_te   <- length(d$team1_test)    # test matches (2016-17)

# colours used in the paper (bayesplot schemes)
c_red   <- "#A25050"   # attack
c_blue  <- "#005B96"   # defence
c_grey  <- "#BFBFBF"   # ribbons
c_red2  <- "#E8312F"   # Fig 8, red
c_blue2 <- "#3A53A4"   # Fig 4(a) and Fig 8, blue

# =====================================================================
# Figure 2: attack / defence effects, 2016-17 teams (median, 50% bars)
# =====================================================================
att <- rstan::extract(fit, "att")$att          # draws x nteams x T
def <- rstan::extract(fit, "def")$def
teams_last <- unique(prep$df$home[prep$df$season == d$T])

effect_summary <- function(arr, label) {
  do.call(rbind, lapply(teams_last, function(tn) {
    x <- arr[, match(tn, teams), ]
    data.frame(team = tn, season = seq_len(d$T), effect = label,
               lo  = apply(x, 2, quantile, 0.25),
               mid = apply(x, 2, median),
               hi  = apply(x, 2, quantile, 0.75))
  }))
}
fig2 <- rbind(effect_summary(att, "Attack"), effect_summary(def, "Defence"))

g2 <- ggplot(fig2, aes(season, mid, group = effect)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = c_grey) +
  geom_line(aes(colour = effect), linewidth = 1, show.legend = FALSE) +
  facet_wrap(~ team) +
  scale_colour_manual(values = c(Attack = c_red, Defence = c_blue)) +
  scale_x_continuous(breaks = c(1, 5, 10), labels = c("07/08", "11/12", "16/17")) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  labs(x = "Seasons", y = "Teams' effects",
       title = "Attack and defense effects (50% posterior bars)",
       subtitle = "for teams of Premier League 2016/2017") +
  theme_classic()
print(g2)
ggsave("fig2.png", g2, width = 5.5, height = 4.2, dpi = 300)

# =====================================================================
# Figure 3: mixing parameter p_m (50% bars)
# AUTHORS' DEFINITION: the 25%, 50% and 75% quantiles are each sorted
# separately (not ordered by the median).
# =====================================================================
p_draws <- as.matrix(fit, pars = "p")
fig3 <- data.frame(game = seq_len(ncol(p_draws)),
                   lo   = sort(apply(p_draws, 2, quantile, 0.25)),
                   mid  = sort(apply(p_draws, 2, median)),
                   hi   = sort(apply(p_draws, 2, quantile, 0.75)))
cat("\nFig 3: posterior median of p_m over matches\n"); print(summary(fig3$mid))

g3 <- ggplot(fig3, aes(game, mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = c_grey) +
  geom_line(colour = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Games", y = "p", title = " ") +
  theme_classic()
print(g3)

# =====================================================================
# Figure 4(a): share of matches with goal difference k, per draw
# =====================================================================
y1r <- rstan::extract(fit, "y1_rep")$y1_rep    # draws x n_tr
y2r <- rstan::extract(fit, "y2_rep")$y2_rep
gd_rep <- y1r - y2r
gd_obs <- d$y1 - d$y2

kk      <- -3:3
share_k <- sapply(kk, function(k) rowMeans(gd_rep == k))     # draws x 7
obs_k   <- sapply(kk, function(k) mean(gd_obs == k))
print(rbind(k = kk, observed = round(obs_k, 3), replicated = round(colMeans(share_k), 3)))

fig4a_rep <- data.frame(k = rep(kk, each = nrow(share_k)), prob = as.vector(share_k))
fig4a_obs <- data.frame(k = kk, prob = obs_k)


g4a <- ggplot() +
  geom_jitter(data = fig4a_rep, aes(k, prob), width = 0.35, height = 0,
              colour = "grey70", alpha = 0.2, size = 1) +
  geom_segment(data = fig4a_obs, aes(x = k - 0.5, xend = k + 0.5, y = prob, yend = prob),
               colour = c_blue2, linewidth = 3, lineend = "butt") +
  scale_x_continuous(breaks = kk) +
  labs(x = "Goal difference", y = "Posterior probabilities") +
  theme_classic()
print(g4a)
g4a <- g4a + theme(aspect.ratio = 0.74)

ggsave("fig4a.png", g4a, width = 5.5, height = 4.2, dpi = 300)

# =====================================================================
# Simulated season 2016-17: Tables 2-3, Figures 6-7
# =====================================================================
y1t <- rstan::extract(fit, "y1_test")$y1_test    # draws x n_te
y2t <- rstan::extract(fit, "y2_test")$y2_test
nd  <- nrow(y1t)
tm  <- sort(unique(c(d$team1_test, d$team2_test)))

# points per team and draw
pts <- matrix(0, nd, length(tm))
for (m in seq_len(n_te)) {
  h <- match(d$team1_test[m], tm); a <- match(d$team2_test[m], tm)
  hw <- y1t[, m] >  y2t[, m]
  dr <- y1t[, m] == y2t[, m]
  aw <- y1t[, m] <  y2t[, m]
  pts[, h] <- pts[, h] + 3 * hw + dr
  pts[, a] <- pts[, a] + 3 * aw + dr
}

# AUTHORS' DEFINITION: ranking on points only (no goal difference);
# ties go to the team with the lower team index.
rank_draw <- t(apply(pts, 1, function(x) rank(-x, ties.method = "first")))
nt        <- length(tm)
rank_prob <- sapply(seq_len(nt), function(r) colMeans(rank_draw == r))
rownames(rank_prob) <- teams[tm]

# observed final table
obs <- prep$df[prep$df$season == d$T, ]
obs_pts <- sapply(teams[tm], function(tn)
  sum(3 * (obs$home == tn & obs$y1 > obs$y2) + 3 * (obs$away == tn & obs$y2 > obs$y1) +
        (obs$home == tn | obs$away == tn) * (obs$y1 == obs$y2)))
obs_rank <- rank(-obs_pts, ties.method = "first")

tab2 <- data.frame(P1 = rank_prob[, 1], P2 = rank_prob[, 2], P3 = rank_prob[, 3],
                   rank = obs_rank, points = obs_pts)
tab3 <- data.frame(P18 = rank_prob[, nt - 2], P19 = rank_prob[, nt - 1], P20 = rank_prob[, nt],
                   rank = obs_rank, points = obs_pts)
cat("\nTable 2 (top 3)\n");     print(round(tab2[order(tab2$rank), ][1:3, ], 4))
cat("\nTable 3 (relegated)\n"); print(round(tab3[order(tab3$rank), ][(nt - 2):nt, ], 4))



# ---- Figure 6: exact-score heatmap for one test match ------------------
# AUTHORS' DEFINITION: P(home = h, away = a) is the product of the two
# marginal predictive probabilities, not the joint frequency.
home_team <- "Hull"
away_team <- "Middlesbrough"
m6 <- which(teams[d$team1_test] == home_team & teams[d$team2_test] == away_team)

p_home <- sapply(0:5, function(h) mean(y1t[, m6] == h))
p_away <- sapply(0:5, function(a) mean(y2t[, m6] == a))
fig6 <- expand.grid(Home = 0:5, Away = 0:5)
fig6$Prob <- as.vector(outer(p_home, p_away))

res6 <- obs[obs$home == home_team & obs$away == away_team, ]
g6 <- ggplot(fig6, aes(Home, Away, fill = Prob)) +
  geom_tile() +
  geom_tile(data = data.frame(Home = res6$y1, Away = res6$y2), aes(Home, Away),
            inherit.aes = FALSE, fill = NA, colour = "red", linewidth = 1.5) +
  scale_fill_gradient(low = "white", high = "black", name = "Prob",
                      limits = c(0, 0.19),
                      breaks = c(0.05, 0.10, 0.15),
                      oob = scales::squish,
                      guide = guide_colourbar(nbin = 20, raster = FALSE)) +
  scale_x_continuous(breaks = seq(0, 5, 2)) +
  scale_y_continuous(breaks = seq(0, 5, 2)) +
  coord_equal() +
  labs(title = paste(home_team, "  ", away_team)) +
  theme_bw() +
  theme(panel.grid = element_blank())
print(g6)

# ---- Figure 7: predicted final points (median, 50% bars) --------------
# AUTHORS' DEFINITION: teams sorted by predicted median (ties keep the
# order in which teams first appear as home team in the test season).
q7 <- apply(pts, 2, quantile, c(0.25, 0.5, 0.75))
fig7 <- data.frame(team = teams[tm], lo = q7[1, ], mid = q7[2, ], hi = q7[3, ], obs = obs_pts)
appear  <- teams[unique(d$team1_test)]
fig7    <- fig7[match(appear, fig7$team), ]
fig7    <- fig7[order(-fig7$mid), ]
fig7$team <- factor(fig7$team, levels = fig7$team)
inside50 <- mean(fig7$obs >= fig7$lo & fig7$obs <= fig7$hi)
cat(sprintf("\nFig 7: share of observed points inside 50%% interval: %.2f\n", inside50))

g7 <- ggplot(fig7, aes(team, mid, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = c_grey) +
  geom_line() +
  labs(x = "Teams", y = "Points") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g7)

# =====================================================================
# Figure 8: expected profit per bookmaker
# Model probabilities: Skellam at the posterior mean of the mixed rates.
# AUTHORS' DEFINITION:
#   Strategy A: one unit on argmax p*odds, profit = max(p*odds) - 1
#   Strategy B: stake = p*odds / (2 * 1.2) (Rue & Salvesen, variance 1.2).
#     The stake for every bookmaker is taken at the outcome chosen for the
#     LAST bookmaker, exactly as in the original loop. Without this, B
#     would be identical to A, since the stake cancels in the profit.
#   Mean and sd are taken over the test matches.
# =====================================================================
th1 <- colMeans(rstan::extract(fit, "theta1_test")$theta1_test)
th2 <- colMeans(rstan::extract(fit, "theta2_test")$theta2_test)
P_model <- t(mapply(skellam_3way, th1, th2))          # n_te x 3 (H, D, A)

profits_AB <- function(P, O) {                        # P, O: matches x bookies x 3
  po    <- P * O
  stake <- po / (2 * 1.2)
  gain  <- stake * po
  nb_   <- dim(O)[2]
  
  prof_A <- apply(po, c(1, 2), max) - 1
  
  best_B <- apply(gain, c(1, 2), max)
  pick_B <- apply(gain, c(1, 2), which.max)
  paid_B <- t(sapply(seq_len(dim(O)[1]), function(m) stake[m, , pick_B[m, nb_]]))
  prof_B <- (best_B - paid_B) / paid_B
  
  list(A = data.frame(est = colMeans(prof_A), sd = apply(prof_A, 2, sd)),
       B = data.frame(est = colMeans(prof_B), sd = apply(prof_B, 2, sd)))
}

odds_te <- prep$odds[prep$test, , , drop = FALSE]
n_bk    <- dim(odds_te)[2]
P_model_arr <- aperm(array(P_model, c(n_te, 3, n_bk)), c(1, 3, 2))   # same probs for every bookie

pr_mod <- profits_AB(P_model_arr, odds_te)
pr_bas <- profits_AB(prep$p_basic[prep$test, , ], odds_te)
pr_shi <- profits_AB(prep$p_shin[prep$test, , ],  odds_te)

bk_lab   <- c(B365 = "Bet365", BW = "Bwin", IW = "Interw.", LB = "L.brokes",
              SB = "S.bet", VC = "VC Bet", WH = "W.Hill")
bk_order <- c("Bet365", "Bwin", "Interw.", "L.brokes", "S.bet", "VC Bet", "W.Hill")
bk_names <- unname(bk_lab[bookies])



# ---- plot in the style of the paper (base R) --------------------------
pick <- function(s) { x <- prof[prof$series == s, ]; x[match(bk_order, x$bookie), ] }
sA <- pick("Str. A"); sB <- pick("Str. B")
sBas <- pick("Basic Str. B"); sShin <- pick("Shin Str. B")
nb <- length(bk_order)
yy <- nb:1                                     # Bet365 on top
hbar <- function(s, y, col) segments(s$est - s$sd, y, s$est + s$sd, y, col = col)

op <- par(mar = c(5, 7, 5, 2))
plot(NA, xlim = c(-0.5, 1.0), ylim = c(0.4, nb + 0.4), axes = FALSE, xlab = "", ylab = "")
axis(3, at = c(-0.5, 0, 0.5, 1.0), labels = c("\u22120.5", "0.0", "0.5", "1.0"))
mtext("Expected profit (estim. +/\u2212 1 sd)", side = 3, line = 2.5)
axis(2, at = yy, labels = bk_order, las = 1, lwd = 0)
abline(v = 0, lty = 3)

hbar(sB, yy + 0.10, c_red2);    points(sB$est, yy + 0.10, pch = 15, col = c_red2, cex = 1.2)
hbar(sA, yy - 0.10, c_blue2);   points(sA$est, yy - 0.10, pch = 15, col = c_blue2)
hbar(sShin, yy - 0.25, c_red2)                                  # Shin: interval only
hbar(sBas, yy - 0.42, c_red2);  points(sBas$est, yy - 0.42, pch = 17, col = c_red2)

legend(x = 0.62, y = nb - 0.6, cex = 0.8,
       legend = c("Str. A", "Str. B", "Basic Str. B", "Shin Str. B"),
       pch = c(15, 15, 17, NA), col = c(c_blue2, c_red2, c_red2, NA))
mtext("Premier League", side = 1, line = 2)
par(op)




# =====================================================================
# Simulated season 2016-17: Table 2 with MCSE
# =====================================================================
library(posterior)

arr1   <- as.array(fit, pars = "y1_test")             # iter x chains x n_te
arr2   <- as.array(fit, pars = "y2_test")
n_iter <- dim(arr1)[1]; n_ch <- dim(arr1)[2]
y1t <- matrix(arr1, n_iter * n_ch, dim(arr1)[3])      # draws x n_te
y2t <- matrix(arr2, n_iter * n_ch, dim(arr2)[3])
nd  <- nrow(y1t)
tm  <- sort(unique(c(d$team1_test, d$team2_test)))

# points per team and draw
pts <- matrix(0, nd, length(tm))
for (m in seq_len(n_te)) {
  h <- match(d$team1_test[m], tm); a <- match(d$team2_test[m], tm)
  hw <- y1t[, m] >  y2t[, m]
  dr <- y1t[, m] == y2t[, m]
  aw <- y1t[, m] <  y2t[, m]
  pts[, h] <- pts[, h] + 3 * hw + dr
  pts[, a] <- pts[, a] + 3 * aw + dr
}

# AUTHORS' DEFINITION: ranking on points only (no goal difference);
# ties go to the team with the lower team index.
rank_draw <- t(apply(pts, 1, function(x) rank(-x, ties.method = "first")))
nt        <- length(tm)
rank_prob <- sapply(seq_len(nt), function(r) colMeans(rank_draw == r))
rownames(rank_prob) <- teams[tm]

# MCSE for P(team j in place r), accounting for autocorrelation
rank_mcse <- sapply(1:3, function(r) sapply(seq_len(nt), function(j)
  mcse_mean(matrix(as.numeric(rank_draw[, j] == r), n_iter, n_ch))))
rownames(rank_mcse) <- teams[tm]

# observed final table
obs <- prep$df[prep$df$season == d$T, ]
obs_pts <- sapply(teams[tm], function(tn)
  sum(3 * (obs$home == tn & obs$y1 > obs$y2) + 3 * (obs$away == tn & obs$y2 > obs$y1) +
        (obs$home == tn | obs$away == tn) * (obs$y1 == obs$y2)))
obs_rank <- rank(-obs_pts, ties.method = "first")

tab2 <- data.frame(P1 = sprintf("%.2f (%.3f)", rank_prob[, 1], rank_mcse[, 1]),
                   P2 = sprintf("%.2f (%.3f)", rank_prob[, 2], rank_mcse[, 2]),
                   P3 = sprintf("%.2f (%.3f)", rank_prob[, 3], rank_mcse[, 3]),
                   rank = obs_rank, points = obs_pts,
                   row.names = teams[tm])
cat("\nTable 2 (top 3), MCSE in parentheses\n")
print(tab2[order(tab2$rank), ][1:3, ])
tab2 <- data.frame(P1 = sprintf("%.4f (%.4f)", rank_prob[, 1], rank_mcse[, 1]),
                   P2 = sprintf("%.4f (%.4f)", rank_prob[, 2], rank_mcse[, 2]),
                   P3 = sprintf("%.4f (%.4f)", rank_prob[, 3], rank_mcse[, 3]),
                   rank = obs_rank, points = obs_pts,
                   row.names = teams[tm])
cat("\nTable 2 (top 3), MCSE in parentheses\n")
print(tab2[order(tab2$rank), ][1:3, ])



draws <- as_draws_df(fit)

summarise_draws(
  subset_draws(draws, variable = c("mu", "home", "mu_att", "mu_def",
                                   "sigma2_att", "sigma2_def", "tau1", "tau2")),
  sd, mcse_sd
)
