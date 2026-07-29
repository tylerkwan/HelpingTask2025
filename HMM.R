## ================================================================
##  HMM Code 
## ================================================================

library(dplyr)
library(depmixS4)
library(ggplot2) 

data_dir <- "."

set.seed(2026)

## ================================================================
## SECTION 1: Loading the data
## ================================================================

data <- read.csv("dataset.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  rename(
    pid              = PID,
    task_design      = 'Task Design',
    trial            = Trial,
    block            = Block,
    policy           = 'Helper Policy',
    effort_intended  = 'Effort Intended',
    effort_actual    = 'Effort Actual',
    partner_effort   = 'Helper Effort',
    satisfaction     = 'Satisfaction Rating',
    trust            = 'Trust Rating',
    points_helper    = 'Points For Helper',
    points_seeker    = 'Points For Help-Seeker',
    inequality = 'Points Inequality',
    pe_self = 'Self Prediction Error')

## ================================================================
## SECTION 2: Realign Trials
## ================================================================

data_adj <- data %>% group_by(pid, block) %>% mutate(
  effort_selection = lead(effort_intended),
  block_trial = row_number(),
  is_closing_row = row_number() == n()
  ) %>% ungroup() %>%
  filter(!is_closing_row) %>%    # drops the closing-evaluation row per block
  arrange(pid, block, block_trial)

z <- function(x) as.numeric(scale(x))

data_adj <- data_adj %>% mutate(
  sat_z            = z(satisfaction),
  trust_z          = z(trust),
  effort_z         = z(effort_selection),
  effort_actual_z  = z(effort_actual),
  partner_effort_z = z(partner_effort),
  points_seeker_z  = z(points_seeker),
  points_helper_z  = z(points_helper),
  pe_z             = z(pe_self),
  inequality_z     = z(inequality),
  points_seeker_rank_z = z(rank(points_seeker, ties.method = "average", na.last = "keep"))
)

cov_cols <- c("partner_effort_z","points_seeker_z","points_helper_z","pe_z","inequality_z","points_seeker_rank_z")
data_adj$covariate_imputed <- Reduce('|', lapply(data_adj[cov_cols], is.na))
cat("Rows with an imputed covariate (from 1 trial with missing effort_actual/points):", sum(data_adj$covariate_imputed), "\n")
for (cc in cov_cols) data_adj[[cc]][is.na(data_adj[[cc]])] <- 0


print(round(cor(data_adj[, c("points_seeker_z","partner_effort_z","points_helper_z","pe_z","inequality_z")]), 2))
ntimes_block <- data_adj %>% count(pid, block) %>% pull(n)

## ================================================================
## SECTION 3: Number of States
## ================================================================

fit_k_multistart <- function(data, ntimes, k, n_starts = 15, sd_floor = 0.15) {
  candidates <- list()
  for (s in seq_len(n_starts)) {
    set.seed(7000 * k + s)
    mod <- depmix(list(effort_z ~ 1, sat_z ~ 1, trust_z ~ 1), data = data,
                    nstates = k, ntimes = ntimes,
                    family = list(gaussian(), gaussian(), gaussian()))
    f <- tryCatch(fit(mod, verbose = FALSE, emcontrol = em.control(maxit = 500, random.start = TRUE)),
                  error = function(e) NULL)
    if (is.null(f)) next
    sds <- unlist(lapply(f@response, function(state) sapply(state, function(r) r@parameters$sd)))
    if (any(is.na(sds))) next
    candidates[[length(candidates) + 1]] <- list(fit = f, min_sd = min(sds), bic = BIC(f))
  }
  valid <- Filter(function(c) c$min_sd >= sd_floor, candidates)
  pool  <- if (length(valid) > 0) valid else candidates
  if (length(pool) == 0) return(NULL)
  best <- pool[[which.min(sapply(pool, function(c) c$bic))]]
  list(fit = best$fit, n_valid = length(valid), n_total = length(candidates))
}

model_diagnostics <- function(fit, data) {
  vit <- posterior(fit, type = "viterbi")
  post_cols <- grep("^S[0-9]+$", names(vit))
  list(certainty = mean(apply(vit[, post_cols, drop = FALSE], 1, max)),
       fewest_subj = min(tapply(data$pid, vit$state, function(s) length(unique(s)))))
}

fits_store <- list(); selection <- data.frame()
for (k in 2:5) {
  r <- fit_k_multistart(data_adj, ntimes_block, k)
  if (is.null(r)) { cat("  k=", k, ": no fit converged\n"); next }
  diag <- model_diagnostics(r$fit, data_adj)
  selection <- rbind(selection, data.frame(
    k = k, non_degenerate_restarts = r$n_valid, BIC = round(BIC(r$fit), 1),
    fewest_subjects_in_any_state = diag$fewest_subj, mean_state_certainty = round(diag$certainty, 3)))
  fits_store[[as.character(k)]] <- r$fit
}

selection$reliable <- selection$non_degenerate_restarts >= 3 & selection$fewest_subjects_in_any_state >= 2
best_k <- selection$k[selection$reliable][which.min(selection$BIC[selection$reliable])]

print(selection)
print(best_k)

## ================================================================
## SECTION 4: IO HMM Models
## ================================================================

fit_transition_model <- function(trans_formula, seed_base, n_starts = 10, sd_floor = 0.15) {
  candidates <- list()
  for (s in seq_len(n_starts)) {
    set.seed(seed_base + s)
    mod <- depmix(list(effort_z ~ 1, sat_z ~ 1, trust_z ~ 1),
                    transition = trans_formula, data = data_adj, nstates = best_k, ntimes = ntimes_block,
                    family = list(gaussian(), gaussian(), gaussian()))
    f <- tryCatch(fit(mod, verbose = FALSE, emcontrol = em.control(maxit = 500, random.start = TRUE)),
                  error = function(e) NULL)
    if (is.null(f)) next
    sds <- unlist(lapply(f@response, function(state) sapply(state, function(r) r@parameters$sd)))
    if (any(is.na(sds))) next
    candidates[[length(candidates)+1]] <- list(fit=f, min_sd=min(sds), bic=BIC(f), ll=as.numeric(logLik(f)))
  }
  valid <- Filter(function(c) c$min_sd >= sd_floor, candidates)
  pool <- if (length(valid) > 0) valid else candidates
  if (length(pool) == 0) return(NULL)
  best <- pool[[which.min(sapply(pool, function(c) c$bic))]]
  list(fit=best$fit, n_valid=length(valid), n_total=length(candidates),
       bic=best$bic, ll=best$ll, npar=attr(logLik(best$fit), "df"))
}

specs <- list(
  "M0_null"                       = ~ 1,
  "M1_effort_z"                   = ~ partner_effort_z,
  "M2_returns_z"                  = ~ points_seeker_z,
  "M3_returns_rank_z"             = ~ points_seeker_rank_z,
  "M4_ineq_z"                     = ~ inequality_z,
  "M5_effort_and_ineq_z"          = ~ partner_effort_z * inequality_z
)

results <- list()
for (nm in names(specs)) {
  cat("  fitting", nm, "...\n")
  results[[nm]] <- fit_transition_model(specs[[nm]], seed_base = 10000 * nchar(nm))
}

comp <- data.frame(
  model = names(results),
  n_valid = sapply(results, function(r) if (is.null(r)) NA else r$n_valid),
  logLik  = sapply(results, function(r) if (is.null(r)) NA else round(r$ll, 1)),
  npar    = sapply(results, function(r) if (is.null(r)) NA else r$npar),
  BIC     = sapply(results, function(r) if (is.null(r)) NA else round(r$bic, 1))
)

comp$reliable <- !is.na(comp$n_valid) & comp$n_valid >= 3
comp <- comp[order(comp$BIC), ]

print(comp, row.names = FALSE)

ll0 <- results[["M0_null"]]$ll; npar0 <- results[["M0_null"]]$npar

for (nm in setdiff(names(results), "M0_null")) {
  r <- results[[nm]]
  if (is.null(r) || !comp[comp$model == nm, "reliable"]) { cat(nm, ": not reliable, skipped\n"); next }
  lr_stat <- 2 * (r$ll - ll0); df_diff <- r$npar - npar0
  p_val <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
  cat(sprintf("%-28s LR=%.2f  df=%d  p=%.4f\n", nm, lr_stat, df_diff, p_val))
}

winner_name <- comp$model[comp$reliable][1]
winner_fit <- results[[winner_name]]$fit

## ================================================================
## SECTION 5: Winning Model
## ================================================================

state_labels <- c("State 1", "State 2", "State 3")
resp <- winner_fit@response
emean <- function(si) mean(sapply(resp[[si]], function(r) unname(r@parameters$coefficients[1])))
state_order <- order(sapply(seq_len(best_k), emean))
relabel <- setNames(state_labels[seq_len(best_k)], state_order)

vit <- posterior(winner_fit, type = "viterbi")
data_adj$state <- factor(relabel[as.character(vit$state)], levels = state_labels[seq_len(best_k)])

profile <- data_adj %>% group_by(state) %>%
  summarise(n_data_adj = n(), n_participants = n_distinct(pid),
            effort = round(mean(effort_z), 2),
            satisfaction = round(mean(sat_z), 2),
            trust = round(mean(trust_z), 2), .groups = "drop")
            
cat("\n--- SECTION 5: winning-model (", winner_name, ") state profile ---\n", sep = "")
print(as.data.frame(profile))

## ================================================================
## SECTION 5b: Uncertainty on the transition matrix (Table 4)
## ================================================================
## WHY: Table 4 reports 9 bare point estimates. A reader can't tell if
## .679 is really more persistent than .518, or whether the .000 cell
## means "impossible" or merely "not observed in 313 trials".
##
## METHOD A (Dirichlet): condition on the Viterbi path, count WITHIN-BLOCK
##   transitions, put a Dirichlet(1,..,1) prior on each row. Conjugate, so
##   each cell's posterior is a Beta -> closed-form, genuine HDIs, instant.
##   Fixes the .000 cell (posterior is NOT a point mass at zero).
##   Limitation: treats the decoded path as certain -> slightly too narrow.

## --- helpers ---
hdi_beta <- function(a, b, prob = 0.95, grid = 4000) {   # NARROWEST interval, not equal-tailed
  lows <- seq(0, 1 - prob, length.out = grid)
  lo <- qbeta(lows, a, b); hi <- qbeta(lows + prob, a, b)
  i <- which.min(hi - lo); c(lo[i], hi[i])
}
trans_of <- function(f, k) {   # extract + relabel by ascending emission mean
  A  <- matrix(unlist(lapply(f@transition, function(x) getpars(x))), nrow = k, byrow = TRUE)
  rr <- f@response
  em <- sapply(seq_len(k), function(si)
    mean(sapply(rr[[si]], function(r) unname(r@parameters$coefficients[1]))))
  A[order(em), order(em), drop = FALSE]
}

## --- METHOD A: Dirichlet posterior ---
seq_id <- paste(data_adj$pid, data_adj$block, sep = "_")
st     <- as.integer(data_adj$state)              # already relabeled in Section 5
cnt    <- matrix(0, best_k, best_k)
for (i in seq_len(nrow(data_adj) - 1))            # never count across a block boundary
  if (seq_id[i] == seq_id[i + 1]) cnt[st[i], st[i + 1]] <- cnt[st[i], st[i + 1]] + 1
dimnames(cnt) <- list(state_labels[seq_len(best_k)], state_labels[seq_len(best_k)])
cat("\n--- Observed within-block transition counts ---\n"); print(cnt)

PRIOR <- 1   # Dirichlet(1,..,1) = uniform on each row, weakly informative
trans_model <- trans_of(winner_fit, best_k)
trans_hdi <- expand.grid(from = state_labels[seq_len(best_k)],
                         to   = state_labels[seq_len(best_k)], stringsAsFactors = FALSE)
trans_hdi$model_prob <- as.vector(trans_model)    
trans_hdi$count <- trans_hdi$post_mean <- trans_hdi$hdi_lo <- trans_hdi$hdi_hi <- NA_real_
for (r in seq_len(nrow(trans_hdi))) {
  i <- match(trans_hdi$from[r], state_labels); j <- match(trans_hdi$to[r], state_labels)
  a_row <- cnt[i, ] + PRIOR
  a <- a_row[j]; b <- sum(a_row) - a_row[j]       # Beta marginal of the Dirichlet
  h <- hdi_beta(a, b, 0.95)
  trans_hdi$count[r]     <- cnt[i, j]
  trans_hdi$post_mean[r] <- a / (a + b)
  trans_hdi$hdi_lo[r] <- h[1]; trans_hdi$hdi_hi[r] <- h[2]
}
cat("\n--- METHOD A: Dirichlet 95% HDIs (use for Table 4) ---\n")
print(trans_hdi, row.names = FALSE, digits = 3)
write.csv(trans_hdi, file.path(data_dir, "transition_hdi_dirichlet.csv"), row.names = FALSE)

## ================================================================
## SECTION 6: Psychometric Correlations
## ================================================================
library(readr)

## design not built yet in data_adj at this point -- derive it here
data_adj$design <- ifelse(data_adj$task_design == "RANDOM", "pilot", "lived")

occupancy <- data_adj %>% count(pid, design, state, .drop = FALSE) %>%
  group_by(pid) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  dplyr::select(pid, design, state, pct) %>%
  tidyr::pivot_wider(names_from = state, values_from = pct, values_fill = 0)
cat("Occupancy (final states):\n"); print(as.data.frame(occupancy))

psych_raw <- read_csv(file.path(data_dir, "Qualtrics.csv"), skip = 2,
                      locale = locale(encoding = "latin1"),
                      show_col_types = FALSE) %>%
  rename(Participant.ID = `Participant ID`) %>%
  mutate(Participant.ID = as.character(Participant.ID)) %>%
  filter(!is.na(Participant.ID) & Participant.ID != "") %>%
  mutate(pid = gsub("\\s+", "", Participant.ID)) %>%
  as.data.frame()

score_row <- function(row) {
  lpfs_total <- sum(as.numeric(row[paste0("LPFS_", 1:12)]))
  pid_domain <- function(items) mean(as.numeric(row[paste0("PID_", items)]))
  pid_na    <- pid_domain(c(1,16, 6,21, 11,26))
  pid_det   <- pid_domain(c(4,19,  9,24, 14,29))
  pid_ant   <- pid_domain(c(2,17,  7,22, 12,27))
  pid_dis   <- pid_domain(c(3,18,  8,23, 13,28))
  pid_anank <- pid_domain(c(5,15, 10,20, 25,30))
  sapas_items <- as.numeric(row[paste0("SAPAS_", 1:8)]); sapas_items[3] <- 1 - sapas_items[3]
  sapas_total <- sum(sapas_items)
  ments <- as.numeric(row[paste0("MentS_", 1:12)]); names(ments) <- 1:12
  ments_self       <- mean(6 - ments[c("2","4","9","10")])
  ments_others     <- mean(ments[c("1","3","5","11")])
  ments_motivation <- mean(ments[c("6","7","8","12")])
  data.frame(pid = row[["pid"]], LPFS_Total = lpfs_total,
             PID_NegAffect = pid_na, PID_Detachment = pid_det, PID_Antagonism = pid_ant,
             PID_Disinhibition = pid_dis, PID_Anankastia = pid_anank, SAPAS_Total = sapas_total,
             MentS_Self = ments_self, MentS_Others = ments_others, MentS_Motivation = ments_motivation)
}
psych_scored <- bind_rows(lapply(seq_len(nrow(psych_raw)), function(i) score_row(psych_raw[i, ])))

merged <- inner_join(occupancy, psych_scored, by = "pid")
cat("\nMerged sample: n =", nrow(merged), "\n")
write.csv(merged, file.path(data_dir, "FINAL_CORRECTED_psychometrics_x_occupancy.csv"), row.names = FALSE)

psych_vars <- c("LPFS_Total","PID_NegAffect","PID_Detachment","PID_Antagonism",
                "PID_Disinhibition","PID_Anankastia","SAPAS_Total",
                "MentS_Self","MentS_Others","MentS_Motivation")

cor_table <- expand.grid(psych_var = psych_vars, state = state_labels[seq_len(best_k)], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(test = list(suppressWarnings(cor.test(merged[[psych_var]], merged[[state]], method = "spearman"))),
         rho = test$estimate, p = test$p.value) %>%
  ungroup() %>% dplyr::select(psych_var, state, rho, p) %>% arrange(psych_var, state)

## --- Holm-Bonferroni across all 30 tests -------------------------
## Holm is a STEP-DOWN procedure: sort p ascending, test the smallest
## against .05/30, next against .05/29, ... STOP at the first failure.
## Same FWER control as Bonferroni but uniformly more powerful, so there
## is no reason to report Bonferroni instead.
## NOTE: Holm's FIRST step uses the SAME threshold as Bonferroni, so if
## the smallest p fails there, Holm rejects nothing either -- identical
## conclusion, better procedure. That is expected here.
cor_table$p_holm  <- p.adjust(cor_table$p, method = "holm")
cor_table$sig_raw <- ifelse(cor_table$p      < 0.05, "*", "")   # uncorrected, reference only
cor_table$sig     <- ifelse(cor_table$p_holm < 0.05, "*", "")   # <- report THIS

## transparent walkthrough: shows exactly why each test does/doesn't survive
holm_walk <- cor_table %>% arrange(p) %>%
  mutate(rank = dplyr::row_number(),
         holm_threshold = 0.05 / (dplyr::n() - rank + 1),
         passes_step    = p <= holm_threshold)
ff <- which(!holm_walk$passes_step)[1]
holm_walk$survives_holm <- if (is.na(ff)) TRUE else
  c(rep(TRUE, ff - 1), rep(FALSE, nrow(holm_walk) - ff + 1))
cat("\n--- Holm-Bonferroni step-down (", nrow(holm_walk), " tests) ---\n", sep = "")
print(as.data.frame(holm_walk[, c("psych_var","state","rho","p","rank",
                                  "holm_threshold","p_holm","survives_holm")]), digits = 3)
cat("Surviving Holm at FWER .05:", sum(holm_walk$survives_holm), "of", nrow(holm_walk), "\n")

print(as.data.frame(cor_table), digits = 2)
write.csv(cor_table, file.path(data_dir, "FINAL_CORRECTED_psychometric_correlations.csv"), row.names = FALSE)

