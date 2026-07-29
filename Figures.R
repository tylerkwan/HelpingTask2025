## ================================================================
## HMM FIGURES
## ================================================================
data_adj <- data_adj %>%
  mutate(pid_clean = ifelse(grepl("^S", pid), sprintf("P%02d", match(pid, sort(unique(pid[grepl("^S", pid)])))), pid))

library(tidyr)
install.packages("ggridges")
library(ggridges)  


pal <- setNames(c("#ef6f6c", "#7cb342", "#3b82f6")[seq_len(best_k)], state_labels[seq_len(best_k)])
pctlab <- function(x) paste0(round(100 * x), "%")

## design (pilot/lived) isn't in data_adj yet -- derive it here
data_adj$design <- ifelse(data_adj$task_design == "RANDOM", "pilot", "lived")

## ================================================================
## FIGURE  1. STATE PROFILES: z-score heatmap (what defines each state)
## ================================================================
hm <- data_adj %>% group_by(state) %>%
  summarise(Effort = mean(effort_z), Satisfaction = mean(sat_z), Trust = mean(trust_z), .groups = "drop") %>%
  pivot_longer(-state, names_to = "channel", values_to = "zv") %>%
  mutate(channel = factor(channel, levels = c("Effort","Satisfaction","Trust")))

fig_profiles <- ggplot(hm, aes(channel, state)) +
  geom_tile(aes(fill = zv), color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%+.2f", zv)), fontface = "bold", size = 6, color = "grey10") +
  scale_fill_gradient2(low = "#ef4444", mid = "#f8fafc", high = "#3b82f6", midpoint = 0,
                        limits = c(-1.5, 1.5), breaks = c(-1, 0, 1),
                        labels = c("-1 SD", "0", "+1 SD")) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(title = "",
       x = "", y = "", fill = "z-score") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.margin = margin(12, 18, 10, 12))

fig_profiles

ggsave(file.path(data_dir, "FIG_01_state_profiles.png"), fig_profiles, width = 9, height = 4.3, dpi = 150)


## ================================================================
## 2. STATE OCCUPANCY -- % time each person spent in each state
## ================================================================
occ <- data_adj %>% count(pid_clean, design, state, .drop = FALSE) %>%
  group_by(pid_clean) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
order_key <- occ %>% filter(state == "State 3") %>% arrange(pct) %>% pull(pid_clean)
order_key <- c(setdiff(unique(occ$pid_clean), order_key), order_key)
occ <- occ %>%
  mutate(pid_label = factor(pid_clean, levels = order_key))

fig_occupancy <- ggplot(occ, aes(pct, pid_label, fill = state)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(pct >= 8, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE), size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = pal) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.04)), breaks = seq(0,100,25)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "",
       x = "% of trials", y = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        panel.grid.major.y = element_blank(), legend.position = "top",
        plot.margin = margin(12, 20, 10, 12))

fig_occupancy
ggsave(file.path(data_dir, "FIG_02_state_occupancy.png"), fig_occupancy, width = 10, height = 6.5, dpi = 150)


## ================================================================
## 3. SELF vs OTHER EVALUATION PLANE
## ================================================================
fig_plane <- ggplot(data_adj, aes(trust, satisfaction, color = state)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_jitter(width = 0.12, height = 0.12, size = 2.6, alpha = 0.75) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = 1:6, limits = c(0.7, 6.3)) +
  scale_y_continuous(breaks = 1:6, limits = c(0.7, 6.3)) +
  labs(title = "",
       x = "Trust in partner (1-6)", y = "Satisfaction with self (1-6)", color = "State") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.margin = margin(12, 18, 10, 12))

fig_plane
ggsave(file.path(data_dir, "FIG_03_self_other_plane.png"), fig_plane, width = 8.5, height = 6.2, dpi = 150)


## ================================================================
## 4. STATES OVER TIME (within block)
## ================================================================
block_pos <- data_adj %>% filter(!is.na(block_trial)) %>%
  count(block_trial, state, .drop = FALSE) %>%
  group_by(block_trial) %>% mutate(prop = n / sum(n)) %>% ungroup()

## optional diagnostic: where does the tail thin out?
block_pos %>% group_by(block_trial) %>% summarise(n_trials = sum(n)) %>% print(n = Inf)

fig_overtime <- ggplot(block_pos,
                       aes(x = block_trial, y = state, height = prop, fill = state)) +
  geom_ridgeline(scale = 2.2, alpha = 0.85, colour = "grey25", linewidth = 0.4) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_discrete(limits = rev, expand = expansion(add = c(0.2, 1.5))) +
  scale_x_continuous(breaks = sort(unique(block_pos$block_trial))) +
  labs(x = "Trial within block", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.margin = margin(12, 18, 10, 12))

fig_overtime
ggsave(file.path(data_dir, "FIG_states_over_time.png"), fig_overtime,
       width = 8.5, height = 4.5, dpi = 300)

## ================================================================
## 5. TRAJECTORIES -- individual state sequence, Lived-Exp participants
## ================================================================
le <- data_adj %>% filter(design == "lived") %>% arrange(pid_clean, block, block_trial)

blocks <- le %>% group_by(pid, block) %>%
  summarise(xmin = min(trial) - 0.5, xmax = max(trial) + 0.5, policy = first(policy), .groups = "drop")

fig_traj <- ggplot(le, aes(trial, state)) +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, fill = policy),
            ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.15) +
  geom_line(aes(group = pid), color = "grey50") +
  geom_point(aes(color = state), size = 3) +
  facet_wrap(~pid, ncol = 1) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = c("CARING" = "#2ca02c", "UNCARING" = "#ff7f0e"), 
                    labels = c("CARING" = "Caring", "UNCARING" = "Uncaring"),
                    name = "Helper behaviour") +
  labs(title = "",
       x = "trial", y = "State") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        strip.text = element_text(face = "bold"), plot.margin = margin(12, 22, 10, 12))
fig_traj

ggsave(file.path(data_dir, "FIG_05_trajectories.png"), fig_traj, width = 9.5, height = 6, dpi = 150)


## ================================================================
## 6. STATE COMPOSITION -- Pooled Pilots vs LE
## ================================================================
comp_state <- data_adj %>% mutate(who = ifelse(design == "pilot", "Pilot participants", pid)) %>%
  count(who, state, .drop = FALSE) %>% group_by(who) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(who = factor(who, levels = c("Pilot participants", "LE1", "LE2")))



fig_composition <- ggplot(comp_state, aes(who, prop, fill = state)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(prop > 0.05, pctlab(prop), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE), fontface = "bold", size = 3.6, color = "white") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = pctlab, expand = expansion(mult = c(0, 0.02))) +
  labs(title = "", subtitle = "",
       x = NULL, y = "Proportion of trials", fill = "State") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.margin = margin(12, 18, 10, 12))

fig_composition

ggsave(file.path(data_dir, "FIG_06_state_composition.png"), fig_composition, width = 8, height = 5.3, dpi = 150)


## ================================================================
## 7. COVARIATE RESULTS -- BIC comparison across all tested specs
## ================================================================
comp_plot <- comp %>% filter(!is.na(BIC))
comp_plot$model <- factor(comp_plot$model, levels = rev(comp_plot$model))


clean_names <- c("M0_null" = "Null (no covariate)",
                 "M1_effort_z" = "Helper's effort",
                 "M2_returns_z" = "Participant returns",
                 "M3_returns_rank_z" = "Participant returns (rank-ordered)",
                 "M4_ineq_z" = "Payoff disparity",
                 "M5_effort_and_ineq_z" = "Helper's effort x Payoff disparity")

fig_covariates <- ggplot(comp_plot, aes(x = BIC, y = model, fill = reliable)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = results[["M0_null"]]$bic, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = BIC), hjust = -0.1, size = 3.2) +
  scale_fill_manual(values = c("TRUE" = "#3b82f6", "FALSE" = "grey80"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(labels = clean_names) +
  labs(title = "",
       x = "BIC", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.margin = margin(12, 24, 10, 12))

fig_covariates

ggsave(file.path(data_dir, "FIG_08_covariate_results.png"), fig_covariates, width = 9, height = 5, dpi = 150)

# NOTE: for figure caption Lower BIC = better. Dashed line = null model (transition ~ 1).


## ================================================================
## 8. TRANSITION MATRIX 
## ================================================================
trans_raw <- matrix(unlist(lapply(winner_fit@transition, function(x) getpars(x))), nrow = best_k, byrow = TRUE)
resp <- winner_fit@response
emean <- function(si) mean(sapply(resp[[si]], function(r) unname(r@parameters$coefficients[1])))
state_order <- order(sapply(seq_len(best_k), emean))
trans <- trans_raw[state_order, state_order]
dimnames(trans) <- list(paste0("from: ", state_labels[seq_len(best_k)]), paste0("to: ", state_labels[seq_len(best_k)]))

trans_df <- as.data.frame(as.table(trans))
names(trans_df) <- c("from", "to", "prob")
trans_df$from <- factor(gsub("^from: ", "", trans_df$from), levels = state_labels[seq_len(best_k)])
trans_df$to   <- factor(gsub("^to: ",   "", trans_df$to),   levels = state_labels[seq_len(best_k)])
diag_df <- subset(trans_df, as.character(from) == as.character(to))

fig_transmat <- ggplot(trans_df, aes(to, from, fill = prob)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_tile(data = diag_df, fill = NA, color = "black", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * prob)), fontface = "bold", size = 4.8) +
  scale_fill_gradient(low = "#f8fafc", high = "#3b82f6", guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = 0.6)) +
  scale_y_discrete(limits = rev, expand = expansion(add = 0.6)) +
  labs(title = "",
       x = "Next trial", y = "Current trial") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.margin = margin(12, 18, 10, 12),
        axis.title.x.top = element_text(margin = margin(b = 12)),
        axis.title.y = element_text(margin = margin(r = 12)))

fig_transmat
ggsave(file.path(data_dir, "FIG_09_transition_matrix.png"), fig_transmat, width = 7, height = 6, dpi = 150)


## ================================================================
## 9. TRANSITION PREDICTIONS -- since the null won, there is no
##    covariate to plot a prediction curve against for the WINNING
##    model. Built instead for the best-performing NON-null candidate,
##    clearly labeled as exploratory context, not the reported result.
## ================================================================
runner_up_name <- comp$model[comp$reliable & comp$model != "M0_null"][1]
runner_up_fit  <- results[[runner_up_name]]$fit
runner_up_terms <- all.vars(specs[[runner_up_name]])

if (length(runner_up_terms) == 1) {
  grid <- seq(-2, 2, by = 0.2)
  pred_df <- do.call(rbind, lapply(seq_len(best_k), function(from_state) {
    do.call(rbind, lapply(grid, function(cv) {
      nd <- matrix(c(1, cv), nrow = 1)
      pr <- predict(runner_up_fit@transition[[from_state]], nd)
      data.frame(from = state_labels[from_state], covariate = cv,
                  to = state_labels[seq_len(best_k)], prob = as.numeric(pr))
    }))
  }))
  pred_df$from <- factor(pred_df$from, levels = state_labels[seq_len(best_k)])
  pred_df$to   <- factor(pred_df$to,   levels = state_labels[seq_len(best_k)])
  fig_predictions <- ggplot(pred_df, aes(covariate, prob, color = to)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~from, labeller = labeller(from = function(x) paste("FROM:", x))) +
    scale_color_manual(values = pal) +
    labs(title = paste(""),
         subtitle = paste0("BIC = ", round(comp$BIC[comp$model==runner_up_name],1),
                            " vs null's ", round(results[["M0_null"]]$bic,1)),
         x = paste0("Rank-ordered participant returns (z-scored)"), y = "P(next state)", color = "Next state") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
          strip.text = element_text(face = "bold"), plot.margin = margin(12, 18, 10, 12))
}

fig_predictions
ggsave(file.path(data_dir, "FIG_11_transition_predictions.png"), fig_predictions, width = 10, height = 4.2, dpi = 150)

## ================================================================
## 10. EMISSION PROBABILITIES -- the fitted Gaussian density each
##     state actually uses to generate effort/satisfaction/trust
##     (not the empirical histogram -- the model's own distributions)
## ================================================================

emissions <- bind_rows(lapply(seq_len(best_k), function(si) {
  bind_rows(lapply(seq_along(resp[[si]]), function(ci) {
    p <- resp[[si]][[ci]]@parameters
    data.frame(state = state_labels[match(si, state_order)], channel = c("Effort","Satisfaction","Trust")[ci],
               mean = unname(p$coefficients[1]), sd = unname(p$sd))
  }))
}))
emissions$state <- factor(emissions$state, levels = state_labels)
emissions$channel <- factor(emissions$channel, levels = c("Effort","Satisfaction","Trust"))

dens <- emissions %>%
  rowwise() %>%
  reframe(x = seq(mean - 4*sd, mean + 4*sd, length.out = 200),
          density = dnorm(x, mean, sd),
          state = state, channel = channel)

fig_emissions <- ggplot(dens, aes(x, density, color = state, fill = state)) +
  geom_line(linewidth = 1) +
  geom_area(alpha = 0.15, position = "identity") +
  facet_wrap(~channel, scales = "fixed", ncol = 3) +
  scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
  labs(title = "",
       x = "z-score", y = "Density", color = "State", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        strip.text = element_text(face = "bold"), legend.position = "top",
        plot.margin = margin(12, 18, 10, 12))

fig_emissions 

# CAPTION: subtitle = "The fitted Gaussian each state uses to generate each channel (z-scored units)",

ggsave(file.path(data_dir, "FIG_07_emission_probabilities.png"), fig_emissions, width = 10, height = 4.5, dpi = 150)

## ================================================================
## 11. Psychometric Correlations
## ================================================================


fig_corr <- ggplot(cor_table, aes(state, psych_var, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sig)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2c6e49", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman\nrho") +
  scale_x_discrete(expand = expansion(add = 0.6)) + scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(title = "",
       subtitle = paste0("* p<.05 (uncorrected)"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        panel.grid = element_blank(), plot.margin = margin(12, 18, 10, 12))

fig_corr

ggsave(file.path(data_dir, "FIG_10_psychometric_correlations.png"), fig_corr, width = 8.5, height = 6.5, dpi = 150)
