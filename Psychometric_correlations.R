suppressMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
data_dir <- "."

## ================================================================
## Load the FINAL (corrected, realigned, M0_null-winning) decoded states
## ================================================================
d <- read.csv(file.path(data_dir, "FINAL_decoded_for_psych.csv"), stringsAsFactors = FALSE)
d$design <- ifelse(d$task_design == "RANDOM", "pilot", "lived")
state_labels <- c("Wary", "Ambivalent", "Cooperative")
d$state <- factor(d$state, levels = state_labels)

occupancy <- d %>% count(pid, design, state, .drop = FALSE) %>%
  group_by(pid) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  select(pid, design, state, pct) %>%
  pivot_wider(names_from = state, values_from = pct, values_fill = 0)
cat("Occupancy (final, corrected states):\n"); print(as.data.frame(occupancy))

## ================================================================
## Load and score psychometrics (validated keys -- reproduce appendix
## descriptives exactly; see prior validation in this project)
## ================================================================
psych_raw <- readr::read_csv(file.path(data_dir, "Qualtrics.csv"), skip = 2,
                               locale = readr::locale(encoding = "latin1"),
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
cat("\nScored psychometrics for", nrow(psych_scored), "participants\n")

## ================================================================
## Merge and correlate (Spearman; exploratory, n=12)
## ================================================================
merged <- inner_join(occupancy, psych_scored, by = "pid")
cat("\nMerged sample: n =", nrow(merged), "\n")
write.csv(merged, file.path(data_dir, "FINAL_CORRECTED_psychometrics_x_occupancy.csv"), row.names = FALSE)

psych_vars <- c("LPFS_Total","PID_NegAffect","PID_Detachment","PID_Antagonism",
                  "PID_Disinhibition","PID_Anankastia","SAPAS_Total",
                  "MentS_Self","MentS_Others","MentS_Motivation")
cor_table <- expand.grid(psych_var = psych_vars, state = state_labels, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(test = list(suppressWarnings(cor.test(merged[[psych_var]], merged[[state]], method = "spearman"))),
         rho = test$estimate, p = test$p.value) %>%
  ungroup() %>% select(psych_var, state, rho, p) %>% arrange(psych_var, state)
cat("\n--- Spearman correlations: personality x state occupancy (FINAL, corrected states) ---\n")
print(as.data.frame(cor_table), digits = 2)
write.csv(cor_table, file.path(data_dir, "FINAL_CORRECTED_psychometric_correlations.csv"), row.names = FALSE)

cor_table$sig <- ifelse(cor_table$p < 0.05, "*", "")
fig <- ggplot(cor_table, aes(state, psych_var, fill = rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sig)), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2c6e49", midpoint = 0,
                        limits = c(-1, 1), name = "Spearman\nrho") +
  scale_x_discrete(expand = expansion(add = 0.6)) + scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(title = "Personality indices x time spent in each state (final, corrected states)",
       subtitle = paste0("Spearman correlations, n=", nrow(merged), "  |  * p<.05 (uncorrected, exploratory)"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        panel.grid = element_blank(), plot.margin = margin(12, 18, 10, 12))
ggsave(file.path(data_dir, "FIG_11_psychometric_correlations.png"), fig, width = 8.5, height = 6.5, dpi = 150)
cat("\nSaved FIG_11_psychometric_correlations.png\n")
