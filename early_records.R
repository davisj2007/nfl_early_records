# ==============================================================================
# early_records.R
#
# Does an early-season record mean anything? Rebuild of the Sep 2022 thread
# (1999-2021, 16 games, empirical win table + binomial GAM playoff table),
# updated to 17 games and the 14-team field, plus two extensions:
#
#   1. Expected final wins by record      (empirical, GAM fallback in thin cells)
#   2. Playoff probability by record      (binomial GAM, reported at 14-team era)
#   3. Persistence by games played        (how much of a record's deviation
#                                           from .500 survives the rest of the
#                                           season; null = 0, "record is truth" = 1)
#   4. LOWFI incremental information      (within a record, does LOWFI separate
#                                           the real teams from the lucky ones?)

# ==============================================================================

# options, packages, and set up ---------
{
  options(scipen = 999) # no scientific notation
  library(tidyverse) # tidy data science
  library(tidymodels) # tidy machine learnings
  library(mgcv) # gam modeling
  library(nflreadr) # nfl data
  library(patchwork) # put plots together
  library(gt) # ggplot style tables

  # you can skip branding or add your own
  JDA_BRANDING  <- getOption("early_records.branding",      "../jda-branding/jda_branding.R")
  # this can be pulled from https://github.com/davisj2007/nfl_lowfi_data
  LOWFI_HISTORY <- getOption("early_records.lowfi_history", "../nfl_models26/outputs/history/lowfi_history.csv")
  LOWFI_SEASONS <- getOption("early_records.lowfi_seasons", 2004:2025)   # all retrodictions; all in-sample
  
  # where you want your outputs saved
  OUT_DIR       <- getOption("early_records.out_dir",       "outputs/early_records")

  # readers will need to comment out
  source(JDA_BRANDING)
  setup_jda_fonts()
  
  # create output directory
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  # paramaters
  SEASONS <- 1999:2025
  N_GAMES <- 17
  THIN_N  <- 10
)


# ------------------------------------------------------------------------------
# 1. Data: one row per team-season-state (after g games, g = 0..n_games)
# ------------------------------------------------------------------------------

sched <- load_schedules(SEASONS) %>%
  mutate(across(c(home_team, away_team), clean_team_abbrs))

playoff_teams <- sched %>%
  filter(game_type != "REG") %>%
  pivot_longer(c(home_team, away_team), values_to = "team") %>%
  distinct(season, team) %>%
  mutate(made_playoffs = 1L)

states <- sched %>%
  filter(game_type == "REG", !is.na(result)) %>%
  select(season, week, home_team, away_team, result) %>%
  pivot_longer(c(home_team, away_team), names_to = "side", values_to = "team") %>%
  mutate(margin = if_else(side == "home_team", result, -result),
         win    = case_when(margin > 0 ~ 1, margin < 0 ~ 0, TRUE ~ 0.5)) %>%
  group_by(season, team) %>%
  arrange(week, .by_group = TRUE) %>%
  mutate(g = row_number(), w = cumsum(win), l = cumsum(1 - win)) %>%
  select(season, team, week, g, w, l) %>%
  group_modify(~ add_row(.x, week = 0L, g = 0L, w = 0, l = 0, .before = 1)) %>%
  mutate(n_games = max(g), final_w = last(w)) %>%
  ungroup() %>%
  left_join(playoff_teams, by = c("season", "team")) %>%
  mutate(made_playoffs = coalesce(made_playoffs, 0L),
         rest_g   = n_games - g,
         rest_w   = final_w - w,
         rest_pct = if_else(rest_g > 0, rest_w / rest_g, NA_real_),
         cur_pct  = if_else(g > 0, w / g, NA_real_),
         W = floor(w), L = floor(l),
         era14 = as.integer(season >= 2020)) %>%
  rename(gweek = week)

# Sanity: 32 teams per season (31 before 2002), 16 or 17 games
states %>%
  distinct(season, team, n_games) %>%
  count(season, n_games) %>%
  print(n = Inf)


# ------------------------------------------------------------------------------
# 2. Models
# ------------------------------------------------------------------------------

# Rest-of-season win rate. Weighted quasibinomial: rest_pct is a proportion of
# rest_g games, and ties make rest_w non-integer.
fit_rest <- states %>%
  filter(rest_g > 0) %>%
  bam(rest_pct ~ te(w, l, k = c(8, 8)),
      family = quasibinomial(), weights = rest_g, data = ., discrete = TRUE)

# Playoffs. With a numeric 0/1 `by`, the era smooth is uncentered: it carries
# the era's level shift AND lets that shift vary by record (the two extra
# slots go to bubble teams, not 0-0 teams). Low k because it rests on six
# seasons.
fit_po <- states %>%
  bam(made_playoffs ~ te(w, l, k = c(8, 8)) + te(w, l, k = c(4, 4), by = era14),
      family = binomial(), data = ., discrete = TRUE)

summary(fit_rest)
summary(fit_po)

# Calibration-in-the-large, 14-team era: at every g, mean P(playoffs) across
# the 2020+ teams should sit near 14/32 = 0.4375.
states %>%
  filter(era14 == 1) %>%
  mutate(p = as.numeric(predict(fit_po, newdata = ., type = "response"))) %>%
  group_by(g) %>%
  summarise(n = n(), mean_pred = mean(p), actual = mean(made_playoffs), .groups = "drop") %>%
  print(n = Inf)


# ------------------------------------------------------------------------------
# 3. The grid: empirical + GAM, both tables
# ------------------------------------------------------------------------------

grid <- expand_grid(W = 0:N_GAMES, L = 0:N_GAMES) %>%
  filter(W + L <= N_GAMES) %>%
  mutate(w = W, l = L, era14 = 1L) %>%
  left_join(
    states %>%
      group_by(W, L) %>%
      summarise(n         = n(),
                rate_emp  = sum(rest_w[rest_g > 0]) / sum(rest_g[rest_g > 0]),
                po_emp_14 = mean(made_playoffs[era14 == 1]),
                n_14      = sum(era14 == 1),
                .groups = "drop"),
    by = c("W", "L")) %>%
  mutate(n         = coalesce(n, 0L),
         rate_gam  = as.numeric(predict(fit_rest, newdata = ., type = "response")),
         po_gam    = as.numeric(predict(fit_po,   newdata = ., type = "response")),
         rem       = N_GAMES - W - L,
         final_emp = if_else(rem == 0, as.numeric(W), W + rem * rate_emp),
         final_gam = W + rem * rate_gam,
         final_nul = W + rem * 0.5,
         thin      = n < THIN_N,
         final_show = if_else(thin | is.na(final_emp), final_gam, final_emp),
         vs_null   = final_show - final_nul)

write_csv(grid %>% select(-w, -l, -era14), file.path(OUT_DIR, "record_grid.csv"))

# Empirical vs GAM where both are well-supported: should be tight
grid %>%
  filter(!thin, rem > 0) %>%
  summarise(mae_final = mean(abs(final_emp - final_gam)),
            max_final = max(abs(final_emp - final_gam)),
            cells = n())


# ------------------------------------------------------------------------------
# 4a. Empirical: the spread of final wins behind each record
#     One panel per games played (1-8), x = losses, y = final wins on a
#     17-game basis (16-game seasons projected at their own rest-of-season
#     rate). Records with a tie are dropped so every box is a clean W-L.
#     Yellow diamond = GAM expectation; dashed line = coin flip the rest of
#     the way. varwidth: box width tracks how many team-seasons sit behind it.
# ------------------------------------------------------------------------------

box_dat <- states %>%
  filter(between(g, 1, 8), W + L == g) %>%
  mutate(final17 = w + rest_pct * (N_GAMES - g),
         panel   = factor(paste("After", g, if_else(g == 1, "game", "games")),
                          levels = paste("After", 1:8, c("game", rep("games", 7)))))

box_ref <- grid %>%
  filter(between(W + L, 1, 8)) %>%
  transmute(g = W + L, L, final_gam, final_nul,
            panel = factor(paste("After", g, if_else(g == 1, "game", "games")),
                           levels = levels(box_dat$panel)))

p_box <- box_dat %>%
  ggplot(aes(L, final17, group = L)) +
  geom_line(data = box_ref, aes(L, final_nul, group = 1),
            color = jda_silver, linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(varwidth = TRUE, fill = jda_navy, color = jda_cyan,
               outlier.color = jda_silver, outlier.size = 0.8, linewidth = 0.5) +
  geom_point(data = box_ref, aes(L, final_gam), inherit.aes = FALSE,
             shape = 23, size = 2.2, fill = jda_yellow, color = jda_black) +
  facet_wrap(~ panel, nrow = 2, scales = "free_x") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(limits = c(0, N_GAMES), breaks = seq(0, N_GAMES, 4)) +
  labs(title = "THE SAME RECORD ISN'T THE SAME OUTCOME",
       subtitle = "Final wins by current record. Yellow: GAM expectation. Dashed: .500 the rest of the way",
       x = "Losses", y = "Final wins (17-game basis)",
       caption =  "Regular seasons 1999-2025, records without ties. Box width scales with team-seasons. Jake Davis Analytics") +
  theme_jda()

save_jda_plot(file.path(OUT_DIR, "record_boxes.png"), p_box, 
              # width = 14, height = 8, dpi = 170, 
              bg = jda_black
              )

# Spread behind each early record, for the text
box_dat %>%
  group_by(g, record = paste0(W, "-", L)) %>%
  summarise(n = n(), median = median(final17),
            p10 = quantile(final17, 0.1), p90 = quantile(final17, 0.9), .groups = "drop") %>%
  filter(g <= 4) %>%
  print(n = Inf)

# The exceptions: every 0-3 or 0-4 start that still made the playoffs (footnote 5)
states %>%
  filter(between(g, 3, 4), W == 0, L == g, made_playoffs == 1) %>%
  group_by(season, team) %>%
  summarise(start = paste0("0-", max(g)), final_w, .groups = "drop") %>%
  distinct() %>%
  print(n = Inf)


# ------------------------------------------------------------------------------
# 4b. The two-panel grid (mirrors the 2022 layout), GAM throughout
# ------------------------------------------------------------------------------

p_wins <- grid %>%
  ggplot(aes(L, W)) +
  geom_tile(aes(fill = final_gam), color = jda_black, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", final_gam), color = final_gam > 10),
            size = 2.6, family = "Orbitron") +
  scale_fill_gradient(low = jda_navy, high = jda_cyan, limits = c(0, N_GAMES), guide = "none") +
  scale_color_manual(values = c(`TRUE` = jda_black, `FALSE` = jda_white), guide = "none") +
  scale_x_continuous(breaks = 0:N_GAMES, expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:N_GAMES, expand = c(0, 0)) +
  coord_equal() +
  labs(title = "Expected final wins", subtitle = "By current record, 17 games",
       x = "Losses", y = "Wins") +
  theme_jda() +
  theme(panel.grid = element_blank())

p_po <- grid %>%
  ggplot(aes(L, W, fill = po_gam)) +
  geom_tile(color = jda_black, linewidth = 0.4) +
  geom_text(aes(label = round(100 * po_gam), color = po_gam > 0.6),
            size = 2.6, family = "Orbitron") +
  scale_fill_gradient(low = jda_navy, high = jda_up, limits = c(0, 1), guide = "none") +
  scale_color_manual(values = c(`TRUE` = jda_black, `FALSE` = jda_white), guide = "none") +
  scale_x_continuous(breaks = 0:N_GAMES, expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:N_GAMES, expand = c(0, 0)) +
  coord_equal() +
  labs(title = "Playoff probability", subtitle = "Percent, by current record, 14-team field",
       x = "Losses", y = NULL) +
  theme_jda() +
  theme(panel.grid = element_blank())

p_grid <- (p_wins | p_po) +
  plot_annotation(
    title   = "EARLY GAMES ARE NOISY...BUT INFORMATIVE",
    caption = paste0("Regular seasons 1999-2025. Wins: GAM on rest-of-season win rate, projected to 17 games.\n",
                     "Playoffs: binomial GAM on wins and losses, with a 14-team-era adjustment. ",
                     "Jake Davis Analytics"),
    theme = theme_jda())

save_jda_plot(file.path(OUT_DIR, "record_grid.png"), p_grid, 
       width = 12, height = 8.5, 
       dpi = 170, bg = jda_black)

save_jda_plot(
# ------------------------------------------------------------------------------
# 5. Persistence: how much of a record's deviation from .500 survives?
#    (rest_pct - .5) = b * (cur_pct - .5); b = 0 coin flip, b = 1 record is truth
#
#    b is the reliability of the record as a read on true strength. Under a
#    beta-binomial, b(g) = g / (g + k): k is the number of .500 games you
#    "add" to a record to regress it. Fit k and draw that curve too.
#    g = 16 dropped: only 17-game seasons have a game left, n ~ 160.
# ------------------------------------------------------------------------------

persist <- states %>%
  filter(between(g, 1, 15), rest_g > 0) %>%
  nest(data = -g) %>%
  mutate(fit  = map(data, ~ lm(I(rest_pct - 0.5) ~ 0 + I(cur_pct - 0.5),
                               weights = rest_g, data = .x)),
         tidy = map(fit, ~ broom::tidy(.x, conf.int = TRUE)),
         n    = map_int(data, nrow)) %>%
  select(g, n, tidy) %>%
  unnest(tidy) %>%
  arrange(g)

k_fit <- persist %>%
  nls(estimate ~ g / (g + k), data = ., start = list(k = 10), weights = 1 / std.error^2)
k_hat <- coef(k_fit)[["k"]]

print(persist %>% select(g, n, estimate, conf.low, conf.high), n = Inf)
broom::tidy(k_fit, conf.int = TRUE) %>% print()
message(sprintf("Implied SD of true win%%: %.3f (%.1f wins over 17)",
                sqrt(0.25 / k_hat), 17 * sqrt(0.25 / k_hat)))

write_csv(persist %>% select(g, n, estimate, std.error, conf.low, conf.high),
          file.path(OUT_DIR, "persistence.csv"))

p_persist <- persist %>%
  ggplot(aes(g, estimate)) +
  geom_hline(yintercept = c(0, 1), color = jda_silver, linetype = "dashed") +
  annotate("text", x = 15, y = c(0.04, 0.96), hjust = 1, vjust = c(0, 1),
           label = c("coin flip the rest of the way", "record is the truth"),
           color = jda_silver, family = "Orbitron", size = 3) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = jda_cyan, alpha = 0.2) +
  geomtextpath::geom_textline(
    data = tibble(g = seq(1, 15, 0.1), estimate = g / (g + k_hat)),
    label = sprintf("add %.0f games of .500", k_hat),
    color = jda_yellow, linetype = "dashed", linewidth = 0.8, size = 3.2,
    hjust = 0.2, vjust = -0.4, family = "Orbitron") +
  geom_line(color = jda_cyan, linewidth = 1.2) +
  geom_point(color = jda_cyan, size = 2.5) +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(title = "HOW MUCH OF YOUR RECORD CARRIES FORWARD",
       subtitle = "Share of the distance from .500 that persists the rest of the way",
       x = "Games played", y = NULL,
       caption = "Weighted least squares by games played, 1999-2025. 95% intervals. Jake Davis Analytics") +
  theme_jda()

save_jda_plot(file.path(OUT_DIR, "persistence.png"), p_persist, 
              #width = 10, height = 6, dpi = 170, 
              bg = jda_black)


# ------------------------------------------------------------------------------
# 6. LOWFI: does the model know which 3-1 teams are real?
#    History `week` is the NFL week (games through that week reflected); the
#    history has no row for a team on its bye, so 15 rows per team in wk 1-16.
#    Expect 480 rows per season once abbreviations are normalized.
# ------------------------------------------------------------------------------

lowfi_raw <- read_csv(LOWFI_HISTORY, show_col_types = FALSE) %>%
  filter(season %in% LOWFI_SEASONS, source == "retrodiction") %>%
  mutate(team = clean_team_abbrs(team))

preseason <- lowfi_raw %>%
  filter(week == 0) %>%
  select(season, team, net0 = net_rating)
has_pre <- nrow(preseason) > 0
if (!has_pre) message("No week-0 rows in the history: skipping the preseason-only model")

lowfi_dat <- lowfi_raw %>%
  filter(between(week, 1, 16)) %>%
  select(season, week, team, net_rating, p_hat) %>%
  left_join(preseason, by = c("season", "team")) %>%
  inner_join(states, by = join_by(season, team, closest(week >= gweek))) %>%
  filter(rest_g > 0, !is.na(net_rating), !is.na(p_hat)) %>%
  mutate(logit_p0 = qlogis(p_hat))

# Join check: every team-season-week should land exactly once
lowfi_dat %>% count(season, week, team) %>% filter(n > 1) %>% nrow() %>%
  { if (. > 0) stop(., " duplicated team-weeks in the LOWFI join", call. = FALSE) }
lowfi_dat %>% count(season) %>% print(n = Inf)

specs <- tibble(
  model   = c("Record", "LOWFI/0", "LOWFI/1", "LOWFI/1 preseason only",
              "Record + LOWFI/0", "Record + LOWFI/1"),
  formula = list(rest_pct ~ cur_pct,
                 rest_pct ~ logit_p0,
                 rest_pct ~ net_rating,
                 rest_pct ~ net0,
                 rest_pct ~ cur_pct + logit_p0,
                 rest_pct ~ cur_pct + net_rating)) %>%
  filter(has_pre | model != "LOWFI/1 preseason only")

# Leave-one-season-out within each week. Skill = 1 - MSE / MSE(coin flip).
lowfi_skill <- lowfi_dat %>%
  nest(data = -week) %>%
  mutate(folds = map(data, ~ group_vfold_cv(.x, group = season))) %>%
  select(week, folds) %>%
  unnest(folds) %>%
  expand_grid(specs) %>%
  mutate(pred = map2(splits, formula, \(s, f)
    assessment(s) %>%
      transmute(season, rest_pct, rest_g,
                pred = predict(glm(f, family = quasibinomial(), weights = rest_g,
                                   data = analysis(s)),
                               newdata = assessment(s), type = "response")))) %>%
  select(week, model, pred) %>%
  unnest(pred) %>%
  filter(!is.na(pred)) %>%
  group_by(week, model) %>%
  summarise(mse  = weighted.mean((rest_pct - pred)^2, rest_g),
            mse0 = weighted.mean((rest_pct - 0.5)^2,  rest_g),
            .groups = "drop") %>%
  mutate(skill = 1 - mse / mse0)

lowfi_skill %>%
  select(week, model, skill) %>%
  pivot_wider(names_from = model, values_from = skill) %>%
  print(n = Inf, width = Inf)
write_csv(lowfi_skill, file.path(OUT_DIR, "lowfi_skill.csv"))

p_skill <- lowfi_skill %>%
  filter(model %in% c("Record", "Record + LOWFI/0", "Record + LOWFI/1", "LOWFI/1 preseason only")) %>%
  ggplot(aes(week, skill, color = model, linetype = model)) +
  geom_hline(yintercept = 0, color = jda_silver, linetype = "dashed") +
  # Label positions per line so the four labels don't stack on each other
  geomtextpath::geom_textline(
    aes(label = model,
        hjust = case_match(model, "Record" ~ 0.12, "Record + LOWFI/0" ~ 0.08,
                           "Record + LOWFI/1" ~ 0.45, "LOWFI/1 preseason only" ~ 0.72)),
    linewidth = 1.2, size = 3.6, vjust = -0.3, family = "Orbitron") +
  scale_color_manual(values = c("Record" = jda_silver,
                                "Record + LOWFI/0" = jda_pink,
                                "Record + LOWFI/1" = jda_cyan,
                                "LOWFI/1 preseason only" = jda_yellow),
                     guide = "none") +
  scale_linetype_manual(values = c("Record" = "solid",
                                   "Record + LOWFI/0" = "solid",
                                   "Record + LOWFI/1" = "solid",
                                   "LOWFI/1 preseason only" = "dashed"),
                        guide = "none") +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "THE RECORD KNOWS SOMETHING, BUT LOWFI KNOWS MORE",
       subtitle = "Skill predicting rest-of-season win rate, vs. a coin flip",
       x = "Week", y = NULL,
       caption = paste0("Leave-one-season-out, ", min(LOWFI_SEASONS), "-", max(LOWFI_SEASONS),
                        " LOWFI retrodictions (in-sample).\nJake Davis Analytics")) +
  theme_jda()

save_jda_plot(file.path(OUT_DIR, "lowfi_skill.png"), p_skill, 
             # width = 10, height = 6, dpi = 170, 
              bg = jda_black)


# ------------------------------------------------------------------------------
# 7. Same record, different teams: week 4 split by LOWFI/1 within record
# ------------------------------------------------------------------------------

split_wk <- 4

tbl_split <- lowfi_dat %>%
  filter(week == split_wk) %>%
  mutate(record = paste0(W, "-", L)) %>%
  group_by(record) %>%
  filter(n() >= 20) %>%
  mutate(half = if_else(net_rating >= median(net_rating), "top", "bottom"),
         final17 = w + rest_pct * (N_GAMES - g)) %>%
  group_by(record, W, half) %>%
  summarise(n = n(), final = mean(final17), po = mean(made_playoffs), .groups = "drop") %>%
  pivot_wider(names_from = half, values_from = c(n, final, po)) %>%
  mutate(gap = final_top - final_bottom) %>%
  arrange(desc(W)) %>%
  select(record, n_top, final_top, po_top, n_bottom, final_bottom, po_bottom, gap)

print(tbl_split)

tbl_split %>%
  gt() %>%
  cols_label(record = "Record",
             n_top = "n", final_top = "Final wins", po_top = "Playoffs",
             n_bottom = "n", final_bottom = "Final wins", po_bottom = "Playoffs",
             gap = "Gap") %>%
  tab_spanner("Top half by LOWFI/1",    c(n_top, final_top, po_top)) %>%
  tab_spanner("Bottom half by LOWFI/1", c(n_bottom, final_bottom, po_bottom)) %>%
  fmt_number(c(final_top, final_bottom, gap), decimals = 1) %>%
  fmt_percent(c(po_top, po_bottom), decimals = 0) %>%
  tab_header(title = "WHAT DOES LOWFI/1 KNOW?",
             subtitle = paste0("After week ", split_wk, ", split at each record's median LOWFI/1 rating")) %>%
  tab_source_note(paste0("Final wins projected to 17 games. ", min(LOWFI_SEASONS), "-", max(LOWFI_SEASONS),
                         " retrodictions; playoffs across both eras. Jake Davis Analytics")) %>%
  gt_theme_jda() %>%
  save_jda_table_substack(file.path(OUT_DIR, "lowfi_split_wk04.png"))

