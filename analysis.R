library(tidyverse)
library(stringr)

## Read in India data
# Proportion samples
dhs4 <- readRDS("dirlecht_dhs_earliest.Rds")
dhs7 <- readRDS("dirlecht_dhs_latest.Rds")
# Realtive risks, order and actual
rfs <- readRDS("dirlecht_risk_factor_order.Rds")
rrs <-  read_csv('Relative_risks_factors_pneumonia.csv')
# Population data
pop <- read_csv("indiau5pop.csv") %>% 
  filter(year %in% c(2000, 2015))
# write_csv(pop, "Outputs/pop_India.csv")
# Names of states - lookup file from DHS codes to Abberviations and state names used
# in population file
state_lkp <- read_csv("lookup_states_india.csv") 
# Hib correction
hib_corr <- readxl::read_excel("Hib_adjustment.xlsx")

# Count samples
n.samples <- nrow(dhs4$`0`)


# Obtain beta parameters for Hib correction
source("scripts/beta_param_quantiles.R")
hib_corr_dstnct <- hib_corr %>% 
  select(clinical_hib_lb, clinical_hib_ub, severe_hib_lb, severe_hib_ub) %>% 
  distinct() %>% 
  gather("key", "value") %>% 
  mutate(severity = if_else(str_detect(key, "severe"), "severe", "clinical")) %>% 
  mutate(key = c(1,2,1,2)) %>% 
  spread(severity, value) %>% 
  select(-key)

beta_params_hib <- map(hib_corr_dstnct, ~ beta.parms.from.quantiles(.x)[c("a", "b")])

## For community and severe, sample from beta
indian_hib_proportion_community <- rbeta(n.samples, beta_params_hib$clinical$a, beta_params_hib$clinical$b)
indian_hib_proportion_severe    <- rbeta(n.samples, beta_params_hib$severe$a, beta_params_hib$severe$b)

# Distinct state data for hib correction
hib_corr <- hib_corr %>% 
  distinct(year, state, state_code, hib3_vec_ind_adj)

# Add state abbreviations to DHS data
## Relabel states with abbreviations 
DHSLKP <- function(dhs) {
  dhs_lkp <- state_lkp$abbrev[state_lkp$dhs == dhs]
  names(dhs_lkp) <- as.character(state_lkp$code[state_lkp$dhs==dhs])
  dhs_lkp
}

dhs_earliest_lkp <- DHSLKP(4)
dhs_latest_lkp <- DHSLKP(7) # Note there are different codes in 5, 6 and 7
names(dhs4) <- dhs_earliest_lkp[names(dhs4)]
names(dhs7) <- dhs_latest_lkp[names(dhs7)]

# Check all present
setdiff(state_lkp$state_pop_name, pop$state)
setdiff(pop$state, state_lkp$state_pop_name)

setdiff(names(dhs4), state_lkp$abbrev)
setdiff(names(dhs7), state_lkp$abbrev)

# Impute for selected states
#Do imputation for states where no survey. As per discussion with Brian Wahl and Harish, these are as follows.

msngs <- "msgn,equiv, abbrev_equiv, abbrev_msng
Andaman & Nicobar Islands, West Bengal, WB, AN
Chandigarh, Punjab, PJ, CH
Dadra & Nagar Haveli, Gujarat, GJ, DN
Daman & Diu,Gujarat, GJ, DD
Lakshadweep, Kerala, KE, LK
Puducherry, Tamil Nadu, TN, PD
Chhattisgarh, Madhya Pradesh, MP, CG
Jharkhand, Bihar, BH, JH
Uttarakhand, Uttar Pradesh, UP, UC"

msngs <- readr::read_csv(msngs) 

# No missing states for DHS 7
intersect(msngs$abbrev_msng, names(dhs4))
intersect(msngs$abbrev_msng, names(dhs7))

ImputeStateFunction <- function(res_states){
  print ("Impute using the following states")
  print(intersect(names(res_states), msngs$abbrev_equiv))
  res_states_for_imput <- res_states[msngs$abbrev_equiv]
  names(res_states_for_imput) <- msngs$abbrev_msng
  c(res_states, res_states_for_imput)
}

dhs4 <- ImputeStateFunction(dhs4)

## Character vector of risk factors
rf_char <- names(rfs)

## Seed the random number generator and the number of samples ----
set.seed(12345)
n.samples <- nrow(dhs4$IA)

## relative risks ----
## in the last analysis.

## Order the risk factors so that they correspond to the
## correct sampled prevalences (i.e. the rrs match the order
## of the risk factors in the data.to.analyse
rrs <- as.data.frame(rrs)
row.names(rrs) <- rrs$`Risk Factor`

# Ensure same order as proportion data
rrs <- rrs[c(
  "lbw" = "Low birth weight (=<2500 g)",
  "inc.imm" = "Incomplete immunization in 1st year of life",
  "non.ebf" = "Non-breastfed exclus. (4 mths)",
  "maln.long" = "Malnutrition (wght-for-age z<-2)",
  "overcrowd" = "Crowding (7 or more persons)",
  "indoor.air.long" = "Use solid fuels (yes)",
  "hiv" = "HIV" ),]
ses1 <- with(rrs, (log(HCI) - log(OR)) / 1.96)
ses2 <- with(rrs, (log(OR) - log(LCI)) / 1.96)
rrs$ses <- pmax(ses1, ses2)

number_rfs <- nrow(rrs)
log.rr.samples <-
  with(rrs, rnorm(number_rfs * n.samples, log(OR), ses))
log.rr.samples <- matrix(log.rr.samples, nrow = n.samples, byrow = TRUE)

## Prepare to loop through creating risk factor combinations
rfs_mtrx <- as.matrix(rfs)
num_rfs_combined <- rowSums(rfs_mtrx)
modifier_combined <- num_rfs_combined * log(0.8)

log.rr.samples.list <- vector(length = 128, mode = "list")

# Loop through creating risk factor combinations
for (choose_row in 1:128){
  multimat <- matrix(rfs_mtrx[choose_row,], nrow = n.samples, ncol = 7, byrow = TRUE)
  log.rr.samples.list[[choose_row]] <- rowSums(log.rr.samples * multimat) + 
    modifier_combined[choose_row]
}
log.rr.samples.matrix <- do.call(cbind, log.rr.samples.list)
rm(ses1, ses2, log.rr.samples, log.rr.samples.list, multimat, rfs_mtrx)


## Normalise prevalence for non-HIV and HIV samples
Normalise <-  function(x, cols){
  # as some sites have no HIV need to set these to zero
  x <- x[,cols]
  x_sum <- rowSums(x)
  a <- x/x_sum
  # print(dim(a))
  # Set NaN to very close to zero
  a[is.nan(a)] <- 0.0000000000001
  a
}
dhs4        <- map(dhs4, Normalise, cols = 1:128)
dhs4_hiv    <- map(dhs4, Normalise, cols = 65:128)
dhs4_no_hiv <- map(dhs4, Normalise, cols = 1:64)
dhs7        <- map(dhs7, Normalise, cols = 1:128)
dhs7_hiv    <- map(dhs7, Normalise, cols = 65:128)
dhs7_no_hiv <- map(dhs7, Normalise, cols = 1:64)

dhs4_by_hiv <- list(all = dhs4, hiv = dhs4_hiv, no_hiv = dhs4_no_hiv)
dhs7_by_hiv <- list(all = dhs7, hiv = dhs7_hiv, no_hiv = dhs7_no_hiv)

## Divide RRs for hiv, no hiv and all hiv
log.rr.samples.matrix_by_hiv <- list(all = log.rr.samples.matrix,
                                     hiv = log.rr.samples.matrix[,65:128],
                                     no_hiv = log.rr.samples.matrix[, 1:64])

### Estimate sum(p_j*rr_j)for both sets of surveys for all, hiv and no hiv
## Prevalence x rr for each rr combination for DHS4
pj_dhs4_by_hiv <- map2(dhs4_by_hiv, log.rr.samples.matrix_by_hiv,
                       function(dhs4, log.rr.samples.matrix) {
                         map(dhs4, function(x) exp(log(x) +log.rr.samples.matrix)) # each state
                       }) # each all, hiv, no_hiv group

## Sum of prevalence x rr for each state
sum_pj_dhs4_by_hiv <- map(pj_dhs4_by_hiv, ~ map(.x, rowSums))

## Prevalence x rr for each rr combination for DHS7
pj_dhs7_by_hiv <- map2(dhs7_by_hiv, log.rr.samples.matrix_by_hiv,
                       function(dhs7, log.rr.samples.matrix) {
                         map(dhs7, function(x) exp(log(x) + log.rr.samples.matrix)) # each state
                       }) # each all, hiv, no_hiv group
## Sum of prevalence x rr for each state
sum_pj_dhs7_by_hiv <- map(pj_dhs7_by_hiv, ~ map(.x, rowSums))

### Estimate Incidence for each state at each year
comm_2000 <- readRDS("Outputs/2000_community_non_severe_regional_unexp_samples.Rds")
comm_2015 <- readRDS("Outputs/2015_community_non_severe_regional_unexp_samples.Rds")

sevr_2015 <- readRDS("Outputs/2015_community_severe_regional_unexp_samples.Rds")
sevr_2000 <- readRDS("Outputs/2000_community_severe_regional_unexp_samples.Rds")

## Summary statistics
function_list <- list(est = function(x) mean(x), lci = function(x) quantile(x, c(0.025)),
                      uci = function(x) quantile(x, c(0.975)))
inc_ests <- map(list(comm_2000, comm_2015, sevr_2000, sevr_2015), ~
                  map_dbl(function_list, function(f) f(.x)))
inc_ests <- do.call(rbind, inc_ests) %>% 
  round(1) %>%  
  as.tibble()
inc_ests <- inc_ests %>% 
  mutate(type = c("comm", "comm", "sevr", "sevr"),
         year = c(2000, 2015, 2000, 2015))
inc_ests <- inc_ests %>% 
  mutate(res = paste0(est, " (", lci, "-", uci, ")")) %>% 
  select(type, res, year)
# Sample so same number in both sets of iterations
num_iters <-  unique(rapply(sum_pj_dhs4_by_hiv, length))
comm_2000 <- sample(comm_2000, num_iters)
comm_2015 <- sample(comm_2015, num_iters)
sevr_2000 <- sample(sevr_2000, num_iters)
sevr_2015 <- sample(sevr_2015, num_iters)

rates_by_hiv4 <- map(sum_pj_dhs4_by_hiv,
                     function(sumpj4s){
                       rate_dhs4_2000 <- map(sumpj4s, ~ comm_2000 * .x)
                       rate_dhs4_2000_sev <- map(sumpj4s, ~ sevr_2000 * .x)
                       list(comm4 = rate_dhs4_2000,
                            sevr4 = rate_dhs4_2000_sev)
                     })


rates_by_hiv7 <- map(sum_pj_dhs7_by_hiv,
                     function(sumpj7s){
                       rate_dhs7_2015 <- map(sumpj7s, ~ comm_2015 * .x)
                       rate_dhs7_2015_sev <- map(sumpj7s, ~ sevr_2015 * .x)
                       list(comm5 = rate_dhs7_2015,
                            sevr5 = rate_dhs7_2015_sev)
                     })


## Correct for HiB vaccination
## Same Hib proprtion in all states,
## But vaccination varies by state
## Only do for 2015
setdiff(names(rates_by_hiv4$all$comm4), hib_corr$state_code)
setdiff(hib_corr$state_code, names(rates_by_hiv4$all$comm4))
## USes state code consistently, but dont have for all INDIA, set this to meanfor coverage so dont subtract anything
hib_corr <- hib_corr %>% 
  add_row(year = 2015, state = "India", state_code = "IA",hib3_vec_ind_adj = mean(hib_corr$hib3_vec_ind_adj) )
## Convert to a list for comparability and sort the list to match the state list
hib_corr_list <- tapply(hib_corr$hib3_vec_ind_adj, hib_corr$state_code, identity, simplify = FALSE)
hib_corr_list <- hib_corr_list[names(rates_by_hiv7$all$comm5)]

rates_by_hiv7 <- map(rates_by_hiv7,
                     function(hiv_status){
                       dropcomm5 <- map2(hiv_status$comm5, hib_corr_list,
                                         function(indian_state_rate, indian_state_hib_vacc){
                                           indian_state_rate - (indian_state_rate * indian_state_hib_vacc * indian_hib_proportion_severe)
                                           
                                         })
                       dropsevr5 <- map2(hiv_status$sevr5, hib_corr_list,
                                         function(indian_state_rate, indian_state_hib_vacc){
                                           indian_state_rate - (indian_state_rate * indian_state_hib_vacc * indian_hib_proportion_severe)
                                         })
                       list(comm5 = dropcomm5,
                            sevr5 = dropsevr5)
                     })
## Plotted comparison, looks fine
# cmpr_hib <- map(names(hib_corr_list), function(states){
#   quantile(rates_by_hiv7$all$comm5[[states]] - rates_by_hiv7_b$all$comm5[[states]], probs = c(0.01, 0.5, 0.99))
# })
# names(cmpr_hib) <- names(hib_corr_list)
# cmpr_hib

## Bind all together
rates_by_hiv <- map2(rates_by_hiv4, rates_by_hiv7, function(x,y) c(x,y))

## Calculate event numbers
india_pop <- pop %>% 
  group_by(year) %>% 
  summarise(u5_pop = sum(u5_pop)) %>% 
  mutate(state = "India")

pop <- bind_rows(pop, india_pop)

pop <- pop %>%
  rename(state_pop_name = state) %>%
  inner_join(state_lkp %>%  distinct(abbrev, state_pop_name)) %>% 
  filter(abbrev != "")

PopLKP <- function (year) {
  lkp <- pop$u5_pop[pop$year == year]
  names(lkp) <- pop$abbrev[pop$year == year]
  lkp
}
pop_2000 <- PopLKP(2000)
pop_2015 <- PopLKP(2015)

pop_2000 <- pop_2000[names(rates_by_hiv$all$comm4)]
pop_2015 <- pop_2015[names(rates_by_hiv$all$comm5)]

pop_by_hiv <- map2(list(all = 1, hiv = 53403/119369000, no_hiv = 1-(53403/119369000)),
                   list(all = 1, hiv = 46087/116176000, no_hiv = 1-(46087/116176000)),
                   ~ list (pop_2000 = .x * pop_2000,
                           pop_2015 = .y * pop_2015))
# rm(pop_2000, pop_2015)

# Events by hiv
events_by_hiv <- map2(rates_by_hiv, pop_by_hiv, function(rates, pops){
  dhs4_events     <- map2(rates$comm4, pops$pop_2000, ~ .y * .x/1000 )
  dhs7_events     <- map2(rates$comm5, pops$pop_2015, ~ .y * .x/1000 )
  dhs4_events_sev <- map2(rates$sevr4, pops$pop_2000, ~ .y * .x/1000 )
  dhs7_events_sev <- map2(rates$sevr5, pops$pop_2015, ~ .y * .x/1000 )
  
  list(dhs4_events = dhs4_events,
       dhs7_events = dhs7_events,
       dhs4_events_sev = dhs4_events_sev,
       dhs7_events_sev = dhs7_events_sev)
})

## Examine change in rates
rate_compare_by_hiv <- map(rates_by_hiv, function(rates){
  in_both <- intersect(names(rates$comm4), names(rates$comm5))
  comm4 <- rates$comm4[in_both]
  comm5 <- rates$comm5[in_both]
  sevr4 <- rates$sevr4[in_both]
  sevr5 <- rates$sevr5[in_both]
  
  log_rate_compare     <- map2(comm5, comm4, ~ log(.x) - log(.y))
  log_rate_compare_sev <- map2(sevr5, sevr4, ~ log(.x) - log(.y))
  
  list(log_rate_compare = log_rate_compare,
       log_rate_compare_sev = log_rate_compare_sev)
})

## Summarise events and rates
RateSum <- function(x, round_lvl = 0, tform = FALSE) {
  res <- c(mean = mean(x),
           q2_5 = quantile(x, probs = 0.025),
           q97_5 = quantile(x, probs = 0.975))
  if(tform) res <- exp(res)
  a <- round(res, round_lvl)
  tibble(mean = a[1], q2_5 = a[2], q97_5 = a[3])
}

## Summarise all results by all groups - and place in a single dataframe
stat_summaries <- list(
  rates = rapply(rates_by_hiv, RateSum, how="replace"),
  events = rapply(events_by_hiv, RateSum, how="replace", round_lvl = 0),
  relative_change = rapply(rate_compare_by_hiv, RateSum, how = "replace", round_lvl = 2, tform = TRUE ))

# MissingNestedList <- function(x){
#   missing <- rapply(x, function(x) any(is.na(x)))
#   missing[missing]
# }
# MissingNestedList(events_by_hiv)

FlatFunctionWide <- function(x) {
  FlatFunction <- function(x) do.call(c, x)
  
  while("list" %in% class(x[[1]])) {
    x <- FlatFunction(x)
  }
  bind_rows(x, .id = "misc")
}
stat_summaries <- FlatFunctionWide(stat_summaries)

stat_summaries2 <- stat_summaries %>% 
  # mutate(misc = str_replace(misc, fixed("events_sev"), "events.sev")) %>% 
  separate(misc, into = c("stat", "sgrp", "severity", "abbrev"), sep = ("\\."), remove = FALSE) %>% 
  mutate(year = if_else(str_detect(misc, "4"), 2000, 2015),
         severity = if_else(str_detect(severity, "_sev|sevr"), "severe", "community")) %>% 
  select(-misc)
table(stat_summaries2$stat, stat_summaries2$sgrp)
stat_summaries2 %>% filter(abbrev == "AP", stat == "rates") %>%  arrange(year, severity)

# Add in credible interval as single value
smrs_all <- stat_summaries2 %>% 
  mutate(res = paste0(mean, " (", q2_5, "-", q97_5, ")")) %>% 
  select(sgrp, severity, abbrev, year, stat, everything()) %>% 
  arrange(severity, stat, abbrev, year,  sgrp)
# saveRDS(smrs_all, "Outputs/india_results_dataframe.Rds")

# Print as CSV file
# write_csv(smrs_all, "outputs/india_states.csv")


## Aggregate rates for regions
region_abbrevs <- list(NortheastState = c("AR", "MN", "MG", "MZ", "NL", "SK", "TR"),
                       Unionterr   =  c("AN", "CH", "DN", "DD",  "LK",  "PD"),
                       Central = c("CG", "MP", "RJ", "UP"), 
                       East = c("BH", "JH", "OR", "WB"), 
                       North = c("CH", "DL", "HR", "HP", "JM", "PJ", "UC"),
                       Northeast = c("AR", "AS", "MN", "MG", "MZ", "NL", "SK", "TR"),
                       South = c("AN", "AP", "KA", "KE", "LK", "PD", "TN"), 
                       West = c("DN", "DD", "GO", "GJ", "MH"))

# Northeast refers to NortheastRegion
# Unionterr refers to Union territories

region_events <- map(region_abbrevs, function(region) {
  res <- events_by_hiv$all
  res[] <- map2(events_by_hiv$no_hiv, names(events_by_hiv$no_hiv), function(measure, measure_name){
    events <- (measure[names(measure) %in% region]) %>% 
      as.data.frame() %>% 
      rowSums()
    if(str_sub(measure_name, 1, 4) == "dhs4") pop_choose <-  pop_2000 else pop_choose <- pop_2015 
    rates <- events / sum(pop_choose[region])
    tibble(events = events, rates = rates)
  })
  res 
})
# saveRDS(region_events, file = "Outputs/Regional rates and events.Rds")

region_events_hiv <- map(region_abbrevs, function(region) {
  res <- events_by_hiv$all
  res[] <- map2(events_by_hiv$hiv, names(events_by_hiv$hiv), function(measure, measure_name){
    events <- (measure[names(measure) %in% region]) %>% 
      as.data.frame() %>% 
      rowSums()
    if(str_sub(measure_name, 1, 4) == "dhs4") pop_choose <-  pop_by_hiv$hiv$pop_2000 else pop_choose <- pop_by_hiv$hiv$pop_2015 
    rates <- events / sum(pop_choose[region])
    tibble(events = events, rates = rates)
  })
  res 
})

QuickSmry <- 
  list(
    mean = function(x) mean(x, na.rm = T),
    q2_5 = function(x) quantile(x, 0.025, na.rm = T),
    q97_5 = function(x) quantile(x, 0.975, na.rm = T)
  )

SummariseRegions <- function(region_events){
  smrs_regs <- map(QuickSmry, ~ rapply(region_events, .x, how = "list"))
  smrs_regs <- FlatFunctionWide(smrs_regs) %>%  t()
  smrs_regs <- smrs_regs %>% 
    as.data.frame() %>% 
    mutate(for_split = rownames(smrs_regs)) %>% 
    as_tibble()
  smrs_regs2 <- smrs_regs %>% 
    mutate(severe_pneum = str_detect(for_split, "events_sev")) %>% 
    separate(for_split, into = c("summary_statistic", "region", "year_survey", "event_or_rate"), 
             sep = ("\\.")) %>% 
    filter(!is.na(region)) %>% 
    mutate(year_survey = str_sub(year_survey, 1, 4)) %>% 
    rename(value = V1) %>% 
    arrange(region, year_survey, event_or_rate, severe_pneum, summary_statistic) %>% 
    spread(summary_statistic, value)
  
  smrs_regs3 <- smrs_regs2 %>% 
    rename(stat = event_or_rate) %>% 
    mutate_at(vars(mean, q2_5, q97_5), function(x) x %>%  as.character() %>%  as.double()) %>% 
    mutate(mean = if_else(stat == "events", round(mean/1000, -1), round(mean*1000)),
           q2_5  = if_else(stat == "events", round(q2_5/1000, -1), round(q2_5*1000)),
           q97_5  = if_else(stat == "events", round(q97_5/1000, -1), round(q97_5*1000)),
           res2 = paste0(mean, " (", q2_5, "-", q97_5, ")" ))
  
  smrs_regs4 <- smrs_regs3 %>% 
    rename(name = region, severity = severe_pneum, year = year_survey) %>% 
    mutate(abbrev = NA,
           severity = if_else(severity, "severe", "community"),
           year = if_else(year == "dhs4", 2000, 2015),
           res = res2) %>% 
    filter(stat %in% c("events", "rates")) %>% 
    unite(spread_by, year, severity, stat) %>% 
    select(spread_by, name, res2) %>% 
    spread(spread_by, res2)
  smrs_regs4
}
smrs_regs4 <- SummariseRegions(region_events)
smrs_regs4_hiv <- SummariseRegions(region_events_hiv)
# saveRDS(smrs_regs4, "Outputs/regional_summarise.Rds")
saveRDS(smrs_regs4_hiv, "Outputs/regional_summarise_hiv.Rds")
## repeat for hiv

# write_csv(smrs_regs2, "Outputs/India_regions.csv")
# region_rates <- region_events
# region_rates[] <- rapply(region_events)


# rm(events_by_hiv, pop_by_hiv, rate_compare_by_hiv, rates_by_hiv)

# Summarise marginal proportions
# Unalbe to share this next file due to DHS agreement
load("../../DemographicHealthSurvey/R extract/India/region aggregated data.Rdata")

varlist <- list("lbw",  "inc.imm", "non.ebf", "indoor.air", "overcrowd", "maln.long", "hiv")

final_abbrevs <- final.agg %>% 
  mutate(dhs = str_sub(ctr.type.2digit, 4, 4) %>%  as.integer()) %>% 
  mutate(code = as.integer(hv024)) %>% 
  select(-hv024, -ctr.type.2digit) %>%  
  inner_join(state_lkp %>%  select(dhs, code, abbrev))

margs <- map(varlist, function(var_choose){
  fa2 <- final_abbrevs %>%
    group_by_("dhs", "abbrev",var_choose) %>% 
    summarise(n = sum(n))  %>% 
    na.omit() %>% 
    mutate(p = n/sum(n)) %>% 
    select_("dhs", "abbrev", var_choose, "p") %>% 
    ungroup()
  fa2 <- fa2 %>% as.data.frame()
  fa2[fa2[, var_choose] ==1, ]
})

names(margs) <- varlist
margs <- map2(margs, varlist, function (x, var_choose) x[, !names(x) %in% var_choose])
margs <- bind_rows(margs, .id = "var")

margs <- margs %>% 
  spread(key = var, value = p)

margs <- margs %>% 
  inner_join(final_abbrevs %>%
               group_by(dhs, abbrev) %>% 
               summarise(n = sum(n),
                         year = mean(year))) 

margs_dhs_earliest <- margs %>%
  filter(dhs == 4) %>% 
  inner_join(state_lkp %>% 
               filter(dhs %in% c(0,4)) %>% 
               select(abbrev, state_pop_name))

margs_dhs_latest <- margs %>%
  filter(dhs == 7) %>% 
  inner_join(state_lkp %>% 
               filter(dhs %in% c(0,7)) %>% 
               select(abbrev, state_pop_name))

## Aggregate for India from each State
IndiaFunction <- function(x, margs_dhs_timepoint, pop){
  margs_dhs_timepoint$pop <- pop[margs_dhs_timepoint$abbrev]
  margs_dhs_timepoint <- margs_dhs_timepoint %>% 
    filter(!is.na(pop))
  margs_var <- margs_dhs_timepoint[, x]
  
  b <- weighted.mean(margs_var, margs_dhs_timepoint$pop, na.rm = TRUE)
  b
}
varlist_india4 <- map(varlist, ~ IndiaFunction(.x, margs_dhs_earliest, pop_2000))
varlist_india7 <- map(varlist, ~ IndiaFunction(.x, margs_dhs_latest,   pop_2015))
varlist_india4 <- bind_cols(varlist_india4) %>% 
  set_names(varlist)
varlist_india7 <- bind_cols(varlist_india7) %>% 
  set_names(varlist)

margs_dhs_earliest <- bind_rows(margs_dhs_earliest,
                                varlist_india4 %>% 
                                  mutate(dhs = 4, abbrev = "IA", state_pop_name = "India",
                                         year = 1999))

margs_dhs_latest <- bind_rows(margs_dhs_latest,
                              varlist_india7 %>% 
                                mutate(dhs = 7, abbrev = "IA", state_pop_name = "India",
                                       year = 2015))


margs <- bind_rows(margs_dhs_earliest, margs_dhs_latest)
print(getwd())
# write_csv(margs, path = paste0("Outputs/risk_factor_data_indian_states.csv"), na = "")

