#------------------
# posterior_summary
#------------------
library(dplyr)

## Read Stan model output (rds file)
results_list <- readRDS("stanfits/all_results_list.rds")

summary_list <- list()
rhat_list    <- list()

for (name in names(results_list)) {
  res <- results_list[[name]]
  
  if (is.list(res) && "posterior" %in% names(res)) {
    post_df <- as.data.frame(res$posterior)} 
  
  fit_obj <- NULL
  if ("fit" %in% names(res)) {
    fit_obj <- res$fit
  }
  
  params <- intersect(c("L","k","x0","phi"), colnames(post_df))
  param_summaries <- lapply(params, function(p) {
    vals <- post_df[[p]]
    data.frame(
      Model     = name,
      Parameter = p,
      Mean      = mean(vals),
      SD        = sd(vals),
      Q2.5      = quantile(vals, 0.025),
      Median    = quantile(vals, 0.5),
      Q97.5     = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    )
  })
  summary_list[[name]] <- bind_rows(param_summaries)
  
  if (!is.null(fit_obj)) {
    fit_sum <- summary(fit_obj)$summary
    rhat_list[[name]] <- data.frame(
      Model     = name,
      Parameter = rownames(fit_sum),
      Rhat      = fit_sum[,"Rhat"],
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- bind_rows(summary_list)
write.csv(summary_df, "results/posterior_summary.csv", row.names = FALSE)

if (length(rhat_list)>0) {
  rhat_df <- bind_rows(rhat_list)
  write.csv(rhat_df, "results/rhat_values.csv", row.names = FALSE)
} 


#-------------------------
# PI_and_reduction_summary 
#-------------------------
library(dplyr)
library(readxl)
library(rstan)
library(tibble)

res_file <- file.path("stanfits/all_results_list.rds")
if (!file.exists(res_file)) stop("Results file not found")
results_list <- readRDS(res_file)

vs <- c(0, 0.4, 0.5215, 0.6)

# Compute PI(v) for each subgroup and each v
PI_summary <- tibble()
for (grp in names(results_list)) {
  post     <- results_list[[grp]]$posterior
  params_df <- data.frame(L = post$L, k = post$k, x0 = post$x0)
  
  for (v in vs) {
    PIv <- apply(params_df, 1, function(params) {
      L  <- params["L"]
      k  <- params["k"]
      x0 <- params["x0"]
      L / (1 + exp(-k * (v - x0)))
    })
    
    PI_summary <- bind_rows(
      PI_summary,
      tibble(
        group  = grp,
        v      = v,
        mean   = mean(PIv),
        lower  = quantile(PIv, 0.025),
        upper  = quantile(PIv, 0.975),
        metric = "PI"
      )
    )
  }
}

# Compute relative reduction for each v in the two age ranges
reduction_summary <- tibble()
for (v in vs) {
  # draws at this v
  PI_fv_16_49  <- apply(data.frame(
    L = results_list[["Rate(16=<V<50)"]]$posterior$L,
    k = results_list[["Rate(16=<V<50)"]]$posterior$k,
    x0 = results_list[["Rate(16=<V<50)"]]$posterior$x0
  ), 1, function(p) p["L"]/(1+exp(-p["k"]*(v-p["x0"]))))
  PI_nfv_16_49 <- apply(data.frame(
    L = results_list[["Rate(16=<UV<50)"]]$posterior$L,
    k = results_list[["Rate(16=<UV<50)"]]$posterior$k,
    x0 = results_list[["Rate(16=<UV<50)"]]$posterior$x0
  ), 1, function(p) p["L"]/(1+exp(-p["k"]*(v-p["x0"]))))
  ratio_16_49  <- 1 - PI_fv_16_49 / PI_nfv_16_49
  
  PI_fv_50p    <- apply(data.frame(
    L = results_list[["Rate(V>=50)"]]$posterior$L,
    k = results_list[["Rate(V>=50)"]]$posterior$k,
    x0 = results_list[["Rate(V>=50)"]]$posterior$x0
  ), 1, function(p) p["L"]/(1+exp(-p["k"]*(v-p["x0"]))))
  PI_nfv_50p   <- apply(data.frame(
    L = results_list[["Rate(UV>=50)"]]$posterior$L,
    k = results_list[["Rate(UV>=50)"]]$posterior$k,
    x0 = results_list[["Rate(UV>=50)"]]$posterior$x0
  ), 1, function(p) p["L"]/(1+exp(-p["k"]*(v-p["x0"]))))
  ratio_50_plus <- 1 - PI_fv_50p / PI_nfv_50p
  
  reduction_summary <- bind_rows(
    reduction_summary,
    tibble(
      group  = "16-49 (FV vs NFV)",
      v      = v,
      mean   = mean(ratio_16_49),
      lower  = quantile(ratio_16_49, 0.025),
      upper  = quantile(ratio_16_49, 0.975),
      metric = "relative_reduction"
    ),
    tibble(
      group  = "50+ (FV vs NFV)",
      v      = v,
      mean   = mean(ratio_50_plus),
      lower  = quantile(ratio_50_plus, 0.025),
      upper  = quantile(ratio_50_plus, 0.975),
      metric = "relative_reduction"
    )
  )
}

# Combine both
combined_summary <- bind_rows(PI_summary, reduction_summary)

write.csv(
  combined_summary,
  file.path("results/combined_PI_and_reduction_summary.csv"),
  row.names = FALSE
)



#-------------------------------------------------
# Averted Infections, with Direct/Indirect averted
#-------------------------------------------------
library(readxl)
library(dplyr)
library(tibble)

res_file <- file.path("stanfits/all_results_list.rds")
results_list    <- readRDS(res_file)

# Population sizes and observed counts
group_sizes <- c(
  "Rate(Age<12)"     = 940404,
  "Rate(12<=Age<16)" = 307473,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(V>=50)"      = 1663770,
  "Rate(UV>=50)"     = 509900
)
actual_numbers <- c(
  "Rate(Age<12)"     = 16348,
  "Rate(12<=Age<16)" = 4722,
  "Rate(16=<UV<50)"  = 34142,
  "Rate(16=<V<50)"   = 8538,
  "Rate(UV>=50)"     = 9217,
  "Rate(V>=50)"      = 6063
)

nfv_groups <- c("Rate(Age<12)", "Rate(12<=Age<16)",
                "Rate(16=<UV<50)", "Rate(UV>=50)")
nf_PI <- lapply(nfv_groups, function(col) {
  post <- results_list[[col]]$posterior
  post$L / (1 + exp(-post$k * (0 - post$x0)))
})
names(nf_PI) <- nfv_groups

fv_match <- list(
  "Rate(16=<V<50)" = "Rate(16=<UV<50)",
  "Rate(V>=50)"    = "Rate(UV>=50)"
)

group_PI <- list()
for (col in names(results_list)) {
  post <- results_list[[col]]$posterior
  PI0  <- post$L / (1 + exp(-post$k * (0 - post$x0)))
  
  if (col %in% nfv_groups) {
    group_PI[[col]] <- PI0
  } else {
    key <- fv_match[[col]]
    group_PI[[col]] <- nf_PI[[key]]
  }
}

averted_summary <- tibble()
overall_total_draws  <- NULL
overall_direct_draws <- NULL

for (col in names(group_PI)) {
  n         <- group_sizes[col]
  exp_risk  <- group_PI[[col]] * n
  obs       <- actual_numbers[col]
  
  # total averted = expected − observed
  total_av  <- exp_risk - obs
  
  # direct averted only for FV groups
  if (col %in% names(fv_match)) {
    key       <- fv_match[[col]]
    fv_draws  <- results_list[[col]]$posterior$L / 
      (1 + exp(-results_list[[col]]$posterior$k * 
                 (0 - results_list[[col]]$posterior$x0)))
    direct_av <- (nf_PI[[key]] - fv_draws) * n
  } else {
    direct_av <- rep(0, length(total_av))
  }
  
  indirect_av <- total_av - direct_av
  
  averted_summary <- bind_rows(averted_summary, tibble(
    group          = col,
    total_mean     = as.integer(round(mean(total_av))),
    total_lower    = as.integer(round(quantile(total_av, 0.025))),
    total_upper    = as.integer(round(quantile(total_av, 0.975))),
    direct_mean    = as.integer(round(mean(direct_av))),
    direct_lower   = as.integer(round(quantile(direct_av, 0.025))),
    direct_upper   = as.integer(round(quantile(direct_av, 0.975))),
    indirect_mean  = as.integer(round(mean(indirect_av))),
    indirect_lower = as.integer(round(quantile(indirect_av, 0.025))),
    indirect_upper = as.integer(round(quantile(indirect_av, 0.975)))
  ))
  
  overall_total_draws  <- if (is.null(overall_total_draws))  total_av  else overall_total_draws  + total_av
  overall_direct_draws <- if (is.null(overall_direct_draws)) direct_av else overall_direct_draws + direct_av
}

overall_indirect_draws <- overall_total_draws - overall_direct_draws
overall_summary <- tibble(
  group           = "Overall",
  total_mean      = as.integer(round(mean(overall_total_draws))),
  total_lower     = as.integer(round(quantile(overall_total_draws, 0.025))),
  total_upper     = as.integer(round(quantile(overall_total_draws, 0.975))),
  direct_mean     = as.integer(round(mean(overall_direct_draws))),
  direct_lower    = as.integer(round(quantile(overall_direct_draws, 0.025))),
  direct_upper    = as.integer(round(quantile(overall_direct_draws, 0.975))),
  indirect_mean   = as.integer(round(mean(overall_indirect_draws))),
  indirect_lower  = as.integer(round(quantile(overall_indirect_draws, 0.025))),
  indirect_upper  = as.integer(round(quantile(overall_indirect_draws, 0.975)))
)

all_averted <- bind_rows(averted_summary, overall_summary)

write.csv(
  all_averted,
  file.path("results/averted_infections_summary.csv"),
  row.names = FALSE
)



#------------------------------------------------------
# Averted Hospitalisation, with Direct/Indirect averted
#------------------------------------------------------
library(readxl)
library(dplyr)
library(tibble)
library(tidyr)

res_file <- file.path("stanfits/all_results_list.rds")
results_list <- readRDS(res_file)

groups        <- names(results_list)
unvax_groups  <- c("Rate(Age<12)", "Rate(12<=Age<16)",
                   "Rate(16=<UV<50)", "Rate(UV>=50)")
map_unvax     <- c("Rate(16=<V<50)" = "Rate(16=<UV<50)",
                   "Rate(V>=50)"    = "Rate(UV>=50)")

group_sizes   <- c(
  "Rate(Age<12)"     =  940404,
  "Rate(12<=Age<16)" =  307473,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(UV>=50)"     =  509900,
  "Rate(V>=50)"      = 1663770
)

infection_hosp <- c(
  "Rate(Age<12)"     = 107/16348,
  "Rate(12<=Age<16)" =  32/4722,
  "Rate(16=<UV<50)"  =1581/34142,
  "Rate(UV>=50)"     =2064/9217,
  "Rate(16=<V<50)"   = 102/8538,
  "Rate(V>=50)"      = 596/6063
)

actual_hosp <- c(
  "Rate(Age<12)"     = 107,
  "Rate(12<=Age<16)" =  32,
  "Rate(16=<UV<50)"  =1581,
  "Rate(UV>=50)"     =2064,
  "Rate(16=<V<50)"   = 102,
  "Rate(V>=50)"      = 596
)

nf_PI <- lapply(unvax_groups, function(grp) {
  post <- results_list[[grp]]$posterior
  post$L / (1 + exp(-post$k * (0 - post$x0)))
})
names(nf_PI) <- unvax_groups

fv_PI    <- list()
group_PI <- list()

for (grp in groups) {
  post <- results_list[[grp]]$posterior
  PI0  <- post$L / (1 + exp(-post$k * (0 - post$x0)))
  
  if (grp %in% unvax_groups) {
    group_PI[[grp]] <- PI0
  } else {
    fv_PI[[grp]]    <- PI0
    base_grp        <- map_unvax[grp]
    group_PI[[grp]] <- nf_PI[[base_grp]]
  }
}

averted_summary        <- tibble()
overall_total_draws    <- NULL
overall_direct_draws   <- NULL
overall_indirect_draws <- NULL

for (grp in groups) {
  n    <- group_sizes[grp]
  IHR0 <- infection_hosp[ if (grp %in% unvax_groups) grp else map_unvax[grp] ]
  IHRi <- infection_hosp[grp]
  obs  <- actual_hosp[grp]
  
  H0       <- group_PI[[grp]] * n * IHR0
  delta_H1 <- group_PI[[grp]] * n * (IHR0 - IHRi)
  delta_I  <- if (grp %in% names(map_unvax)) {
    (nf_PI[[ map_unvax[grp] ]] - fv_PI[[grp]]) * n
  } else {
    rep(0, length(H0))
  }
  delta_H2 <- delta_I * IHRi
  direct   <- delta_H1 + delta_H2
  
  total_av <- H0 - obs
  indirect <- total_av - direct
  
  averted_summary <- bind_rows(averted_summary, tibble(
    group          = grp,
    total_mean     = as.integer(round(mean(total_av))),
    total_lower    = as.integer(round(quantile(total_av, 0.025))),
    total_upper    = as.integer(round(quantile(total_av, 0.975))),
    direct_mean    = as.integer(round(mean(direct))),
    direct_lower   = as.integer(round(quantile(direct, 0.025))),
    direct_upper   = as.integer(round(quantile(direct, 0.975))),
    indirect_mean  = as.integer(round(mean(indirect))),
    indirect_lower = as.integer(round(quantile(indirect, 0.025))),
    indirect_upper = as.integer(round(quantile(indirect, 0.975)))
  ))
  
  overall_total_draws    <- if (is.null(overall_total_draws))    total_av    else overall_total_draws    + total_av
  overall_direct_draws   <- if (is.null(overall_direct_draws))   direct      else overall_direct_draws   + direct
  overall_indirect_draws <- if (is.null(overall_indirect_draws)) indirect    else overall_indirect_draws + indirect
}

overall_indirect_draws <- overall_total_draws - overall_direct_draws
overall_df <- tibble(
  group          = "Overall",
  total_mean     = as.integer(round(mean(overall_total_draws))),
  total_lower    = as.integer(round(quantile(overall_total_draws, 0.025))),
  total_upper    = as.integer(round(quantile(overall_total_draws, 0.975))),
  direct_mean    = as.integer(round(mean(overall_direct_draws))),
  direct_lower   = as.integer(round(quantile(overall_direct_draws, 0.025))),
  direct_upper   = as.integer(round(quantile(overall_direct_draws, 0.975))),
  indirect_mean  = as.integer(round(mean(overall_indirect_draws))),
  indirect_lower = as.integer(round(quantile(overall_indirect_draws, 0.025))),
  indirect_upper = as.integer(round(quantile(overall_indirect_draws, 0.975)))
)

averted_long <- averted_summary %>%
  pivot_longer(
    cols      = -group,
    names_to  = c("measure", ".value"),
    names_sep = "_"
  )

overall_long <- overall_df %>%
  select(group, everything())

all_averted <- bind_rows(averted_long, overall_long)

write.csv(
  all_averted,
  file.path("results/averted_hospitalisations_summary.csv"),
  row.names = FALSE
)




#--------------------------------------------
# Averted Deaths, and Direct/Indirect averted
#--------------------------------------------
library(readxl)
library(dplyr)
library(tibble)
library(tidyr)

res_file     <- "stanfits/all_results_list.rds"
results_list <- readRDS(res_file)

groups        <- names(results_list)
unvax_groups  <- c("Rate(Age<12)", "Rate(12<=Age<16)",
                   "Rate(16=<UV<50)", "Rate(UV>=50)")
map_unvax     <- c("Rate(16=<V<50)" = "Rate(16=<UV<50)",
                   "Rate(V>=50)"    = "Rate(UV>=50)")

group_sizes   <- c(
  "Rate(Age<12)"     =  940404,
  "Rate(12<=Age<16)" =  307473,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(UV>=50)"     =  509900,
  "Rate(V>=50)"      = 1663770
)

infection_dea <- c(
  "Rate(Age<12)"     = 1/16348,
  "Rate(12<=Age<16)" = 1/4722,
  "Rate(16=<UV<50)"  = 24/34142,
  "Rate(UV>=50)"     = 391/9217,
  "Rate(16=<V<50)"   = 2/8538,
  "Rate(V>=50)"      = 195/6063
)

actual_dea <- c(
  "Rate(Age<12)"     = 1,
  "Rate(12<=Age<16)" = 1,
  "Rate(16=<UV<50)"  = 24,
  "Rate(UV>=50)"     = 391,
  "Rate(16=<V<50)"   = 2,
  "Rate(V>=50)"      = 195
)

nf_PI <- lapply(unvax_groups, function(grp) {
  post <- results_list[[grp]]$posterior
  post$L / (1 + exp(-post$k * (0 - post$x0)))
})
names(nf_PI) <- unvax_groups

fv_PI    <- list()
group_PI <- list()

for (grp in groups) {
  post <- results_list[[grp]]$posterior
  PI0  <- post$L / (1 + exp(-post$k * (0 - post$x0)))
  
  if (grp %in% unvax_groups) {
    group_PI[[grp]] <- PI0
  } else {
    fv_PI[[grp]]    <- PI0
    base_grp        <- map_unvax[grp]
    group_PI[[grp]] <- nf_PI[[base_grp]]
  }
}

averted_summary        <- tibble()
overall_total_draws    <- NULL
overall_direct_draws   <- NULL
overall_indirect_draws <- NULL

for (grp in groups) {
  n    <- group_sizes[grp]
  IFR0 <- infection_dea[ if (grp %in% unvax_groups) grp else map_unvax[grp] ]
  IFRi <- infection_dea[grp]
  obs  <- actual_dea[grp]
  
  D0       <- group_PI[[grp]] * n * IFR0
  delta_D1 <- group_PI[[grp]] * n * (IFR0 - IFRi)
  delta_I  <- if (grp %in% names(map_unvax)) {
    (nf_PI[[ map_unvax[grp] ]] - fv_PI[[grp]]) * n
  } else {
    rep(0, length(D0))
  }
  delta_D2 <- delta_I * IFRi
  direct   <- delta_D1 + delta_D2
  
  total_av <- D0 - obs
  indirect <- total_av - direct
  
  averted_summary <- bind_rows(averted_summary, tibble(
    group          = grp,
    total_mean     = as.integer(round(mean(total_av))),
    total_lower    = as.integer(round(quantile(total_av, 0.025))),
    total_upper    = as.integer(round(quantile(total_av, 0.975))),
    direct_mean    = as.integer(round(mean(direct))),
    direct_lower   = as.integer(round(quantile(direct, 0.025))),
    direct_upper   = as.integer(round(quantile(direct, 0.975))),
    indirect_mean  = as.integer(round(mean(indirect))),
    indirect_lower = as.integer(round(quantile(indirect, 0.025))),
    indirect_upper = as.integer(round(quantile(indirect, 0.975)))
  ))
  
  overall_total_draws    <- if (is.null(overall_total_draws))    total_av    else overall_total_draws    + total_av
  overall_direct_draws   <- if (is.null(overall_direct_draws))   direct      else overall_direct_draws   + direct
  overall_indirect_draws <- if (is.null(overall_indirect_draws)) indirect    else overall_indirect_draws + indirect
}

overall_indirect_draws <- overall_total_draws - overall_direct_draws
overall_df <- tibble(
  group          = "Overall",
  total_mean     = as.integer(round(mean(overall_total_draws))),
  total_lower    = as.integer(round(quantile(overall_total_draws, 0.025))),
  total_upper    = as.integer(round(quantile(overall_total_draws, 0.975))),
  direct_mean    = as.integer(round(mean(overall_direct_draws))),
  direct_lower   = as.integer(round(quantile(overall_direct_draws, 0.025))),
  direct_upper   = as.integer(round(quantile(overall_direct_draws, 0.975))),
  indirect_mean  = as.integer(round(mean(overall_indirect_draws))),
  indirect_lower = as.integer(round(quantile(overall_indirect_draws, 0.025))),
  indirect_upper = as.integer(round(quantile(overall_indirect_draws, 0.975)))
)

averted_long <- averted_summary %>%
  pivot_longer(
    cols      = -group,
    names_to  = c("measure", ".value"),
    names_sep = "_"
  )

overall_long <- overall_df %>%
  select(group, everything())

all_averted <- bind_rows(averted_long, overall_long)

write.csv(
  all_averted,
  file.path("results/averted_deaths_summary.csv"),
  row.names = FALSE
)





#----------------------------------------------
# Avertable Infections (Estimate at V = 52.15%)
#----------------------------------------------
res_file     <- "stanfits/all_results_list.rds"
results_list <- readRDS(res_file)
groups       <- names(results_list)

v_target            <- 0.5215
group_sizes         <- c(
  "Rate(Age<12)"     =  940404,
  "Rate(12<=Age<16)" =  307473,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(UV>=50)"     =  509900,
  "Rate(V>=50)"      = 1663770
)
observed_infections <- c(
  "Rate(Age<12)"     = 16348,
  "Rate(12<=Age<16)" =  4722,
  "Rate(16=<UV<50)"  = 34142,
  "Rate(16=<V<50)"   =  8538,
  "Rate(UV>=50)"     =  9217,
  "Rate(V>=50)"      =  6063
)

get_stats <- function(x) {
  ci <- quantile(x, c(0.025, 0.975))
  tibble(
    mean  = mean(x),
    lower = unname(ci[1]),
    upper = unname(ci[2])
  )
}

group_dfs <- lapply(groups, function(grp) {
  post        <- results_list[[grp]]$posterior
  pi_draws    <- post$L / (1 + exp(-post$k * (v_target - post$x0)))
  est_draws   <- pi_draws * group_sizes[grp]
  ave_draws   <- observed_infections[grp] - est_draws
  
  pi_stats   <- get_stats(pi_draws)
  est_stats  <- get_stats(est_draws)
  ave_stats  <- get_stats(ave_draws)
  
  tibble(
    group            = grp,
    PI_mean          = round(pi_stats$mean,   4),
    PI_lower         = round(pi_stats$lower,  4),
    PI_upper         = round(pi_stats$upper,  4),
    estimated_mean   = as.integer(round(est_stats$mean)),
    estimated_lower  = as.integer(round(est_stats$lower)),
    estimated_upper  = as.integer(round(est_stats$upper)),
    avertable_mean   = as.integer(round(ave_stats$mean)),
    avertable_lower  = as.integer(round(ave_stats$lower)),
    avertable_upper  = as.integer(round(ave_stats$upper))
  )
})

group_summary <- bind_rows(group_dfs)

all_est_draws <- Reduce(`+`, lapply(groups, function(grp) {
  post      <- results_list[[grp]]$posterior
  pi_draws  <- post$L / (1 + exp(-post$k * (v_target - post$x0)))
  pi_draws * group_sizes[grp]
}))
all_obs       <- sum(observed_infections)
all_ave_draws <- all_obs - all_est_draws

tot_est_stats <- get_stats(all_est_draws)
tot_ave_stats <- get_stats(all_ave_draws)

total_row <- tibble(
  group            = "Total",
  PI_mean          = NA_real_,
  PI_lower         = NA_real_,
  PI_upper         = NA_real_,
  estimated_mean   = as.integer(round(tot_est_stats$mean)),
  estimated_lower  = as.integer(round(tot_est_stats$lower)),
  estimated_upper  = as.integer(round(tot_est_stats$upper)),
  avertable_mean   = as.integer(round(tot_ave_stats$mean)),
  avertable_lower  = as.integer(round(tot_ave_stats$lower)),
  avertable_upper  = as.integer(round(tot_ave_stats$upper))
)

final_df <- bind_rows(group_summary, total_row)

write.csv(
  final_df,
  file.path("results/avertable_infections_summary.csv"),
  row.names = FALSE
)



#----------------------------------------------------
# Avertable Hospitalisations (Estimate at V = 52.15%)
#----------------------------------------------------
res_file     <- "stanfits/all_results_list.rds"
results_list <- readRDS(res_file)
groups       <- names(results_list)

v_target            <- 0.5215

group_sizes         <- c(
  "Rate(Age<12)"     =  940404,
  "Rate(12<=Age<16)" =  307473,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(UV>=50)"     =  509900,
  "Rate(V>=50)"      = 1663770
)
infection_hosp_rate <- c(
  "Rate(Age<12)"     = 107/16348,
  "Rate(12<=Age<16)" =  32/4722,
  "Rate(16=<UV<50)"  = 1581/34142,
  "Rate(16=<V<50)"   =  102/8538,
  "Rate(UV>=50)"     = 2064/9217,
  "Rate(V>=50)"      =  596/6063
)
observed_hosp <- c(
  "Rate(Age<12)"     =  107,
  "Rate(12<=Age<16)" =   32,
  "Rate(16=<UV<50)"  = 1581,
  "Rate(16=<V<50)"   =  102,
  "Rate(UV>=50)"     = 2064,
  "Rate(V>=50)"      =  596
)

summary_stats <- function(x) {
  ci <- quantile(x, c(0.025, 0.975))
  tibble(
    mean  = mean(x),
    lower = unname(ci[1]),
    upper = unname(ci[2])
  )
}

PI_target_list <- lapply(groups, function(grp) {
  post <- results_list[[grp]]$posterior
  post$L / (1 + exp(-post$k * (v_target - post$x0)))
})
names(PI_target_list) <- groups

hosp_scaled_list <- mapply(function(pi_draws, grp) {
  pi_draws * group_sizes[grp] * infection_hosp_rate[grp]
}, PI_target_list, groups, SIMPLIFY = FALSE)

group_summary <- bind_rows(lapply(groups, function(grp) {
  hs    <- hosp_scaled_list[[grp]]
  st_h  <- summary_stats(hs)
  av    <- observed_hosp[grp] - hs
  st_av <- summary_stats(av)
  
  tibble(
    group            = grp,
    mean             = as.integer(round(st_h["mean"],  0)),
    lower            = as.integer(round(st_h["lower"], 0)),
    upper            = as.integer(round(st_h["upper"], 0)),
    avertable_mean   = as.integer(round(st_av["mean"],  0)),
    avertable_lower  = as.integer(round(st_av["lower"], 0)),
    avertable_upper  = as.integer(round(st_av["upper"], 0))
  )
}))

all_hs     <- Reduce(`+`, hosp_scaled_list)
st_tot_h   <- summary_stats(all_hs)
all_av     <- sum(observed_hosp) - all_hs
st_tot_av  <- summary_stats(all_av)

total_row <- tibble(
  group            = "Total",
  mean             = as.integer(round(st_tot_h["mean"],  0)),
  lower            = as.integer(round(st_tot_h["lower"], 0)),
  upper            = as.integer(round(st_tot_h["upper"], 0)),
  avertable_mean   = as.integer(round(st_tot_av["mean"],  0)),
  avertable_lower  = as.integer(round(st_tot_av["lower"], 0)),
  avertable_upper  = as.integer(round(st_tot_av["upper"], 0))
)

final_df  <- bind_rows(group_summary, total_row)

write.csv(
  final_df,
  file.path("results/avertable_hospitalisations_summary.csv"),
  row.names = FALSE
)



#------------------------------------------
# Avertable Deaths (Estimate at V = 52.15%)
#------------------------------------------
res_file     <- "stanfits/all_results_list.rds"
results_list <- readRDS(res_file)
groups       <- names(results_list)

v_target            <- 0.5215

group_sizes         <- c(
  "Rate(Age<12)"     =  940404,
  "Rate(12<=Age<16)" =  307473,
  "Rate(16=<UV<50)"  = 1323856,
  "Rate(16=<V<50)"   = 1694560,
  "Rate(UV>=50)"     =  509900,
  "Rate(V>=50)"      = 1663770
)
infection_dea_rate <- c(
  "Rate(Age<12)"     = 1/16348,
  "Rate(12<=Age<16)" = 1/4722,
  "Rate(16=<UV<50)"  = 24/34142,
  "Rate(16=<V<50)"   = 2/8538,
  "Rate(UV>=50)"     = 391/9217,
  "Rate(V>=50)"      = 195/6063
)
observed_dea <- c(
  "Rate(Age<12)"     =  1,
  "Rate(12<=Age<16)" =  1,
  "Rate(16=<UV<50)"  =  24,
  "Rate(16=<V<50)"   =  2,
  "Rate(UV>=50)"     =  391,
  "Rate(V>=50)"      =  195
)

summary_stats <- function(x) {
  ci <- quantile(x, c(0.025, 0.975))
  tibble(
    mean  = mean(x),
    lower = unname(ci[1]),
    upper = unname(ci[2])
  )
}

PI_target_list <- lapply(groups, function(grp) {
  post <- results_list[[grp]]$posterior
  post$L / (1 + exp(-post$k * (v_target - post$x0)))
})

names(PI_target_list) <- groups

dea_scaled_list <- mapply(function(pi_draws, grp) {
  pi_draws * group_sizes[grp] * infection_dea_rate[grp]
}, PI_target_list, groups, SIMPLIFY = FALSE)

group_summary <- bind_rows(lapply(groups, function(grp) {
  ds    <- dea_scaled_list[[grp]]
  st_d  <- summary_stats(ds)
  av    <- observed_dea[grp] - ds
  st_av <- summary_stats(av)
  
  tibble(
    group            = grp,
    mean             = as.integer(round(st_d["mean"],  0)),
    lower            = as.integer(round(st_d["lower"], 0)),
    upper            = as.integer(round(st_d["upper"], 0)),
    avertable_mean   = as.integer(round(st_av["mean"],  0)),
    avertable_lower  = as.integer(round(st_av["lower"], 0)),
    avertable_upper  = as.integer(round(st_av["upper"], 0))
  )
}))

all_ds     <- Reduce(`+`, dea_scaled_list)
st_tot_d   <- summary_stats(all_ds)
all_av     <- sum(observed_dea) - all_ds
st_tot_av  <- summary_stats(all_av)

total_row <- tibble(
  group            = "Total",
  mean             = as.integer(round(st_tot_d["mean"],  0)),
  lower            = as.integer(round(st_tot_d["lower"], 0)),
  upper            = as.integer(round(st_tot_d["upper"], 0)),
  avertable_mean   = as.integer(round(st_tot_av["mean"],  0)),
  avertable_lower  = as.integer(round(st_tot_av["lower"], 0)),
  avertable_upper  = as.integer(round(st_tot_av["upper"], 0))
)

final_df  <- bind_rows(group_summary, total_row)

write.csv(
  final_df,
  file.path("results/avertable_deaths_summary.csv"),
  row.names = FALSE
)

