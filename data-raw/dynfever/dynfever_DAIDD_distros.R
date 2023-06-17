
library(ICI3D)
library(data.table)

sampn <- 500

system.time(dynfever_DAIDD_distros <-
  seq(0, 1, by = 0.05) |> lapply(\(coverage) {
  rbind(
    (dynfever_sample(
      n = sampn, parms = dynfever_params(vax_dogs = coverage),
      ref.seed = 1:sampn
      ) |> dynfever_summarize()
    )$summary[state == "C"][, vax_pop := "d" ],
    (dynfever_sample(
      n = sampn, parms = dynfever_params(vax_humans = coverage),
      ref.seed = 1:sampn
    ) |> dynfever_summarize()
    )$summary[state == "C"][, vax_pop := "h" ]
  )[, vax_cov := coverage ]
}) |> rbindlist())

dynfever_DAIDD_distros$state <- NULL

dynfever_DAIDD_distros |> setkey(vax_cov, vax_pop, pop, sample) |>
  setcolorder()

useclinic::use_data(dynfever_DAIDD_distros, overwrite = TRUE)
