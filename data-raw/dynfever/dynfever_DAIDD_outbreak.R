
library(ICI3D)
library(data.table)

#' setup parameters -
params <- list(
  list(), list(), list(), list(),
  list(vax_dog = 0.4), list(vax_dog = 0.5),
  list(vax_human = 0.5, vax_dog = 0.2),
  list(vax_human = 0.8)
) |> lapply(\(p) do.call(dynfever_params, args = p))

#' setup the seeds for dynfever simulations
seed.dt <- data.table(sample = 1L:8L, seed = c(8, 2, 9, 1, 6, 15, 14, 4))

dynfever_DAIDD_outbreak <- dynfever_sample(
  n = length(params), parms = params, ref.seed = seed.dt$seed
)

useclinic::use_data(
  dynfever_DAIDD_outbreak,
  overwrite = TRUE
)
