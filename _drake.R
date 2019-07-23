# This file serves the r_*() functions (e.g. r_make()) documented at
# https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity # nolint
# and
# https://ropensci.github.io/drake/reference/r_make.html

# Load your packages and supporting functions into your session.
# If you use supporting scripts like the ones below,
# you will need to supply them yourself. Examples:
# https://github.com/wlandau/drake-examples/tree/master/main/R
#
#
source("R/config.R")
source("R/plan.R")      # Create your drake plan.
# Define your custom code as a bunch of functions.    # Create your drake plan.

# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().

future::plan(future::multiprocess)
make(plan,
    # parallelism = data_config$parallelism,
     memory_strategy = data_config$memory_strategy,
     garbage_collection = TRUE,
#     jobs = data_config$jobs,
     caching = data_config$caching,
     verbose = 4,
     cache = data_config$cache,prework = quote(future::plan(future::multiprocess)))
