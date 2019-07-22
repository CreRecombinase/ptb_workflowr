  # This file serves the r_*() functions (e.g. r_make()) documented at
 # Define your custom code as a bunch of functions.


source("R/config.R")
source("R/plan.R")      # Create your drake plan.

# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().

make(plan,
     parallelism = data_config$parallelism,
     memory_strategy = data_config$memory_strategy,
     garbage_collection = TRUE,
     jobs = data_config$jobs,
     caching = data_config$caching,
     verbose = 4,
     cache = data_config$cache,prework="future::plan(future::multiprocess)")
