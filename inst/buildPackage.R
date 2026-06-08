library(devtools)
library(usethis)

usethis::create_package("C:/Users/vandambatesp/Documents/GitHub/acousticMove")

#
# options(
#   usethis.description = list(
#     "Authors@R" = utils::person(
#       "Paul", "van Dam-Bates",
#       email = "paul.vandambates@gmail.com",
#       role = c("aut", "cre")
#     ),
#     License = "MIT + file LICENSE"
#   )
# )
# use_description(fields = list(), check_name = TRUE, roxygen = TRUE)

use_package("RTMB")
use_package("R6")
# use_package("dplyr")
# use_package("ggplot2")

use_build_ignore(c("buildPackage.R"))

usethis::use_tidy_description()
devtools::document()
devtools::check()
