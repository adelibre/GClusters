library(tidyverse)
library(stringr)


# Read header
hdr <- scan(file = "data/GCS_table.txt", skip = 36, nlines = 1, what = "character")

# Manualy set header names
# I tried to automate this, but couldn't find a way
header_name <- c( "Galaxy", "OtherID", "RA", "DEC", "Type" , "D", "D_u",
                  "Method", "A_V", "M_V_T", "M_V_T_u", "M_K", "M_K_u",
                  "N_GC", "N_GC_u", "Source", "sig_e", "sig_e_u", "R_e",
                  "R_e_u", "lg_M_d", "lg_M_d_u", "lg_M_G", "lg_M_G_u",
                  "lg_M_B", "lg_M_B_up", "lg_M_B_low")

# Read data
galaxy <- read.table(file = "data/GCS_table.txt", header = FALSE,
                  skip = 39, nrows = 422, na.strings = "nd",
                  col.names = header_name)

# There is still a problem in the colomn "N_GC" !
# Instead of a numeric, the authors wrote a ">3". So this lead to conversion
# of column into factor instead of numeric
# first relevel to a number
levels(galaxy$N_GC)[levels(galaxy$N_GC)==">3"] <- "4"
# Then convert to a colomn of numeric
galaxy <- galaxy %>% mutate_at(14, parse_number)


save(galaxy, file = "rda/galaxy.rda")

