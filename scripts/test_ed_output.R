library(redr)
library(PEcAn.ED2)
library(PEcAnRTM)
import::from(here, here)

img_path <- NULL
edr_exe_path <- "/projectnb/dietzelab/ashiklom/ED2/EDR/build/ed_2.1-opt"

test_dir <- here("test")
out_dir <- here("test-edr")
dir.create(out_dir, showWarnings = FALSE)

list.files(file.path(test_dir, "out"))

ed2in_orig <- read_ed2in(file.path(test_dir, "ED2IN"))

s <- setup_edr(
  ed2in = ed2in_orig,
  output_dir = out_dir,
  datetime = ISOdatetime(1983, 07, 01, 12, 00, 00)
)

p <- prospect(defparam("prospect_5"), 5)

pfts <- pft_lookup$pft_name
npft <- length(pfts)
spectra_list <- rep(list(p), npft)
trait.values <- vector("list", 6)
names(trait.values) <- names(spectra_list) <- pfts

a <- EDR(
  img_path = img_path,
  ed2in_path = s,
  spectra_list = spectra_list,
  trait.values = trait.values,
  edr_exe_path = edr_exe_path
)
