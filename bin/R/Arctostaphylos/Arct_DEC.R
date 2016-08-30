## Cambiando directorio de salida
setwd("~/MEGAsync/bio_comparadaII/result/Arctostaphylos/Biogeo/")

### Calculando ML para diferentes modelos biogeograficos para Arctostaphylos

library(BioGeoBEARS)
library(snow)
library(parallel)
library(optimx)
library(FD)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully

## arbol 
trfn <- "~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/output/Arctostaphylos.new"
tr <-  read.tree("~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/output/Arctostaphylos.new")

#plot(tr)
#title("Example Arctostaphylos")
#axisPhylo()

## datos geograficos
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = "~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/Arctostaphylos_geo.dat")
geogfn <- "~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/Arctostaphylos_geo.dat"
max_range_size <- 3



## Corriendo DEC
BioGeoBEARS_Arct_DEC = define_BioGeoBEARS_run()
# locacion de la filogenia
BioGeoBEARS_Arct_DEC$trfn = trfn
# locacion de los datos geograficos
BioGeoBEARS_Arct_DEC$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_Arct_DEC$max_range_size = max_range_size

# Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_Arct_DEC$min_branchlength = 0.000001
# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_Arct_DEC$include_null_range = TRUE
# shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_Arct_DEC$speedup = TRUE
# if FALSE, use optim() instead of optimx()
BioGeoBEARS_Arct_DEC$use_optimx = TRUE
BioGeoBEARS_Arct_DEC$num_cores_to_use = 1

BioGeoBEARS_Arct_DEC$force_sparse = FALSE  


BioGeoBEARS_Arct_DEC = readfiles_BioGeoBEARS_run(BioGeoBEARS_Arct_DEC)

BioGeoBEARS_Arct_DEC$return_condlikes_table = TRUE
BioGeoBEARS_Arct_DEC$calc_TTL_loglike_from_condlikes_table = TRUE
# get ancestral states from optim run
BioGeoBEARS_Arct_DEC$calc_ancprobs = TRUE



# This contains the model object
#BioGeoBEARS_Arct_DEC$BioGeoBEARS_model_object

# This table contains the parameters of the model 
#BioGeoBEARS_Arct_DEC$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_Arct_DEC)



# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Arctocstaphylos_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_Arct_DEC)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

