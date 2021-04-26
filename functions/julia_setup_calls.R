#Julia setup
 Sys.setenv(JULIA_NUM_THREADS = "8")
 
 julia <- JuliaCall::julia_setup("/julia-1.5.4/bin") 


if (exists("check.packages") && check.packages){
  julia_install_package_if_needed("Random")
  julia_install_package_if_needed("Plots")
  julia_install_package_if_needed("Turing")
  julia_install_package_if_needed("Distributions")
  julia_install_package_if_needed("StatsPlots")
  julia_install_package_if_needed("DataFrames")
  julia_install_package_if_needed("MCMCChains")
}

julia_library("Random")
julia_library("Plots")
julia_library("Turing")
julia_library("Distributions")
julia_library("StatsPlots")
julia_library("DataFrames")
julia_library("MCMCChains")