using Random
using Plots
using Turing
using Distributions
using StatsPlots
using StatsBase
using DataFrames
using MCMCChains
using CSV
using FillArrays
#using TableView

Random.seed!(99);
nchains = 4
niter = 10000
nburnini = 5000

cd("/Users/rod/Dropbox/Coronavirus/EPCOVID_V2/functions/Model/")

include("betabinomial.jl")
include("betabinomialFAST.jl")

#PCR
delta = 0.86
gamma = 0.9999


n = [165, 879, 473]
z = [0, 12, 3]

sim = betabinomial(n, z, 1, 1, delta, gamma)
hmcsample = Turing.sample(sim, HMC(0.01, 5), 100) #starting point



