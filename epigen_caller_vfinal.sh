#!/bin/bash

# replace the path here with your directory containing the scripts
cd YOUR_WORKING_DIRECTORY

#### Simulating different epimutation rates #####

## making strings for model parameters
m="0,0.001,0.005,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5"
r2="0.001,0.1,0.5"
s_e="0.01,0.1,0.5"
mu="0,0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5"
phi="minus1,0,0.5,1" # minus = - but getopt doesn't like the - character

## making function_string string
# the simulations use transition rates when performing epimutation, rather than mu and phi
# so we use this function to calulate the transmission rates from mu and phi
function_string='param_mat[,"i_2"]=param_mat[,"i_1"]*(1-param_mat[,"c_1"]);param_mat[,"c_2"]=param_mat[,"i_1"]*(1+param_mat[,"c_1"]);param_mat[,"i_1"]=param_mat[,"c_2"];param_mat[,"c_1"]=param_mat[,"i_2"]'

## Run script
Rscript Epigen_sim2.6.1.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-14 \
    --log_gen 1000 \
    --pre_contact_mig T \
    --param_split "c_1" \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"

#### Simulating different epimutation rates with pre-selection epimutation ####

## model parameters stay the same

## function_string stays the same

## Run pre-selection epimutation script
Rscript Epigen_sim2.6.1_pre-sel.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-14 \
    --log_gen 1000 \
    --pre_contact_mig T \
    --param_split "c_1" \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"

#### Simulating different adaptive skews #####

## making strings for model parameters
m="0,0.001,0.005,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5"
r2="0.5"
s_e="0.5"
mu="0.5"
phi="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"

## function_string stays the same

## Run script
Rscript Epigen_sim2.6.1.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-14 \
    --log_gen 1000 \
    --pre_contact_mig T \
    --param_split "c_1" \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"

#### Simulationm recording RI over time ####

## making strings for model parameters
m="0.1"
r2="0.001"
s_e="0.5"
mu="0,0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5"
# mu="0.01" # this works as expected
phi="1"
    
## function_string stays the same

## Run script
Rscript Epigen_sim2.6.1.R \
    \
    --max_gen 1001000 \
    --max_burnin 1e6 \
    --equil_threshold 1e-14 \
    --burnin_equil_loci "E" \
    --equil_loci "" \
    --log_gen 1000 \
    --pre_contact_mig T \
    --param_split "i_1" \
    --record_loci "B" \
    --record_gen 1 \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"
    