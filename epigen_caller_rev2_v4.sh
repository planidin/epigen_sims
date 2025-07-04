#!/bin/bash

# replace the path here with your directory containing the scripts
cd /Users/anon/Desktop/MS5/drafts/proc_B_rev2/rev2_code

#### Simulating different epimutation rates #####

## making strings for model parameters
m="0.000001,0.001,0.005,0.01,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5"
r2="0,0.001,0.1,0.5"
s_e="0.01,0.05,0.1,0.5"
mu="0,0.02,0.1,0.2,0.3,0.4,0.5"
phi="minus1,0,0.5,1"

## making function_string string
# the simulations use transition rates when performing epimutation, rather than mu and phi
# so we use this function to calulate the transmission rates from mu and phi
function_string='param_mat[,"i_2"]=param_mat[,"i_1"]*(1-param_mat[,"c_1"]);param_mat[,"c_2"]=param_mat[,"i_1"]*(1+param_mat[,"c_1"]);param_mat[,"i_1"]=param_mat[,"c_2"];param_mat[,"c_1"]=param_mat[,"i_2"]'

## Run script
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
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

#### Simulating different epimutation rates with secondary contact ####

## Run script with no pre-contact migration (--pre_contact_mig F)
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
    --log_gen 1000 \
    --pre_contact_mig F \
    --param_split "c_1" \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"

#### Simulating different epimutation rates with pre-selection epimutation ####

## Run script with pre-selection epimutation (--pre_sel_epi T)
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
    --log_gen 1000 \
    --pre_contact_mig T \
    --pre_sel_epi T \
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
m="0.000001,0.001,0.005,0.01,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5"
r2="0.001,0.1,0.5"
s_e="0,0.01,0.1,0.5"
mu="0.02,0.1,0.5"
phi="0,0.25,0.5,0.75,1"

## function_string stays the same

## Run script
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
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

#### Simulating different strength of selection #####
# only used for epimutation-selection balance plots

## making strings for model parameters
m="0"
r2="0.5"
s_e="0.000001,0.001,0.0025,0.005,0.0075,0.01,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5"
mu="0,0.02,0.1,0.2,0.3,0.4,0.5"
phi="0,0.25,0.5,0.75,1"

## function_string stays the same

## Run script
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
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

#### Simulating different epimutation rates over recombination rates ####

## making strings for model parameters
m="0.000001,0.001,0.1,0.5"
r2="0,0.000001,0.001,0.0025,0.005,0.0075,0.01,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5"
s_e="0.01,0.1,0.5"
mu="0,0.02,0.1,0.2,0.3,0.4,0.5"
phi="0,0.5,0.9,1"

## Run script
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
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


#### Simulating different epimutation rates and recording LD ####

m="0.001,0.1,0.5"
r2="0,0.000001,0.001,0.0025,0.005,0.0075,0.01,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5"
s_e="0.01,0.1,0.5"
mu="0,0.02,0.1,0.2,0.3,0.4,0.5"
phi="0,1"

## Run script recording LD (--record_gen 1) & (--record_LD "Dnorm_eB_3")
Rscript Epigen_sim2.6.3.R \
    \
    --max_gen 4e6 \
    --max_burnin 2e6 \
    --equil_threshold 1e-12 \
    --log_gen 1000 \
    --pre_contact_mig T \
    --record_gen 1 \
    --record_LD "Dnorm_eB_3" \
    --param_split "c_1" \
    --m $m \
    --r2 $r2 \
    --s.e $s_e \
    --i_1 $mu \
    --c_1 $phi \
    \
    --function_string "$function_string"

#### Simulating different periods of allopatry ####

m="0.01,0.1,0.5"
r2="0,0.001,0.1,0.5"
s_e="0.01,0.1,0.5"
mu="0,0.02,0.1,0.2,0.3,0.4,0.5"
phi="1"

for i in $(seq 0 1 20) $(seq 25 5 95) $(seq 100 100 1000)
do
    echo $i
    Rscript Epigen_sim2.6.3.R \
        \
        --max_gen 4e6 \
        --max_burnin 2e6 \
        --allo_gen $i \
        --equil_threshold 1e-12 \
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
done
