#!/usr/bin/Rscript

###############################################################
## v2.6.1 - Implement LD calculation                         ##
## removed induct_start and converted e_X to c_X             ##
## merged secondary_contactC and injection_migrationC.       ##
## TO DO: implement LD calculation in Cpp                    ##
## TO DO: set up checking for equilibrium each X generations ##
## TO DO: merge record loci and record LD, make record loci  ##
## by deme possible                                          ##
## TO DO: properly implement non-random B introduction.      ##
## TO DO: implement population size                          ##
## TO DO: output and log directory                           ##
###############################################################

## set up a clean work space
rm(list=ls())
sim_num = "sim2.6.1"

if(interactive()){
    # change to your working directory
    setwd("YOUR_WORKING_DIRECTORY")
}

## begin writing to log file and create matching .csv output
run_date_time = format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
out_text_file = paste0("RI_",sim_num,"_",run_date_time,".log")
out_csv_file = paste0("RI_",sim_num,"_",run_date_time,".csv")
cat("\nLogging model progress at: ",getwd(),"/",out_text_file,"\n\n", sep = "")
sink(out_text_file)

#### Load Packages and functions ####

cat("\n## Loading packages ##\n\n")

# loading doParallel package, for multithreading simulations
if(!require(doParallel)) {
    stop("'doParallel' package is required")
} else {
    cat("'doParallel' loaded\n")
}
# loading Rcpp package, to run faster C++ versions of simulation functions
if(!require(Rcpp)) {
    stop("'Rcpp' package is required")
} else {
    cat("'Rcpp' loaded\n")
}
# loading EpiSims package, which contains the custom C++ functions
# NOTE: it is necessary to save the custom C++ functions as part of a package
# so that the custom functions can be quickly loaded in parallel for multithreading
if(!require(EpiSims)) {
    stop("'EpiSims' package is required")
} else {
    cat("'EpiSims' loaded\n")
}
# loading getopt package, for passing model parameters to the Rscript from bash
if(!require(getopt)) {
    stop("'getopt' package is required")
} else {
    cat("'getopt' loaded\n")
}

# loading custom R functions from another file
# NOTE: functions were put in another .R file for organizational purposes
cat("sourcing functions from 'epigen_funcs_R1.3.R'\n")
source("epigen_funcs_R1.3.R")

### If the EpiSims package is not installed 
### it can be compiled from source and loaded with these

# Rcpp.package.skeleton("EpiSims", example_code = F, cpp_files = "epigen_funcs1.4.cpp")
# install.packages("./EpiSims",repos=NULL, type = "source")

### or sourced with this, though sourcing Rcpp functions does not work with doParallel

# sourceCpp('epigen_funcs1.4.cpp')

#### Setting simulation parameters ####

# specifying a getopt object to get simulation parameters from bash
spec = matrix(c(
    'help',              'h',    0, "logical",
    'max_gen',           'mg',   1, "character",
    'max_burnin',        'mb',   1, "character",
    'locus_order',       'lo',   1, "character",
    'equil_loci',        'el',   1, "character",
    'burnin_equil_loci', 'bl',   1, "character",
    'equil_threshold',   'et',   1, "character",
    'record_loci',       'rl',   1, "character",
    'record_LD',         'rd',   1, "character",
    'record_gen',        'rg',   1, "character",
    'post_contact_mig',  'am',   1, "logical",
    'pre_contact_mig',   'bm',   1, "logical",
    'inject_delay',      'id',   1, "character",
    'log_gen',           'lg',   1, "character",
    'param_split',       'ps',   1, "character",
    'm',                 'm',    1, "character",
    'r1',                'r1',   1, "character",
    'r2',                'r2',   1, "character",
    's.a',               'sa',   1, "character",
    's.e',               'se',   1, "character",
    'h.a',               'ha',   1, "character",
    'h.e',               'he',   1, "character",
    'i_1',               'i1',   1, "character",
    'i_2',               'i2',   1, "character",
    'c_1',               'c1',   1, "character",
    'c_2',               'c2',   1, "character",
    'i_gen_1',           'ig1',  1, "character",
    'i_gen_2',           'ig2',  1, "character",
    'c_gen_1',           'cg1',  1, "character",
    'c_gen_2',           'cg2',  1, "character",
    'h.i_gen_1',         'hi1',  1, "character",
    'h.i_gen_2',         'hi2',  1, "character",
    'h.c_gen_1',         'hc1',  1, "character",
    'h.c_gen_2',         'hc2',  1, "character",
    'epistat',           'epi',  1, "character",
    'function_string',   'fs',   1, "character"
),byrow=T,ncol = 4)

opt = getopt(spec)

help_message = "

### META-PARAMETERS ###

    --max_gen | mg (Max number of generations)
        Parameter sets that do not reach post-contact equilibrium before this many generations will be output
        and have RI calculated as-is
        Default: 2e5  

    --max_burnin | mb (Max number of burn-in generations)
        The number of generations until contact of neutral B allele. Parameter sets that do not reach equilibrium 
        before this many generations point will undergo contact of the neutral locus. After pre-contact 
        equilibrium is reached, the generation is set to max_burnin such that the maximum number of post-contact
        generations is max_gen - max_burnin.
        Default: 1e5
    
    --locus_order | lo (Locus order)
        Determines the order of the three loci for recombination. The format of this option must be the three locus
        names seperated by commas e.g., 'A,B,E'. For example, if the locus order is A,B,E the neutral locus is between
        the genetic locus under selection and the epigenetic locus. In this case r1 determines the rate of
        recombination between A and B, and r2 determines the rate of recombination between B and E.
        Default: 'A,B,E'
    
    --equil_loci | el (Loci measured for equilibrium)
        Determine which loci are measured for reaching post-contact equilibrium. The format of this option must be
        the locus names seperated by commas e.g., 'B,E'. When all loci listed reach
        equilibrium, simulation of that parameter set is stopped and RI is measured
        Default: 'B,E'
    
    --burnin_equil_loci | bl (Loci measured for burn-in equilibirum)
        Determine which loci are measured for reaching pre-contact equilibrium. The format of this option must be
        the locus names seperated by commas e.g., 'A,E'. When all loci listed reach
        equilibrium, contact of the neutral B allele occurs.
        Default: 'E' 
    
    --equil_threshold | et (Equilibrium threshold)
        Threshold of difference in allele frequency between generations to halt a parameter set.
        Set to 0 to run all parameter sets to the maximum number of generations
        Default: 1e-8
    
    --record_loci | rl (Loci with frequencies recorded each generation)
        Determine which loci have their frequencies recorded during simulations. The format of this option 
        must be the locus names seperated by commas e.g., 'B,E'.
        Default: ''
    
    --record_LD | rd (Types of linkage disequilibrium recorded each generation)
        Determine which loci pairs have their linkage disequilibrium (LD) calulated and recorded
        during simulations. A single instance of recording LD is specified as X_Y_Z. where X is the LD metric,
        Y is the pair of loci and Z is the deme. Possible LD metrics are: 'D' (basic LD),
        'Dnorm' (normalized LD i.e. D'), and 'cov' (covariation). Possible locus pairs are: 'AE', 'AB', and 'EB'
        (order of the loci does not matter). Possible demes are '1', '2', and '0', where zero averages
        over both demes. For example, calculating normalized LD between locus E and B in deme 1 is 'Dnorm_EB_1'.
        To calculate multiple types of LD, seperate multiple instances by commans e.g., 'Dnorm_EB_1,Dnorm_EB_2'.
        Default: ''
      
    --record_gen | rg (Record allele frequency each X generations)
        Set how frequently allele frequencies and LD are recorded. Set to zero to not record. After a parameter
        set is halted recordings are filled with NA. WARNING: makes potentially HUGE output files, so I have
        hard coded in some limits e.g., after 10000 generations, recording cannot be more frequent than once
        every 1000 generations (see lines XXX-XXX).
        Default: 0
    
    --post_contact_mig | am (Enables migration after contact of the neutral B allele)
        If true, there is migration between demes after contact of the neutral B allele.
        Set to false and a small m to model analytical equations from [Anonymous 2022]
        Default: T
    
    --pre_contact_mig | bm (Enables migration before contact of the neutral B allele)
        If true, there is migration between demes before contact of the neutral B allele.
        Set to F to model secondary contact.
        Default: T
    
    --inject_delay | id (Delay contact of the neutral allele X generations)
        Determine the number of generations after contact (after burn-in generations are complete)
        to introduce the neutral B allele. Use in conjunction with pre_contact_mig=F to simulate 
        contact at different stages between secondary and primary contact.
        Default: 0
      
    --log_gen | lg (Write model progress to log file each X generations)
        Writes the number of parameter sets and generation of each worker to the log file every X generations.
        Set to zero to not write model progress to log file
        Default: 1000
    
    --param_split | ps (Divide parameter list by unique combinations of these parameters)
        Based on this option, the parameter list is divided into a list-of-lists, each containing a unique
        value of the parameter(s) specified. Each list within the list-of-lists will be used to calculate
        all possible parameter combinations, which are then simulated in parallel. Functionally, the number
        of unique combinations of the specified parameters determines the maximum number of workers that can
        be run in parallel. Multiple parameters can be passed to this option if they are seperated by commas
        e.g., --param_split='m,s.e' will split by unique combinations of values for migration and selection at
        the epigenetic locus. For example, say you run a simulation with m='0.001,0.1,0.5', s.e='0.1,0.5' and all
        other parameters have only one value. Then by using --param_split='m', you simulate m=0.001, m=0.1 and m=0.5
        each on their own threads, where each thread runs simulations for s.e=0.1 and s.e=0.5. If you use 
        --param_split='m,s.e', each combination of values of m and s.e will be run seperately, in this case
        allowing for up to 6 workers at a time. If the number of lists is greater than the maximum number of
        threads, new workers will be launched as previous workers finish.
        Default: 'm'
    
### SIMULATION PARAMETERS ###
  
    All model parameters are given as a list of values seperated by commas e.g. '0.1,0.2,0.3'.
    From the model parameters given, all unique combinations of model parameters are generated and simulated.
    For example, if you give m='0.01,0.2' and s.e='0.1,0.5', four simulations will be run i.e.,
    m=0.01 & s.e=0.1, m=0.2 & s.e=0.1, m=0.01 & s.e=0.5, and m=0.2 & s.e=0.5
    
    --m | m (Migration rate)
        The proportion of individuals which migrate between demes each generation.
        Default: '0.01'
    
    --r1 | r1 (Recombination rate between the 1st and 2nd locus)
        The proportion of the population where alleles at the 1st locus swap strands.
        Non-interfering recombination is assumed, such that the 2nd locus recombines
        away from the 1st and 3rd loci at rate r1*r2
        Default: '0.5'
    
    --r2 | r2 (Recinbination rate between the 2nd and 3rd locus)
        The proportion of the population where alleles at the 3rd locus swap strands.
        Non-interfering recombination is assumed, such that the 2nd locus recombines
        away from the 1st and 3rd loci at rate r1*r2
        Default: '0.5'
    
    --s.a | sa (Selection at the genetic locus)
        The strength of divergent selection on the genetic locus (A). The A allele is favored
        in deme one, with additive fitnesses of AA=1, Aa=1-s.a*h.a and aa=1-s.a. The a allele is
        favored in deme two, with with additive fitnesses of aa=1, Aa=1-s.a*h.a and AA=1-s.a.
        These values assume no selection or epistasis from the epigenetic locus.
        Default: '0'
    
    --s.e | se (Selection at the epigenetic locus)
        The strength of divergent selection on the epigenetic locus (E). The E allele is favored
        in deme one, with additive fitnesses of EE=1, Ee=1-s.e*h.e and ee=1-s.e. The e allele is
        favored in deme two, with with additive fitnesses of ee=1, Ee=1-s.e*h.e and EE=1-s.e.
        These values assume no selection or epistasis from the epigenetic locus.
        Default: '0.5'
    
    --h.a | ha (Dominance of genetic A allele)
        See s.a for how h.a affects fitness.
        Default: '0.5'
      
    --h.e | he (Dominance of epigenetic E allele)
        See s.e for how h.e affects fitness.
        Default: '0.5'
    
    --i_1 | i1 (Rate of e to E transitions in deme 1)
        NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --i_2 | i2 (Rate of e to E transitions in deme 2)
        NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --c_1 | c1 (Rate of E to e transitions in deme 1)
        NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --c_2 | c2 (Rate of E to e transitions in deme 2)
        NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --i_gen_1 | ig1 (Genetic influence on the rate of e to E transitions in deme 1)
        Percent of induction determined by genetic allele in habitat 1. NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --i_gen_2 | ig2 (Genetic influence on the rate of e to E transitions in deme 2)
        Percent of induction determined by genetic allele in habitat 1. NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --c_gen_1 | cg1 (Genetic influence on the rate of E to e transitions in deme 1)
        Percent of induction determined by genetic allele in habitat 1. NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --c_gen_2 | cg2 (Genetic influence on the rate of E to e transitions in deme 2)
        Percent of induction determined by genetic allele in habitat 1. NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
    --h.i_gen_1 | hi1 (Dominance of genetic influence on e to E transitions in deme 1)
        NEED TO ADD MORE DESCRIPTION
        Default: '0.5'
    
    --h.i_gen_2 | hi2 (Dominance of genetic influence on e to E transitions in deme 2)
        NEED TO ADD MORE DESCRIPTION
        Default: '0.5'
    
    --h.c_gen_1 | hc1 (Dominance of genetic influence on E to e transitions in deme 1)
        NEED TO ADD MORE DESCRIPTION
        Default: '0.5'
    
    --h.c_gen_2 | hc2 (Dominance of genetic influence on E to e transitions in deme 2)
        NEED TO ADD MORE DESCRIPTION
        Default: '0.5'

    --epistat | epi (Epistasis coefficient)
        Coefficient for the interactive effects of genetic and epigenetic alleles in determining
        fitness. NEED TO ADD MORE DESCRIPTION
        Default: '0'
    
### Parameter manipulation function ###

    If the default behaviour of simulating all unique combinations of parameters is not wanted
    this function can be used to subset parameter combinations or to calculate parameters from
    other values passed into the parameter matrix.

    --function_string | fs (Function to change parameter matrix on each worker)
        The string written here is exectued on each worker individually, after the parameter
        matrix is calculated from the parameter list, but before any derived parameters such as
        fitness are calculated
        Default: 'paste(\"no function applied\")'

### Saved Parameter manipulation functions ###

    Some useful parameter manipulation functions. Note that there the function strings are
    divided onto multiple lines for readability but when they are written in bash they
    need to be on a single line. In bash, the function_string must be incased in single quotes ('')
    with double quotes for quotes within R. The double quotes should not be backslash escaped in bash.

    No modifications
        --function_string='paste(\"no function applied\")'
               
    Neutral epimutation
        --function_string='param_mat = subset(param_mat, param_mat[,\"i_1\"] == param_mat[,\"c_1\"] &
                           param_mat[,\"i_2\"] == param_mat[,\"c_2\"] &
                           param_mat[,\"c_1\"] == param_mat[,\"c_2\"])'
    Pass values of epimutation rate to i_1, then set all other i/c parameters to 0.
    This function then sets all i/c parameters to be the same as i_1, creating only parameter
    combinations of neutral epimutation.
  
    Adaptively-biased epimutation
        --function_string = 'param_mat[,\"c_2\"] = param_mat[,\"i_1\"]'
    Pass values of epimutation rate to i_1, then set all other i/c parameters to 0.
    This function then sets all c_2 to be the same as i_1, creating only parameter
    combinations of adaptive epimutation.  
  
    Reparameterization in terms of epimutation rate (mu) and adaptive induction (gamma)
        --function_string='param_mat[,\"i_2\"]=param_mat[,\"i_1\"]*(1-param_mat[,\"c_1\"]);
                           param_mat[,\"c_2\"]=param_mat[,\"i_1\"]*(1+param_mat[,\"c_1\"]);
                           param_mat[,\"i_1\"]=param_mat[,\"c_2\"];
                           param_mat[,\"c_1\"]=param_mat[,\"i_2\"]'
    Pass values of epimutation rate to i_1 and pass values of adaptive induction to c_1,
    then set all other i/c parameters to 0. This function calculates i/c values in deme 2
    using i_1 and c_1 as epimutation rate and adaptive indunction respectively. Then i_1 and c_1
    are overwritten to mirror values in deme 2.
                
"

# print the help message if it is asked for
if(!is.null(opt$help)) {
    sink()
    cat(getopt(spec,usage=T))
    cat(help_message)
    q(status=1)
}

# if running simulations as an interactive R session, set parameter values here
if(interactive()) {
  
    opt$max_gen           = 1000500    
    opt$max_burnin        = 1e6    
    opt$locus_order       = "A,B,E"
    opt$equil_loci        = ""  
    opt$burnin_equil_loci = "E"    
    opt$equil_threshold   = 1e-6   
    opt$record_loci       = "B"
    opt$record_LD         = ""
    opt$record_gen        = 1     
    opt$post_contact_mig  = T   
    opt$pre_contact_mig   = T
    opt$inject_delay      = 0
    opt$log_gen           = 1000   
    opt$param_split       = "i_1"    
    opt$m                 = "0.1" 
    opt$r1                = "0.5"  
    opt$r2                = "0.001"  
    opt$s.a               = "0"    
    opt$s.e               = "0.5"  
    opt$h.a               = "0.5"  
    opt$h.e               = "0.5"  
    # opt$i_1               = "0,0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5"
    opt$i_1               = "0.01"
    opt$i_2               = "0"    
    opt$c_1               = "1"    
    opt$c_2               = "0"    
    opt$i_gen_1           = "0"    
    opt$i_gen_2           = "0"    
    opt$c_gen_1           = "0"    
    opt$c_gen_2           = "0"    
    opt$h.i_gen_1         = "0.5"  
    opt$h.i_gen_2         = "0.5"  
    opt$h.c_gen_1         = "0.5"  
    opt$h.c_gen_2         = "0.5"  
    opt$epistat           = "0"    
    opt$function_string   = "param_mat[,\"i_2\"]=param_mat[,\"i_1\"]*(1-param_mat[,\"c_1\"]);param_mat[,\"c_2\"]=param_mat[,\"i_1\"]*(1+param_mat[,\"c_1\"]);param_mat[,\"i_1\"]=param_mat[,\"c_2\"];param_mat[,\"c_1\"]=param_mat[,\"i_2\"]"
    # opt$function_string   = "param_mat[,c(\"c_1\",\"c_2\",\"i_2\")]=param_mat[,\"i_1\"]"
    # opt$function_string   = "param_mat[\"_2\"]=param_mat[,\"i_1\"]"
    # opt$function_string = ""
  
# otherwise parameter values which are not given in bash are
# filled in with the default values
} else { 
  
    if(is.null(opt$max_gen)          ) {opt$max_gen           = 2e5    }
    if(is.null(opt$max_burnin)       ) {opt$max_burnin        = 1e5    }
    if(is.null(opt$locus_order)      ) {opt$locus_order       = "A,B,E"}
    if(is.null(opt$equil_loci)       ) {opt$equil_loci        = "B,E"  }
    if(is.null(opt$burnin_equil_loci)) {opt$burnin_equil_loci = "E"    }
    if(is.null(opt$equil_threshold)  ) {opt$equil_threshold   = 1e-8   }
    if(is.null(opt$record_loci)      ) {opt$record_loci       = ""     }
    if(is.null(opt$record_LD)        ) {opt$record_LD         = ""     }
    if(is.null(opt$record_gen)       ) {opt$record_gen        = 0      }
    if(is.null(opt$induct_start)     ) {opt$induct_start      = F      }
    if(is.null(opt$post_contact_mig) ) {opt$post_contact_mig  = T      }
    if(is.null(opt$pre_contact_mig)  ) {opt$pre_contact_mig   = T      }
    if(is.null(opt$inject_delay)     ) {opt$inject_delay      = 0      }
    if(is.null(opt$log_gen)          ) {opt$log_gen           = 1000   }
    if(is.null(opt$param_split)      ) {opt$param_split       = "m"    }
    if(is.null(opt$m)                ) {opt$m                 = "0.01" }
    if(is.null(opt$r1)               ) {opt$r1                = "0.5"  }
    if(is.null(opt$r2)               ) {opt$r2                = "0.5"  }
    if(is.null(opt$s.a)              ) {opt$s.a               = "0"    }
    if(is.null(opt$s.e)              ) {opt$s.e               = "0.5"  }
    if(is.null(opt$h.a)              ) {opt$h.a               = "0.5"  }
    if(is.null(opt$h.e)              ) {opt$h.e               = "0.5"  }
    if(is.null(opt$i_1)              ) {opt$i_1               = "0"    }
    if(is.null(opt$i_2)              ) {opt$i_2               = "0"    }
    if(is.null(opt$c_1)              ) {opt$c_1               = "0"    }
    if(is.null(opt$c_2)              ) {opt$c_2               = "0"    }
    if(is.null(opt$i_gen_1)          ) {opt$i_gen_1           = "0"    }
    if(is.null(opt$i_gen_2)          ) {opt$i_gen_2           = "0"    }
    if(is.null(opt$c_gen_1)          ) {opt$c_gen_1           = "0"    }
    if(is.null(opt$c_gen_2)          ) {opt$c_gen_2           = "0"    }
    if(is.null(opt$h.i_gen_1)        ) {opt$h.i_gen_1         = "0.5"  }
    if(is.null(opt$h.i_gen_2)        ) {opt$h.i_gen_2         = "0.5"  }
    if(is.null(opt$h.c_gen_1)        ) {opt$h.c_gen_1         = "0.5"  }
    if(is.null(opt$h.c_gen_2)        ) {opt$h.c_gen_2         = "0.5"  }
    if(is.null(opt$epistat)          ) {opt$epistat           = "0"    }
    if(is.null(opt$function_string)  ) {opt$function_string   = "paste(\"no function applied\")"}
  
}

#### Convert args meta-parameters to R objects ####

cat("\n## Setting meta-parameters ##\n\n")

# this set of brackets is just so all of these lines run together 
# during an interactive session
{
  
    # as.numeric is used here and below to coerce values written in scientific
    # notation in bash into numeric values in R
    max_gen = as.numeric(opt$max_gen)
    cat("Setting maximum generations to ",max_gen,"\n")
    
    burnin = as.numeric(opt$max_burnin)
    cat("Setting maximum burnin generations to ",burnin,"\n")
    
    # strsplit and unlist and used here and below to convert strings of comma
    # separated values in bash into vectors in R
    locus_order = unlist(strsplit(opt$locus_order,","))
    # checks that the string provided for locus_order creates a valid locus order
    if(sum(locus_order %in% "A") == 1  &
        sum(locus_order %in% "B") == 1  &
        sum(locus_order %in% "E") == 1) {
        cat("Setting locus order to",locus_order,"\n")
    } else {
        stop("Invalid locus order given, please provide\n
              a character vector containing one of each: A B E")
    }
    
    equilibrium_loci = unlist(strsplit(opt$equil_loci,","))
    cat("Setting equilibrium_loci to",equilibrium_loci,"\n")
    
    burnin_equilibrium_loci = unlist(strsplit(opt$burnin_equil_loci,","))
    cat("Setting burnin_equilibrium_loci to",burnin_equilibrium_loci,"\n")
    
    equilibrium_threshold = as.numeric(opt$equil_threshold)
    cat("Setting equilibrium_threshold to",equilibrium_threshold,"\n")
    
    record_loci = unlist(strsplit(opt$record_loci,","))
    record_LD = unlist(strsplit(opt$record_LD,","))
    record_gen = as.numeric(opt$record_gen)
    cat("Recording allele frequency of ",record_loci," each ",record_gen,"generations\n")
    cat("Recording linkage disequilibrium of ",record_LD," each ",record_gen,"generations\n")
    
    post_contact_mig = opt$post_contact_mig
    cat("Setting post-contact migration to",post_contact_mig,"\n")
    
    pre_contact_mig = opt$pre_contact_mig
    cat("Setting pre-contact migration to",pre_contact_mig,"\n")
    
    inject_delay = as.numeric(opt$inject_delay)
    cat("Setting neutral allele injection delay to",inject_delay,"\n")
    
    progress_gen = as.numeric(opt$log_gen)
    cat("Outputting model progress to log file each ",progress_gen,"generations \n")
    
    split_by = unlist(strsplit(opt$param_split,","))
    cat("Splitting parameter matricies by ",split_by," \n")
}

#### Convert args simulation parameters to an R list ####

cat("\n## Setting simulation parameters ##\n\n")

param_list = list(
    # as with the meta-parameters, strsplit, unlist and as.numeric are
    # used to convert strings of comma separated values into vectors in R
    # and coerce scientific notation into numeric values in R
    m =  as.numeric(unlist(strsplit(opt$m,","))),
    r1 = as.numeric(unlist(strsplit(opt$r1,","))),
    r2 = as.numeric(unlist(strsplit(opt$r2,","))),
    s.a = as.numeric(unlist(strsplit(opt$s.a,","))),
    s.e = as.numeric(unlist(strsplit(opt$s.e,","))),
    h.a = as.numeric(unlist(strsplit(opt$h.a,","))),
    h.e = as.numeric(unlist(strsplit(opt$h.e,","))),
    # gsub swaps out minus for - here, as a hack to be able to pass
    # negative values with getopt, since the - symbol is not parsed correctly
    i_1 = as.numeric(unlist(strsplit(gsub("minus","-",opt$i_1),","))),
    i_2 = as.numeric(unlist(strsplit(gsub("minus","-",opt$i_2),","))),
    c_1 = as.numeric(unlist(strsplit(gsub("minus","-",opt$c_1),","))),
    c_2 = as.numeric(unlist(strsplit(gsub("minus","-",opt$c_2),","))),
    i_gen_1 = as.numeric(unlist(strsplit(opt$i_gen_1,","))),
    i_gen_2 = as.numeric(unlist(strsplit(opt$i_gen_2,","))),
    c_gen_1 = as.numeric(unlist(strsplit(opt$c_gen_1,","))),
    c_gen_2 = as.numeric(unlist(strsplit(opt$c_gen_2,","))),
    h.i_gen_1 = as.numeric(unlist(strsplit(opt$h.i_gen_1,","))),
    h.i_gen_2 = as.numeric(unlist(strsplit(opt$h.i_gen_2,","))),
    h.c_gen_1 = as.numeric(unlist(strsplit(opt$h.c_gen_1,","))),
    h.c_gen_2 = as.numeric(unlist(strsplit(opt$h.c_gen_2,","))),
    epistat = as.numeric(unlist(strsplit(opt$epistat,",")))
)

## outputs model parameters to the log file
for(x in 1:length(param_list)) {
    cat(names(param_list)[x],"=",paste(param_list[[x]], collapse=', ' ),"\n")
}

#### Converting custom function to R string ####

function_string = opt$function_string
cat("\nApplying custom parameter function:\n\n")
cat(function_string,"\n")

#### Breaking parameter list into a list-of-lists for parallelization ####

# get all of the possible parameter value combinations
# for parameters in split_by
split_rows = names(param_list) %in% split_by
split_combos = expand.grid(param_list[split_rows])

# create a list of parameter lists, where each item in the list contains
# only one unique value from split_combos
param_list_list = list()
for(i in 1:dim(split_combos)[1]) {
    
    # add all parameters that are not in split_by to the sub-list
    param_list_list[[i]] = param_list[!(names(param_list) %in% split_by)]
    
    # add the unique parameter values of split_combo to the sub-list
    for(j in 1:dim(split_combos)[2]) {
        param_list_list[[i]][split_by[j]] = split_combos[i,j]
    }
    
}

#### Setting up Parallelization ####

cat("\n## Setting up parallelization ##\n\n")

## get number of cores
# -1 to save a core for other processes
n.cores = parallel::detectCores() - 1

## set number of workers to the smaller of
# n.cores or the the length of the parameter list-of-lists
if(length(param_list_list) < n.cores) {
    n.workers = length(param_list_list)
} else {
    n.workers = n.cores
}

## create a parallelization cluster
my.cluster <- parallel::makeCluster(
    n.workers, 
    type = "PSOCK",
    outfile = out_text_file
)

## connect the cluster to doParallel and check that it worked
doParallel::registerDoParallel(cl = my.cluster)
if(!foreach::getDoParRegistered()) {
    stop("Paralell clusters not registered")
}

cat("Parallelizing over ", foreach::getDoParWorkers(), " clusters\n")

#### Run Simulation ####

cat("\n## Running simulation loop ##\n\n")

# loop through the parameter list-of-lists in parallel
# .packages tells doParallel to load the Rcpp and EpiSims on each worker
# so that they can use the functions written in C++
out_mat_list = foreach(i = 1:length(param_list_list),
                       .packages = c("Rcpp","EpiSims")) %dopar% {
    
    #### Set up inputs ####
    
    in_list = param_list_list[[i]] # list of parameter values to be simulated
    allele_combos = make_allele_combos() # generate vector of possible states
    cur_gen = 1 # set generation to 1
    
    ## outputs model parameters to log file
    log_message = paste("\nInitializing worker",i,"with parameters:\n\n")
    for(x in 1:length(param_list)) {
        log_message = paste(log_message,
                            names(in_list)[x],
                            "=",
                            paste(in_list[[x]],collapse=', ' ),
                            "\n")
    }
    cat(log_message,"\n")
    
    #### Build starting objects ####
    
    # create a matrix of all combinations of 
    # parameters in in_list
    param_mat = expand.grid(in_list)
    
    # evaluate parameter matrix manipulation function
    eval(parse(text = function_string))
    
    ## add calculated columns to the parameter matrix
    # calculate fitness of different (epi)genotypes 
    # based on s.e, s.a, h.e, h.a and epistat
    param_mat = add_fitness(param_mat,allele_combos)
    # calculate rates of transitions between different epigenetic states
    # based on the parameters listed in the function call
    param_mat = add_epiparams(param_mat,
                              param_mat[,"i_1"],
                              param_mat[,"i_2"],
                              param_mat[,"c_1"],
                              param_mat[,"c_2"],
                              param_mat[,"i_gen_1"],
                              param_mat[,"i_gen_2"],
                              param_mat[,"c_gen_1"],
                              param_mat[,"c_gen_2"],
                              param_mat[,"h.i_gen_1"],
                              param_mat[,"h.i_gen_2"],
                              param_mat[,"h.c_gen_1"],
                              param_mat[,"h.c_gen_2"])
    
    # add_fitness() and add_epiparams() convert param_mat 
    # to a data.frame, so we convert it back to a matrix 
    # so that it can be passed to C++
    param_mat = as.matrix(param_mat)
    
    ## Create the allele frequency matrix
    # TO DO: simplify allele_start()
    # creates a matrix of (epi)allele frequencies
    # deme 1 starts fixed for EE_AA_bb
    # deme 2 starts fixed for ee_aa_bb
    freq_mat = allele_start(n_param_rows = nrow(param_mat),
                            allele_combos)
    
    ## Create a list of column indices used by model functions
    # allows the passing of certain columns of freq_mat and param_mat 
    # to C++ to save time, rather than passing the entirely of both
    # matrices each time a function is evaluated
    index_list = make_index_list(param_mat,
                                 freq_mat,
                                 locus_order = locus_order)
    
    ## Create matrices for model functions
    # using the pre-generated indices, create matrices of just the
    # required param_mat columns for selection and epimutation
    selection_mat = param_mat[,index_list[["freq2fit"]],drop = F]
    c_mat = param_mat[,index_list[["freq2c"]],drop = F]
    i_mat = param_mat[,index_list[["freq2i"]],drop = F]
    
    ## Create a matrix for outputting parameter sets to during dynamic halting
    # the matrix is freq_mat filled with zeros plus a column for storing the 
    # generation that equilibrium was reached
    out_mat = matrix(0,
                     nrow = nrow(freq_mat),
                     ncol = ncol(param_mat)+ncol(freq_mat)+1)
    colnames(out_mat) = c(colnames(param_mat),
                          colnames(freq_mat),
                          "equil_gen")
    out_mat[,"equil_gen"] = NA
    
    ## Create objects for calculating equilibrium
    # equilibrium_mat stores the frequencies of each (epi)allele and whether they have changed
    # between generations, for each parameter combination
    equilibrium_mat = matrix(0, nrow = nrow(freq_mat), ncol = 9) # 3 columns for each locus, A B E 
    colnames(equilibrium_mat) = c("prop_B_1","prop_B_2","equil_B",
                                  "prop_A_1","prop_A_2","equil_A",
                                  "prop_E_1","prop_E_2","equil_E")
    # cumulative_equilibrium_mat will store the previous generations (epi)allele frequencies
    # but for now it is set to be equal to equilibrium_mat so that the two matricies
    # have the same dimension
    cumulative_equilibrium_mat = equilibrium_mat
    # create a vector to store which rows to output when equilibrium is reached
    output_rows = rep(0,nrow(freq_mat))
    # create a vector to store which rows have been output
    cumulative_output_rows = logical(length = nrow(freq_mat))
    # create a vector to store how many generations to reach 
    # equilibrium during burn-in generations
    burnin_equil_gen = rep(NA,nrow(out_mat))
    
    ## Create a matrix for recording allele frequencies each generation
    # if record_gen == 0, this matrix is not created
    if(record_gen != 0) {
        record_freqs_mat = matrix(nrow = nrow(freq_mat),ncol = 0)
        thin_record_gen = record_gen # create a variable to store more sparse recording in later generations
    }
    
    #### Loop through generations
    
    # loop until max_gen number of generations
    while (cur_gen < max_gen) {
        
        # output the simulation to the log file each
        # progress_gen generations
        if ((cur_gen %% progress_gen) == 0) {
            cat("Worker", i, "with", nrow(freq_mat),"param sets is at gen =",cur_gen,"\n")
        }
        
        # when (epi)allele frequencies reach equilibrium for a certain set 
        # of parameters during burn-in, that row is removed from all of the 
        # matrices used to calculate the next generation, to stop simulating
        # parameter sets which have already reached equilibrium.
        # Therefore, after the burn-in generations are complete, all of the matrices have
        # to be 'recreated' to simulate all parameter combinations after burn-in
        if (cur_gen == burnin) {
            
            # output all rows of the parameter and (epi)allele frequency matrices 
            # that did not reach equilibrium to out_mat
            out_mat[!cumulative_output_rows,] = cbind(param_mat[,,drop = F],
                                                      freq_mat[,,drop = F],
                                                      rep(NA,sum(!cumulative_output_rows)))
            # output all rows of the matrix recording when equilibrium was reached
            # to cumulative_equilibrium_mat
            cumulative_equilibrium_mat[!cumulative_output_rows,] = equilibrium_mat
            
            # recreate (epi)allele frequency, parameter and equilibrium matricies from
            # output matricies
            freq_mat = out_mat[,colnames(freq_mat),drop = F]
            param_mat = out_mat[,colnames(param_mat),drop = F]
            equilibrium_mat = cumulative_equilibrium_mat
            
            # reset output matrix
            out_mat[,] = 0
            out_mat[,"equil_gen"] = NA
            
            # recreate matrices which are subsets of param_mat
            # and output row vectors
            selection_mat = param_mat[,index_list[["freq2fit"]],drop = F]
            c_mat = param_mat[,index_list[["freq2c"]],drop = F]
            i_mat = param_mat[,index_list[["freq2i"]],drop = F]
            cumulative_output_rows = logical(length = nrow(freq_mat))
            output_rows = rep(0,nrow(freq_mat))
            
        }
         
        # this if statement determines when contact of the
        # neutral B allele occurs, it tests if the loop is at
        # inject_delay generations after burn-in
        if(cur_gen == (burnin + inject_delay)) {
            
            # write that contact occurred in the log file
            cat("Worker", i, "neutral locus contact at cur_gen =",
                burnin,"+",inject_delay,"\n")
            
            # contact of the B allele
            freq_mat = injection_migrationC(freq_mat,
                                            m = param_mat[,"m"],
                                            index_list[["bb1_cols"]],
                                            index_list[["BB2_cols"]],
                                            index_list[["bb2_cols"]])
            
        # if not the contact generation, migration occurs as normal
        } else if (((cur_gen < burnin) & pre_contact_mig) |
                   ((cur_gen >= burnin) & post_contact_mig)) {
            
            # applies migration life-history step, see epigen_funcs1.4.cpp for details
            freq_mat = migrationC(freq_mat,
                                  param_mat[,"m"])
            
        }
        
        # applies selection life-history step, see epigen_funcs1.4.cpp for details
        freq_mat = selectionC(freq_mat,
                              selection_mat)
        
        # applies recombination life-history step, see epigen_funcs1.4.cpp for details
        freq_mat = recombinationC(freq_mat,
                                  index_list[["no_reco"]],
                                  index_list[["reco"]],
                                  index_list[["just_r1"]],
                                  index_list[["just_r2"]],
                                  index_list[["both_r"]],
                                  param_mat[,"r1"],
                                  param_mat[,"r2"])
        
        # calculates haplotype frequencies (gamete frequencies), see epigen_funcs1.4.cpp for details
        prop_list = calc_propsC(freq_mat,
                                index_list[["strand1_num"]],
                                index_list[["strand2_num"]])
        
        # applies mating life-history step, see epigen_funcs1.4.cpp for details
        freq_mat = matingC(prop_list,
                           index_list[["vect2"]],
                           index_list[["freq_cols"]])
        
        # applies epimutation life-history step, see epigen_funcs1.4.cpp for details
        freq_mat = epimutationC(freq_mat,
                                c_mat,
                                i_mat,
                                index_list[["ee_cols"]],
                                index_list[["eE_cols"]],
                                index_list[["Ee_cols"]],
                                index_list[["EE_cols"]],
                                index_list[["ee_to_eE_cols"]],
                                index_list[["ee_to_Ee_cols"]],
                                index_list[["ee_to_EE_cols"]],
                                index_list[["eE_to_ee_cols"]],
                                index_list[["eE_to_Ee_cols"]],
                                index_list[["eE_to_EE_cols"]],
                                index_list[["Ee_to_ee_cols"]],
                                index_list[["Ee_to_eE_cols"]],
                                index_list[["Ee_to_EE_cols"]],
                                index_list[["EE_to_ee_cols"]],
                                index_list[["EE_to_eE_cols"]],
                                index_list[["EE_to_Ee_cols"]])
        
        ## this is a debugging function that checks that
        ## (epi)allele frequencies are within the range of possible values
        # freq_check(freq_mat,
        #           all_eq=T,
        #           to_print=F)
        
        ## calculate frequencies of (epi)alleles and determine if equilibrium has been reached
        # if there is a non-zero equilibrium threshold
        if(equilibrium_threshold != 0) {
            
            # if during burn-in, calculate equilibrium of burnin_equilibrium_loci
            if(cur_gen < burnin & length(burnin_equilibrium_loci) != 0) {
                
                # calculate (epi)allele frequencies and check if equilibrium is reached
                equilibrium_mat = multi_equilibrium(burnin_equilibrium_loci)
                # update output_rows vector with rows that have reached equilibrium
                output_rows = apply(equilibrium_mat[,paste0("equil_",burnin_equilibrium_loci),drop = F],1,function(x) {!(0 %in% x)})
                
            # if after burn-in, calculate equilibrium of equilibrium_loci
            } else if (cur_gen >= (burnin+inject_delay) & length(equilibrium_loci) != 0){
                # ^ this else if is so that there is no halting prior to the delayed
                # injection of the neutral B allele
                
                equilibrium_mat = multi_equilibrium(equilibrium_loci)
                output_rows = apply(equilibrium_mat[,paste0("equil_",equilibrium_loci),drop = F],1,function(x) {!(0 %in% x)})
            
            } else {} # nothing happens between burnin and burnin_inject_delay generations
            
        }
        
        ## output rows that have reached equilibrium
        if(sum(output_rows) > 0) {
            
            # use the index of all rows that have already been output (cumulative_output_rows)
            # and the index of rows to be output this generation (output_rows)
            # to find where to output rows into the out_mat this generation
            output_rows_to_out_mat = rep(F,length(cumulative_output_rows))
            output_rows_to_out_mat[!cumulative_output_rows][output_rows] = TRUE
            
            # update cumulative_output_rows with the rows output this generation
            cumulative_output_rows[!cumulative_output_rows][output_rows] = TRUE
            
            # output from parameter sets that have reached equilibrium to out_mat
            # and cumulative_equilibrium_mat
            out_mat[output_rows_to_out_mat,] = cbind(param_mat[output_rows,,drop = F],
                                                     freq_mat[output_rows,,drop = F],
                                                     rep(cur_gen,sum(output_rows)))
            cumulative_equilibrium_mat[output_rows_to_out_mat,] = equilibrium_mat[output_rows,]
            
            # if during burn-in generations save the generation which equilibrium was reach
            # in the burnin_equil_gen vector. This is done because out_mat gets overwriten
            # after burn-in is complete and we want to save these values to output at the end
            # of the simulation
            if(cur_gen < burnin) {
                burnin_equil_gen[output_rows_to_out_mat] = out_mat[output_rows_to_out_mat,"equil_gen"]
            }
            
            # update model matrices to take out rows that have reached equilibrium
            freq_mat = freq_mat[!output_rows,,drop = F]
            prop_list$haplo_freq_pop1 = prop_list$haplo_freq_pop1[!output_rows,,drop = F]
            prop_list$haplo_freq_pop2 = prop_list$haplo_freq_pop2[!output_rows,,drop = F]
            param_mat = param_mat[!output_rows,,drop = F]
            equilibrium_mat = equilibrium_mat[!output_rows,,drop = F]
            selection_mat = selection_mat[!output_rows,,drop = F]
            c_mat = c_mat[!output_rows,,drop = F]
            i_mat = i_mat[!output_rows,,drop = F]
            
            # test if all parameter sets have finished
            if (nrow(freq_mat) == 0 & nrow(param_mat) == 0) {
                
                # if burn-in has finished, jump to burnin generations
                # to start simulating post-contact
                if (cur_gen < burnin) {
                    
                    # set current generation to the max burn-in generations
                    cur_gen = burnin - 1 # -1 since it is incremented before the next loop iteration
                
                # otherwise end the simulation  
                } else {
                    
                    # calculate RI, see epigen_funcs1.4.cpp for details
                    # ncol(param_mat) is added to indexes since the indexes
                    # give the columns of freq_mat but out_mat has param_mat columns
                    # followed by freq_mat columns
                    RI = RI_calcC(out_mat,
                                  out_mat[,"m"],
                                  index_list[["BB1_cols"]]+ncol(param_mat),
                                  index_list[["BB2_cols"]]+ncol(param_mat),
                                  index_list[["B_het1_cols"]]+ncol(param_mat),
                                  index_list[["B_het2_cols"]]+ncol(param_mat))
                    
                    # add RI and burnin_equil_gen to out_mat
                    out_mat = cbind(out_mat,RI,burnin_equil_gen)
                    colnames(out_mat)[ncol(out_mat)] = "burnin_equil_gen"
                    
                    # if recording over time, add recorded values to out_mat
                    if(record_gen != 0) {
                        out_mat = cbind(out_mat,record_freqs_mat)
                    }
                    
                    return(out_mat) # pass out_mat back from the worker
                    break("Worker",i,"finished at cur_gen =",cur_gen) # stop this worker
                    
                }
            }
        }
                                
        # record frequencies each record_gen generations
        # only after contact of the neutral B allele
        if(record_gen != 0 & cur_gen >= (burnin+inject_delay)) {
            
            # record less frequently in later generations to avoid making too large an output file
            if(cur_gen == (burnin + inject_delay + 100) & record_gen < 10) {thin_record_gen = 10}
            if(cur_gen == (burnin + inject_delay + 1000) & record_gen < 100) {thin_record_gen = 100}
            if(cur_gen == (burnin + inject_delay + 10000) & record_gen < 1000) {thin_record_gen = 1000}
            
            if(cur_gen %% thin_record_gen == 0) {
                
                # since (epi)alle frequencies are already calculated for
                # testing if equilibrium is reached, (epi)allele frequencies are only
                # calculated here, if they were not already calculated
                # either since equilibrium is not being calculated (equilibrium_threshold == 0)
                # or there are loci whoes frequencies are being recorded but whos frequencies
                # are not being tested for equilibrium
                if(equilibrium_threshold == 0 |
                   sum(record_loci %in% equilibrium_loci) != length(record_loci)) {
                    unique_record_loci = record_loci[!(record_loci %in% equilibrium_loci)]
                    equilibrium_mat = multi_equilibrium(unique_record_loci)
                }
                
                # save (epi)allele frequencies
                record_freqs_mat = multi_record_gen(record_loci,
                                                    record_freqs_mat,
                                                    equilibrium_mat)
                
                # calculate and save LD metrics
                record_freqs_mat = calc_LD(prop_list,
                                           record_LD,
                                           record_freqs_mat,
                                           index_list[["LD_list"]])
                
            }
        }
         
        cur_gen = cur_gen + 1 # increment the number of generations
         
    }
    
    # output any remaining parameter sets that have not reached 
    # equilibrium by the max generation number
    out_mat[!cumulative_output_rows,] = cbind(param_mat[,,drop = F],
                                                    freq_mat[,,drop = F],
                                                    rep(cur_gen,nrow(param_mat)))
    
    # the functions applied here are the same as those applied when 
    # all parameter sets have reached equilibrium
    
    RI = RI_calcC(out_mat,
                  out_mat[,"m"],
                  index_list[["BB1_cols"]]+ncol(param_mat),
                  index_list[["BB2_cols"]]+ncol(param_mat),
                  index_list[["B_het1_cols"]]+ncol(param_mat),
                  index_list[["B_het2_cols"]]+ncol(param_mat))
    
    out_mat = cbind(out_mat,RI,burnin_equil_gen)
    colnames(out_mat)[ncol(out_mat)] = "burnin_equil_gen"
    
    if(record_gen != 0) {
        out_mat = cbind(out_mat,record_freqs_mat)
    }
    
    return(out_mat)
    # write that the worker is done to the log file
    cat("Worker",i,"finished at cur_gen =",cur_gen,"\n")
    
}

cat("\n## Stopping cluster ##\n\n")
parallel::stopCluster(cl = my.cluster)

#### Merge Output List to a data.frame ####

## count the number of columns in each output matrix
## and save the largest number of columns.
# Since the simulations are halted when all parameter sets
# have reached equilibrium, output matricies whose parameter sets
# reach equilibrium faster will have fewer columns. So, we save
# the number of columns of the largest output matrix and pad
# the smaller matrices with NAs so that they have the same number
# of columns and all the output matricies can be merged
ncol_vect = sapply(out_mat_list,ncol)
biggest_out_mat = out_mat_list[[which(ncol_vect == max(ncol_vect))[1]]]

for(i in 1:length(out_mat_list)) {
    append_col_count = ncol(biggest_out_mat) - ncol_vect[i] 
    temp_matrix = matrix(NA,nrow(biggest_out_mat),append_col_count)
    out_mat_list[[i]] = cbind(out_mat_list[[i]],temp_matrix)
    colnames(out_mat_list[[i]]) = colnames(biggest_out_mat)
}

# merge the output matricies
out_mat = do.call(rbind,out_mat_list)

#### Write output to .csv ####

cat("Ouputting results to: ",getwd(),out_csv_file,"\n",sep = "")
write.csv(out_mat,out_csv_file)

## Output runtime to log-file and close log file
cat("Runtime:\n")
cat(proc.time(),"\n")
sink()

