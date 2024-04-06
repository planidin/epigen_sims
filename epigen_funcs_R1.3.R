
##############################################
## v1.3 - Added indexing for calculating LD ##
## converted e_X to c_X                     ##
##############################################

#### Load Packages ####

library(EpiSims)

#### Data Matrix Construction Functions ####

# creates a vector of (epi)genotypes using
# some string manipulation
make_allele_combos = function() {
    
    # create a vector of the possible haplotypes
    a_vect = rep(c("a","A"), each = 4,times = 1)
    e_vect = rep(c("e","E"), each = 2,times = 2)
    b_vect = rep(c("b","B"), each = 1,times = 4)
    allele_vect = paste0(a_vect,e_vect,b_vect)
    
    allele_mat = NULL
    
    # calculate all possible diploid combinations
    for(i in allele_vect) {
        
        pairs = NULL
        
        # pair the current haplotype (strand1, plus)
        # with each other haplotye (strand2, minus)
        for(j in 1:nchar(i)) {
            
            plus = substr(i,j,j)
            minus = substr(allele_vect,j,j)
            pairs = cbind(pairs,paste0(plus,minus)) # append the combo to the pairs matrix
            
        }
        
        # turn each row in the pairs matrix into a single string
        allele_col = apply(pairs,1,function(x) paste(x,collapse = "_"))
        
        # append the string to the (epi)alleles matrix
        allele_mat = cbind(allele_mat,allele_col)
        
    }
    
    # take just the lower triangular portion of the 
    # allele combos matrix since it its redundant with the
    # upper trangular portion (we don't care which strand is stand 1)
    allele_combos = allele_mat[lower.tri(allele_mat, diag = T)]
    
    return(allele_combos)
    
}

# calculate the matrix of fitnesses based on selection parameters
add_fitness = function(param_mat,
                       allele_combos,
                       s = param_mat[,c("s.a","s.e")],
                       h = param_mat[,c("h.a","h.e")],
                       epistat = param_mat[,"epistat"]) {
    
    # get the possible combinations of loci under selection 
    # by removing the B locus from the vector of possible 
    # (epi)genotypes and keeping unique combinations after
    fitness_combos = gsub("_bb","",allele_combos,ignore.case = T)
    fitness_combos = unique(fitness_combos)
    
    # uppercase alleles are always favored in pop1
    # lowercase alleles are always favored in pop2
    
    sh = s*h # calulate heterozygote fitness
    
    # get the locus names from the fitness and dominance columns
    # e.g. s.e tells us its the fitness of e individuals in deme 1
    s_loci = sapply(1:length(strsplit(colnames(s),"\\.")), function(x) strsplit(colnames(s),"\\.")[[x]][2])
    h_loci = sapply(1:length(strsplit(colnames(h),"\\.")), function(x) strsplit(colnames(h),"\\.")[[x]][2])
    
    # check that the names are valid
    # then create the alternate (epi)allele with toupper()
    if(sum(s_loci != h_loci) > 0) {
        stop("Selection and dominance loci don't match")
    } else if (sum(nchar(s_loci) - 1) != 0) {
        stop("Locus symbol is something other than a single letter")
    } else {
        allele1 = s_loci
        allele2 = toupper(s_loci)
    }
    
    # calculate the fitness for each unique possible fitness
    for(i in fitness_combos) {
        
        # lc stands for lowercase
        
        # calculate the number of lowercase letters (allele1) for each
        # (epi)gentype, for calculating epistasis
        num_lc = sapply(s_loci, function(x) sum(gregexpr(x,i)[[1]] != -1))
        num_lc = as.data.frame(t(num_lc))
        
        # create fitness data.frame filled with zeros
        pop1_fit = numeric(length(s[,1]))
        pop2_fit = numeric(length(s[,1]))
        
        # calculate additive fitness with epistasis
        # with epistat = 0, fitness is additive within and between loci
        # epistasis is calculated as epistat times the number of
        # interactions between maladapted alleles 
        # e.g. an aa_ee individual has two a alleles each interacting
        # with two e epialleles, so there are 4 epistatic interactions
        if(num_lc$a == 2 & num_lc$e == 2) {
            pop1_fit = 1 - (s[,"s.a"] + s[,"s.e"] + 4*epistat)
            pop2_fit = 1
        } else if (num_lc$a == 1 & num_lc$e == 2) {
            pop1_fit = 1 - (sh[,"s.a"] + s[,"s.e"] + 2*epistat)
            pop2_fit = 1 - (sh[,"s.a"] - 2*epistat)
        } else if (num_lc$a == 0 & num_lc$e == 2) {
            pop1_fit = 1 - s[,"s.e"]
            pop2_fit = 1 - s[,"s.a"]
        } else if (num_lc$a == 2 & num_lc$e == 1) {
            pop1_fit = 1 - (s[,"s.a"] + sh[,"s.e"] + 2*epistat)
            pop2_fit = 1 - (sh[,"s.e"] - 2*epistat)
        } else if (num_lc$a == 1 & num_lc$e == 1) {
            pop1_fit = 1 - (sh[,"s.a"] + sh[,"s.e"] + epistat)
            pop2_fit = 1 - (sh[,"s.a"] + sh[,"s.e"] - epistat)
        } else if (num_lc$a == 0 & num_lc$e == 1) {
            pop1_fit = 1 - sh[,"s.e"]
            pop2_fit = 1 - (s[,"s.a"] + sh[,"s.e"])
        } else if (num_lc$a == 2 & num_lc$e == 0) {
            pop1_fit = 1 - s[,"s.a"]
            pop2_fit = 1 - s[,"s.e"]
        } else if (num_lc$a == 1 & num_lc$e == 0) {
            pop1_fit = 1 - sh[,"s.a"]
            pop2_fit = 1 - (sh[,"s.a"] + s[,"s.e"])
        } else if (num_lc$a == 0 & num_lc$e == 0) {
            pop1_fit = 1
            pop2_fit = 1 - (s[,"s.a"] + s[,"s.e"])
        }
        
        # merge values for each population into one matrix
        # and append this to the parameter matrix
        fit_cols = cbind(pop1_fit,pop2_fit)
        colnames(fit_cols) = paste0("s_",i,"_pop",c(1,2))
        param_mat = cbind(param_mat,fit_cols)
        
    }
    
    return(param_mat)
}

# calculate transition rates for each genotype at the A
# locus, depending on base epimutation rate, 
# genetic effect (X_gen_Y) and dominance (h.X_gen_Y)
add_epiparams = function(param_mat,
                         i_1 = param_mat[,"i_1"],
                         i_2 = param_mat[,"i_2"],
                         c_1 = param_mat[,"c_1"],
                         c_2 = param_mat[,"c_2"],
                         i_gen_1 = param_mat[,"i_gen_1"],
                         i_gen_2 = param_mat[,"i_gen_2"],
                         c_gen_1 = param_mat[,"c_gen_1"],
                         c_gen_2 = param_mat[,"c_gen_2"],
                         h.i_gen_1 = param_mat[,"h.i_gen_1"],
                         h.i_gen_2 = param_mat[,"h.i_gen_2"],
                         h.c_gen_1 = param_mat[,"h.c_gen_1"],
                         h.c_gen_2 = param_mat[,"h.c_gen_2"]) {
    
    # when a gen effect is zero, the transition rate is the same
    # for all genotypes
    
    # otherwise the A allele can increase a given 
    # transition rate, calulated like additive fitness
    
    i_aa_1 = i_1*(1-i_gen_1)
    i_aA_1 = i_1*(1-i_gen_1) + i_1*i_gen_1*h.i_gen_1
    i_AA_1 = i_1*(1-i_gen_1) + i_1*i_gen_1
    
    i_aa_2 = i_2*(1-i_gen_2)
    i_aA_2 = i_2*(1-i_gen_2) + i_2*i_gen_2*h.i_gen_2
    i_AA_2 = i_2*(1-i_gen_2) + i_2*i_gen_2
    
    c_aa_1 = c_1*(1-c_gen_1)
    c_aA_1 = c_1*(1-c_gen_1) + c_1*c_gen_1*h.c_gen_1
    c_AA_1 = c_1*(1-c_gen_1) + c_1*c_gen_1
    
    c_aa_2 = c_2*(1-c_gen_2)
    c_aA_2 = c_2*(1-c_gen_2) + c_2*c_gen_2*h.c_gen_1
    c_AA_2 = c_2*(1-c_gen_2) + c_2*c_gen_2
    
    # merge the values together and append to the parameter matrix
    param_mat = cbind(param_mat,
                      i_aa_1,i_aA_1,i_AA_1,
                      i_aa_2,i_aA_2,i_AA_2,
                      c_aa_1,c_aA_1,c_AA_1,
                      c_aa_2,c_aA_2,c_AA_2)
    
}

# create the (epi)genotype frequency matrix with starting frequencies
allele_start = function(n_param_rows,
                        allele_combos = allele_combos) {
    
    # start both populations with just bb with the highest fitness
    # (epi)genotype
    
    # create a matrix of zeros to store (epi)genotype frequencies
    # based on the vector of possible (epi)genotypes
    freq_mat = matrix(0,n_param_rows,2*length(allele_combos))
    # create column names for deme 1 and deme 2
    colnames(freq_mat) = c(paste("freq",allele_combos,1,sep = "_"),
                           paste("freq",allele_combos,2,sep = "_"))
    
    # fix each deme for the highest fitness (epi)genotype
    freq_mat[,"freq_AA_EE_bb_1"] = 1 # favored alleles pop 1
    freq_mat[,"freq_aa_ee_bb_2"] = 1 # favored alleles pop 2
    
    return(freq_mat)
    
}

#### Trouble Shooting Functions ####

# this function checks that frequencies in deme 1 and 2 are adding up to 1
freq_check = function(IN, all_eq = F, to_print = T) {
    
    # calculate the total frequency in each deme
    tot_freq_1 = apply(matrix(IN[,colnames(IN)[grepl("freq_",colnames(IN)) &
                                               grepl("1",colnames(IN))]],ncol = 36),1,sum)
    tot_freq_2 = apply(matrix(IN[,colnames(IN)[grepl("freq_",colnames(IN)) &
                                               grepl("2",colnames(IN))]],ncol = 36),1,sum)
    
    # prints out the frequencies with many digits to check for
    # floating point errors
    if(to_print) {
        cat(sprintf("%.54f\n",head(tot_freq_1)))
        cat(sprintf("%.54f\n",head(tot_freq_2)))
    } else {}
    
    # use all.equal instead of == due to rounding errors
    # in the floating point arithmetic making the totals
    # slightly different than one, even when everything is
    # working properly
    if(all_eq) {
        
        # if not summing to 1, stop the simulation
        if(!all.equal(sum(tot_freq_1),length(tot_freq_1))) {
            stop("allele frequencies in pop 1 do not sum to 1!")
        }
        if(!all.equal(sum(tot_freq_2),length(tot_freq_2))) {
            stop("allele frequencies in pop 2 do not sum to 1!")
        }
        
    } else {
        
        # if not summing to 1, stop the simulation
        if(sum(tot_freq_1 != 1) > 0) {
            stop("allele frequencies in pop 1 do not sum to 1!")
        }
        if(sum(tot_freq_2 != 1) > 0) {
            stop("allele frequencies in pop 2 do not sum to 1!")
        }
        
    }
    
}

#### String Manipulation Functions ####

# this function reverses the strand order of column names for freq_mat
# for example freq_Aa_Ee_bB_2 is converted to freq_aA_eE_Bb_2
strand_reverse = function(x) {
 
    # index string positions, swapping the strands for each locus
    # freq_Aa_Ee_bB_2 is 15 characters long
    flip_seq = c(1:5,7,6,8,10,9,11,13,12,14,15)
    out_vect = character(length = 15)
    a = 1
    for(i in flip_seq) {
        out_vect[a] = substr(x,i,i) # fill a vector with the characters in the new sequence
        a = a + 1
    }
    
    # collapse the vector into a string
    out_string = paste(out_vect, sep="", collapse="") 
    
    return(out_string)
    
}

# like strand reverse but only applied to a single locus
locus_reco = function(x,locus) {
    
    # set position to flip according to selected locus
    if(locus == "A") {
        start_pos = 6
    } else if(locus == "E") {
        start_pos = 9
    } else if(locus == "B") {
        start_pos = 12
    } else {
        stop("Invalid locus provided, must be A E or B")
    }
    
    flip_seq = c(1:(start_pos-1),
                 start_pos+1,
                 start_pos,
                 (start_pos+2):15)
    out_vect = character(length = 15)
    a = 1
    for(i in flip_seq) {
        out_vect[a] = substr(x,i,i)
        a = a + 1
    }
    
    out_string = paste(out_vect, sep="", collapse="") 
    
    return(out_string)
    
}

#### Index Creation Functions ####

# create the freq2fit index to use in selectionC function
# this tells which columns in the fitness matrix
# correspond to which columns in the frequency matrix
make_selection_index = function(fit_cols,freq_cols) {
    
    # intitalize the freq2fit map with frequency column names
    freq2fit = character(length = length(freq_cols))
    names(freq2fit) = freq_cols
    
    for(i in fit_cols) {
        
        # for each fitness column, pull out the substring
        # corresponding to A and E alleles e.g. Aa_EE
        sel = gsub("s_","",i)
        sel = gsub("_pop.","",sel)
        # get the deme
        if(grepl("pop1",i)) {
            pop = 1
        } else {
            pop = 2
        }
        
        # if a frequency column name matches a fitness allele combo and deme
        # replace that column name with the fitness column name
        # this creates a fitness for each frequency column
        freq2fit[grepl(sel,freq_cols) & grepl(pop,freq_cols)] = fit_cols[grepl(sel,fit_cols) & grepl(pop,fit_cols)]
    }
    
    # output the freq2fit vector as a list to merge with other indicies later
    out_list = list(freq2fit)
    names(out_list) = "freq2fit"
    
    return(out_list)
    
}

# create indicies to use in recombinationC using string manipulation
make_recombination_index = function(freq_cols,locus_order) {
    
    # swap strands for locus 1
    # new (epi)genotype tahat are not a valid (epi)genotype e.g. aA_ee_Bb
    # swap trands for all loci to change the "reference" strand
    # e.g. aA_ee_Bb becomes Aa_ee_bB
    just_r1 = sapply(freq_cols,function(x) locus_reco(x,locus_order[1]))
    just_r1[!(just_r1 %in% freq_cols)] = sapply(just_r1[!(just_r1 %in% freq_cols)],strand_reverse)
    
    # apply the same to locus 3
    just_r2 = sapply(freq_cols,function(x) locus_reco(x,locus_order[3]))
    just_r2[!(just_r2 %in% freq_cols)] = sapply(just_r2[!(just_r2 %in% freq_cols)],strand_reverse)
  
    # apply the same to locus 2
    both_r = sapply(freq_cols,function(x) locus_reco(x,locus_order[2]))
    both_r[!(both_r %in% freq_cols)] = sapply(both_r[!(both_r %in% freq_cols)],strand_reverse)
    
    # NOTE: this method is slightly inefficient since each event i.e. r1, r2 and both_r
    # only applies to 24 genotypes each. However it is easier to apply them all to the
    # 32 genotypes where at least one event effects that genotype, rather than making
    # unique indexes for each
    
    # get columns which are not affected by any recombination
    no_reco = freq_cols[(just_r1 == freq_cols) &
                        (just_r2 == freq_cols) &
                        (both_r == freq_cols)]
    
    # get columns that are affected by at least one recombination event
    reco = freq_cols[(just_r1 != freq_cols) |
                     (just_r2 != freq_cols) |
                     (both_r != freq_cols)]
    
    # remove columns that are unaffected by recombination from the
    # vectors of post-recombination (epi)genotypes
    just_r1 = just_r1[freq_cols %in% reco]
    just_r2 = just_r2[freq_cols %in% reco]
    both_r = both_r[freq_cols %in% reco]
    
    # inititate output vectors
    reco_num = just_r1_num = just_r2_num = both_r_num = numeric(length = length(reco))
    
    # find the column in the frequency matrix that correspond to the
    # recombination event, converting to numeric indices for C++
    for(i in 1:length(reco)) {
        # - 1 for C++ indexing
        reco_num[i] = which(freq_cols == reco[i]) - 1
        just_r1_num[i] = which(freq_cols == just_r1[i]) - 1
        just_r2_num[i] = which(freq_cols == just_r2[i]) - 1
        both_r_num[i] = which(freq_cols == both_r[i]) - 1
    }
    
    # get the indices of column not affected by recombination
    no_reco_num = numeric(length = length(no_reco))
    
    for(i in 1:length(no_reco)) {
        no_reco_num[i] = which(freq_cols == no_reco[i]) - 1
    }
    
    # merge indicies into a list 
    out_list = list(no_reco_num,reco_num,just_r1_num,just_r2_num,both_r_num)
    names(out_list) = c("no_reco","reco","just_r1","just_r2","both_r")
    
    return(out_list)
    
}

# make indicies to use in epimutationC function
make_epimutation_index = function(freq_cols) {
    
    # get the frequency matrix columns corresponding to each
    # epigenotype
    ee_cols = freq_cols[grepl("ee",freq_cols)]
    eE_cols = freq_cols[grepl("eE",freq_cols)]
    Ee_cols = freq_cols[grepl("Ee",freq_cols)]
    EE_cols = freq_cols[grepl("EE",freq_cols)]
    
    # convert the column names based on each
    # possible transition in epigenetic state
    ee_to_eE_cols = gsub("ee","eE",ee_cols)
    ee_to_Ee_cols = gsub("ee","Ee",ee_cols)
    ee_to_EE_cols = gsub("ee","EE",ee_cols)
    
    eE_to_ee_cols = gsub("eE","ee",eE_cols)
    eE_to_Ee_cols = gsub("eE","Ee",eE_cols)
    eE_to_EE_cols = gsub("eE","EE",eE_cols)
    
    Ee_to_ee_cols = gsub("Ee","ee",Ee_cols)
    Ee_to_eE_cols = gsub("Ee","eE",Ee_cols)
    Ee_to_EE_cols = gsub("Ee","EE",Ee_cols)
    
    EE_to_ee_cols = gsub("EE","ee",EE_cols)
    EE_to_eE_cols = gsub("EE","eE",EE_cols)
    EE_to_Ee_cols = gsub("EE","Ee",EE_cols)
    
    # flipping the post-epimutation indicides which are not in freq_cols
    # to make them valid
    ee_to_Ee_cols[!(ee_to_Ee_cols %in% freq_cols)] = sapply(ee_to_Ee_cols[!(ee_to_Ee_cols %in% freq_cols)],strand_reverse)
    eE_to_ee_cols[!(eE_to_ee_cols %in% freq_cols)] = sapply(eE_to_ee_cols[!(eE_to_ee_cols %in% freq_cols)],strand_reverse)
    eE_to_Ee_cols[!(eE_to_Ee_cols %in% freq_cols)] = sapply(eE_to_Ee_cols[!(eE_to_Ee_cols %in% freq_cols)],strand_reverse)
    eE_to_EE_cols[!(eE_to_EE_cols %in% freq_cols)] = sapply(eE_to_EE_cols[!(eE_to_EE_cols %in% freq_cols)],strand_reverse)
    EE_to_Ee_cols[!(EE_to_Ee_cols %in% freq_cols)] = sapply(EE_to_Ee_cols[!(EE_to_Ee_cols %in% freq_cols)],strand_reverse)
    
    # converting to numeric C++ indicies
    ee_cols_num = numeric() 
    eE_cols_num = numeric() 
    Ee_cols_num = numeric() 
    EE_cols_num = numeric() 
    ee_to_eE_cols_num = numeric() 
    ee_to_Ee_cols_num = numeric() 
    ee_to_EE_cols_num = numeric() 
    eE_to_ee_cols_num = numeric() 
    eE_to_Ee_cols_num = numeric() 
    eE_to_EE_cols_num = numeric() 
    Ee_to_ee_cols_num = numeric() 
    Ee_to_eE_cols_num = numeric() 
    Ee_to_EE_cols_num = numeric() 
    EE_to_ee_cols_num = numeric() 
    EE_to_eE_cols_num = numeric() 
    EE_to_Ee_cols_num = numeric()
    
    # did this weird if statement thing to avoid making different loops
    # for each different length of vector, 24 is the longest index vector
    for(i in 1:24) {
        if(i <= length(ee_cols)) {ee_cols_num[i] = which(freq_cols == ee_cols[i]) - 1}
        if(i <= length(eE_cols)) {eE_cols_num[i] = which(freq_cols == eE_cols[i]) - 1}
        if(i <= length(Ee_cols)) {Ee_cols_num[i] = which(freq_cols == Ee_cols[i]) - 1}
        if(i <= length(EE_cols)) {EE_cols_num[i] = which(freq_cols == EE_cols[i]) - 1}
        if(i <= length(ee_to_eE_cols)) {ee_to_eE_cols_num[i] = which(freq_cols == ee_to_eE_cols[i]) - 1}
        if(i <= length(ee_to_Ee_cols)) {ee_to_Ee_cols_num[i] = which(freq_cols == ee_to_Ee_cols[i]) - 1}
        if(i <= length(ee_to_EE_cols)) {ee_to_EE_cols_num[i] = which(freq_cols == ee_to_EE_cols[i]) - 1}
        if(i <= length(eE_to_ee_cols)) {eE_to_ee_cols_num[i] = which(freq_cols == eE_to_ee_cols[i]) - 1}
        if(i <= length(eE_to_Ee_cols)) {eE_to_Ee_cols_num[i] = which(freq_cols == eE_to_Ee_cols[i]) - 1}
        if(i <= length(eE_to_EE_cols)) {eE_to_EE_cols_num[i] = which(freq_cols == eE_to_EE_cols[i]) - 1}
        if(i <= length(Ee_to_ee_cols)) {Ee_to_ee_cols_num[i] = which(freq_cols == Ee_to_ee_cols[i]) - 1}
        if(i <= length(Ee_to_eE_cols)) {Ee_to_eE_cols_num[i] = which(freq_cols == Ee_to_eE_cols[i]) - 1}
        if(i <= length(Ee_to_EE_cols)) {Ee_to_EE_cols_num[i] = which(freq_cols == Ee_to_EE_cols[i]) - 1}
        if(i <= length(EE_to_ee_cols)) {EE_to_ee_cols_num[i] = which(freq_cols == EE_to_ee_cols[i]) - 1}
        if(i <= length(EE_to_eE_cols)) {EE_to_eE_cols_num[i] = which(freq_cols == EE_to_eE_cols[i]) - 1}
        if(i <= length(EE_to_Ee_cols)) {EE_to_Ee_cols_num[i] = which(freq_cols == EE_to_Ee_cols[i]) - 1}
    }
    
    # get columns to make a matrix of genotype dependent transitions rates 
    # with X_gen_Y paramaters e.g. converts freq_aA_EE_bB_1 to 
    # c_aA_1 in freq2c and i_aA_1 in freq2i
    freq2c = paste0("c_",substr(freq_cols,6,7),substr(freq_cols,14,15))
    freq2i = paste0("i_",substr(freq_cols,6,7),substr(freq_cols,14,15))
    
    # merge indices into a list and output
    out_list = list(ee_cols_num,
                    eE_cols_num,
                    Ee_cols_num,
                    EE_cols_num,
                    ee_to_eE_cols_num,
                    ee_to_Ee_cols_num,
                    ee_to_EE_cols_num,
                    eE_to_ee_cols_num,
                    eE_to_Ee_cols_num,
                    eE_to_EE_cols_num,
                    Ee_to_ee_cols_num,
                    Ee_to_eE_cols_num,
                    Ee_to_EE_cols_num,
                    EE_to_ee_cols_num,
                    EE_to_eE_cols_num,
                    EE_to_Ee_cols_num,
                    freq2c,
                    freq2i)
    
    names(out_list) = c("ee_cols",
                        "eE_cols",
                        "Ee_cols",
                        "EE_cols",
                        "ee_to_eE_cols",
                        "ee_to_Ee_cols",
                        "ee_to_EE_cols",
                        "eE_to_ee_cols",
                        "eE_to_Ee_cols",
                        "eE_to_EE_cols",
                        "Ee_to_ee_cols",
                        "Ee_to_eE_cols",
                        "Ee_to_EE_cols",
                        "EE_to_ee_cols",
                        "EE_to_eE_cols",
                        "EE_to_Ee_cols",
                        "freq2c",
                        "freq2i")
    
    return(out_list)
    
}

# make a list of indicies to be used by all functions in the simulation
# this function uses the above index creation functions plus makes some
# simplier indicies, merging everything into one big list
make_index_list = function(param_mat,
                           freq_mat,
                           locus_order = locus_order) {
    
    # get column names from the frequency and parameter matrix
    freq_cols = colnames(freq_mat)[grepl("freq_",colnames(freq_mat))]
    freq_cols1 = freq_cols[grepl("1",freq_cols)]
    freq_cols2 = freq_cols[grepl("2",freq_cols)]
    fit_cols = colnames(param_mat)[grepl("s_",colnames(param_mat))]
    
    # get columns that correspond to each genotype at the B locus
    bb1_cols = freq_cols1[grep("bb",freq_cols1)]
    BB1_cols = freq_cols1[grep("BB",freq_cols1)]
    bb2_cols = freq_cols2[grep("bb",freq_cols2)]
    BB2_cols = freq_cols2[grep("BB",freq_cols2)]
    B_het1_cols = freq_cols1[!(freq_cols1 %in% bb1_cols) & !(freq_cols1 %in% BB1_cols)]
    B_het2_cols = freq_cols2[!(freq_cols2 %in% bb2_cols) & !(freq_cols2 %in% BB2_cols)]
    
    # get columns that correspond to each genotype at the A locus
    aa1_cols = freq_cols1[grep("aa",freq_cols1)]
    AA1_cols = freq_cols1[grep("AA",freq_cols1)]
    aa2_cols = freq_cols2[grep("aa",freq_cols2)]
    AA2_cols = freq_cols2[grep("AA",freq_cols2)]
    A_het1_cols = freq_cols1[!(freq_cols1 %in% aa1_cols) & !(freq_cols1 %in% AA1_cols)]
    A_het2_cols = freq_cols2[!(freq_cols2 %in% aa2_cols) & !(freq_cols2 %in% AA2_cols)]
    
    # get columns that correspond to each epigenotype at the E locus
    ee1_cols = freq_cols1[grep("ee",freq_cols1)]
    EE1_cols = freq_cols1[grep("EE",freq_cols1)]
    ee2_cols = freq_cols2[grep("ee",freq_cols2)]
    EE2_cols = freq_cols2[grep("EE",freq_cols2)]
    E_het1_cols = freq_cols1[!(freq_cols1 %in% ee1_cols) & !(freq_cols1 %in% EE1_cols)]
    E_het2_cols = freq_cols2[!(freq_cols2 %in% ee2_cols) & !(freq_cols2 %in% EE2_cols)]
    
    # convert colnames to numerics, -1 for C++ indexing
    bb1_cols = which(freq_cols %in% bb1_cols) - 1
    BB1_cols = which(freq_cols %in% BB1_cols) - 1
    bb2_cols = which(freq_cols %in% bb2_cols) - 1
    BB2_cols = which(freq_cols %in% BB2_cols) - 1
    B_het1_cols = which(freq_cols %in% B_het1_cols) - 1
    B_het2_cols = which(freq_cols %in% B_het2_cols) - 1
    
    aa1_cols = which(freq_cols %in% aa1_cols) - 1
    AA1_cols = which(freq_cols %in% AA1_cols) - 1
    aa2_cols = which(freq_cols %in% aa2_cols) - 1
    AA2_cols = which(freq_cols %in% AA2_cols) - 1
    A_het1_cols = which(freq_cols %in% A_het1_cols) - 1
    A_het2_cols = which(freq_cols %in% A_het2_cols) - 1
    
    ee1_cols = which(freq_cols %in% ee1_cols) - 1
    EE1_cols = which(freq_cols %in% EE1_cols) - 1
    ee2_cols = which(freq_cols %in% ee2_cols) - 1
    EE2_cols = which(freq_cols %in% EE2_cols) - 1
    E_het1_cols = which(freq_cols %in% E_het1_cols) - 1
    E_het2_cols = which(freq_cols %in% E_het2_cols) - 1
    
    # create a vector for duplicate pairing
    # of non-identical gametes
    # VectOf2s in form_gametesC()
    vect2 = numeric()
    for(i in 8:1) {
        vect2 = c(vect2,1,rep(2,i-1))
    }
    vect2 = c(vect2,vect2)
    
    # the the string corresponding to each strand and convert to
    # a numeric index
    strand1 = paste0(substr(freq_cols1,6,6),substr(freq_cols1,9,9),substr(freq_cols1,12,12))
    strand2 = paste0(substr(freq_cols1,7,7),substr(freq_cols1,10,10),substr(freq_cols1,13,13))
    strand1_num = as.numeric(as.factor(strand1)) - 1 # -1 for C++ indexing
    strand2_num = as.numeric(as.factor(strand2)) - 1 # -1 for C++ indexing
    
    # create a list of allele pairings for LD calculation
    # unique(strand1) and unique(strand2) are the same, so using strand1
    LD_list = list()
    for(i in c("A","E","B","a","e","b")) {
        LD_list[[length(LD_list) + 1]] = grep(i,unique(strand1))
        names(LD_list)[length(LD_list)] = i
        for(j in c("A","E","B","a","e","b")) {
            if(tolower(i) != tolower(j)) {
                LD_list[[length(LD_list) + 1]] = which(grepl(i,unique(strand1)) & grepl(j,unique(strand1)))
                names(LD_list)[length(LD_list)] = paste0(i,j)
            }
        }
    }
    
    # create index list from other functions
    index1 = make_epimutation_index(freq_cols)
    index2 = make_selection_index(fit_cols,freq_cols)
    index3 = make_recombination_index(freq_cols,locus_order)
    
    # merge all indicies made in this function into one list and name it
    out_list = list(freq_cols, freq_cols1, freq_cols2, fit_cols,
                    bb1_cols, BB1_cols, bb2_cols, BB2_cols, B_het1_cols, B_het2_cols,
                    aa1_cols, AA1_cols, aa2_cols, AA2_cols, A_het1_cols, A_het2_cols,
                    ee1_cols, EE1_cols, ee2_cols, EE2_cols, E_het1_cols, E_het2_cols,
                    vect2, strand1_num, strand2_num,
                    LD_list)
    names(out_list) = c("freq_cols", "freq_cols1", "freq_cols2", "fit_cols",
                        "bb1_cols", "BB1_cols", "bb2_cols", "BB2_cols", "B_het1_cols", "B_het2_cols",
                        "aa1_cols", "AA1_cols", "aa2_cols", "AA2_cols", "A_het1_cols", "A_het2_cols",
                        "ee1_cols", "EE1_cols", "ee2_cols", "EE2_cols", "E_het1_cols", "E_het2_cols",
                        "vect2", "strand1_num", "strand2_num",
                        "LD_list")
    
    # merge the current index list with the other index lists
    out_list = c(out_list,index1,index2,index3)
    
    return(out_list)
    
}

#### Model functions ####

# this function selectively applies calc_equilibriumC()
# to the loci that are passed to it, allowing equilibrium for
# any combination of the loci to be calculated
multi_equilibrium = function(equilibrium_loci = equilibrium_loci) {
    
    if("B" %in% equilibrium_loci) {
        # populate the equilibrium matrix with calc_equilibriumC
        # see epigen_funcs1.4.cpp for details
        equilibrium_mat[,c("prop_B_1","prop_B_2","equil_B")] = 
          calc_equilibriumC(freq_mat,
                            index_list[["BB1_cols"]],
                            index_list[["BB2_cols"]],
                            index_list[["B_het1_cols"]],
                            index_list[["B_het2_cols"]],
                            equilibrium_mat[,c("prop_B_1","prop_B_2","equil_B"),drop = F],
                            equilibrium_threshold)
    }
    
    if("A" %in% equilibrium_loci) {
        equilibrium_mat[,c("prop_A_1","prop_A_2","equil_A")] = 
          calc_equilibriumC(freq_mat,
                            index_list[["AA1_cols"]],
                            index_list[["AA2_cols"]],
                            index_list[["A_het1_cols"]],
                            index_list[["A_het2_cols"]],
                            equilibrium_mat[,c("prop_A_1","prop_A_2","equil_A"),drop = F],
                            equilibrium_threshold)
    }
    
    if("E" %in% equilibrium_loci) {
        equilibrium_mat[,c("prop_E_1","prop_E_2","equil_E")] = 
          calc_equilibriumC(freq_mat,
                            index_list[["EE1_cols"]],
                            index_list[["EE2_cols"]],
                            index_list[["E_het1_cols"]],
                            index_list[["E_het2_cols"]],
                            equilibrium_mat[,c("prop_E_1","prop_E_2","equil_E"),drop = F],
                            equilibrium_threshold)
    }
    
    return(equilibrium_mat)
    
}

# record the frequencies of selected loci
multi_record_gen = function(record_loci = record_loci,
                            record_freqs_mat = record_freqs_mat,
                            equilibrium_mat = equilibrium_mat) {
    
    if("B" %in% record_loci) {
        # create a vector of NAs for all parameter combinations
        prop_B_1_vect = rep(NA,nrow(record_freqs_mat))
        # save equilibirum mat for rows that are still being simulated
        prop_B_1_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_B_1"]
        # append this vector to the recording matrix and name it based on the generaiton
        record_freqs_mat = cbind(record_freqs_mat,prop_B_1_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_B_1_gen_",cur_gen)
        # same for deme 2
        prop_B_2_vect = rep(NA,nrow(record_freqs_mat))
        prop_B_2_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_B_2"]
        record_freqs_mat = cbind(record_freqs_mat,prop_B_2_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_B_2_gen_",cur_gen)
    }
    
    if("A" %in% record_loci) {
        prop_A_1_vect = rep(NA,nrow(record_freqs_mat))
        prop_A_1_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_A_1"]
        record_freqs_mat = cbind(record_freqs_mat,prop_A_1_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_A_1_gen_",cur_gen)
        prop_A_2_vect = rep(NA,nrow(record_freqs_mat))
        prop_A_2_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_A_2"]
        record_freqs_mat = cbind(record_freqs_mat,prop_A_2_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_A_2_gen_",cur_gen)
    }
    
    if("E" %in% record_loci) {
        prop_E_1_vect = rep(NA,nrow(record_freqs_mat))
        prop_E_1_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_E_1"]
        record_freqs_mat = cbind(record_freqs_mat,prop_E_1_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_E_1_gen_",cur_gen)
        prop_E_2_vect = rep(NA,nrow(record_freqs_mat))
        prop_E_2_vect[!cumulative_output_rows] = equilibrium_mat[,"prop_E_2"]
        record_freqs_mat = cbind(record_freqs_mat,prop_E_2_vect)
        colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("prop_E_2_gen_",cur_gen)
    }
    
    return(record_freqs_mat)
    
}



# calulate LD between a pair of (epi)alleles
# can record LD, normalized LD and correlation
calc_LD = function(prop_list,
                   record_LD,
                   record_freqs_mat,
                   LD_list) {
    
    # loop through each LD metric which has been asked to be calculated
    for (i in record_LD) {
        
        measure = strsplit(i,"_")[[1]][1] # get the LD metric to calculate
        alleles = strsplit(i,"_")[[1]][2] # get alleles to calculate LD
        allele1 = substr(alleles,1,1)
        allele2 = substr(alleles,2,2)
        pop = strsplit(i,"_")[[1]][3] # get the deme to calculate LD in, 0 is both
        
        # use the prop_list to get the frequencies of each haplotype
        # selecting the population based on the record_LD string
        if (pop == 1) {
            haplo_freq_total = prop_list$haplo_freq_pop1
        } else if (pop == 2) {
            haplo_freq_total = prop_list$haplo_freq_pop2
        } else if (pop == 0) {
            haplo_freq_total = (prop_list$haplo_freq_pop1 + prop_list$haplo_freq_pop2) / 2
        } else {
            stop("Invalid population for calculating LD, most be one of: 1 2 0")
        }
        
        if (measure %in% c("D","Dnorm","cor")) {
            
            # calculate when alleles are together (p12)
            # and their total frequencies (p1 and p2)
            p12 = apply(haplo_freq_total[,LD_list[[alleles]], drop = F],1,sum)
            p1 = apply(haplo_freq_total[,LD_list[[allele1]], drop = F],1,sum)
            p2 = apply(haplo_freq_total[,LD_list[[allele2]], drop = F],1,sum)
            D = p12 - p1*p2 # calulate LD
            
            if (measure == "D") {
                
                # output LD for rows (parameter sets) that are still being simulated
                D_col = rep(NA,nrow(record_freqs_mat))
                D_col[!cumulative_output_rows] = D
                record_freqs_mat = cbind(record_freqs_mat,D_col)
                colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("D_",alleles,"_",pop,"_gen_",cur_gen)
                
            } else if (measure == "Dnorm") {
                
                # calculate the max possible LD pased on allele frequencies
                Dmax = ifelse(D > 0,
                              pmin(p1*(1-p2),(1-p1)*p2),
                              pmin(p1*p2,(1-p1)*(1-p2)))
                
                # calulate normalized LD, maintaining the sign of D
                Dnorm = D/Dmax
                
                # output LD for rows (parameter sets) that are still being simulated
                Dnorm_col = rep(NA,nrow(record_freqs_mat))
                Dnorm_col[!cumulative_output_rows] = Dnorm
                record_freqs_mat = cbind(record_freqs_mat,Dnorm_col)
                colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("Dnorm_",alleles,"_",pop,"_gen_",cur_gen)
                
            } else if (measure == "cor") {
                
                # calculate correlation
                cor12 = (D^2)/(p1*(1-p1)*p2*(1-p2))
                
                # output LD for rows (parameter sets) that are still being simulated
                cor_col = rep(NA,nrow(record_freqs_mat))
                cor_col[!cumulative_output_rows] = cor12
                record_freqs_mat = cbind(record_freqs_mat,cor_col)
                colnames(record_freqs_mat)[ncol(record_freqs_mat)] = paste0("cor_",alleles,"_",pop,"_gen_",cur_gen)
                
            } else {}
            
        } else {
            stop("Invalid LD measure provided, must be one of: D Dnorm cor")
        }
        
    }
    
    return(record_freqs_mat)
    
}
