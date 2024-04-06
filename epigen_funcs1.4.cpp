// epigen_funcs.cpp
#include <Rcpp.h>
using namespace Rcpp;

// v1.4 made injection_migrationC more genralized by allowing non-zero bb1_cols values
// converted e_X to c_X
// TO DO: implement LD calculation in C++
// TO DO: see if I can make everything faster using itterators

// migration between demes at rate m
// [[Rcpp::export]]                   
NumericMatrix migrationC(NumericMatrix freq_mat, NumericVector m) {

    // writing to a new output matrix is required since each genotype is
    // being written to and used to calculate the matching genotype in the
    // other population
    // NOTE: could just store genotypes for one pop. temporarily but this way is easier to read

    // make an output matrix and pass it the column names from the input matrix
    NumericMatrix freq_mat_out(freq_mat.nrow(),freq_mat.ncol());
    colnames(freq_mat_out) = colnames(freq_mat);

    // there are 36 (epi)genotypes within each deme, so we loop over each one
    for(int i = 0; i < 36; i++) {
        
        // 'i' selects from deme 1 and 'i+36' selects from deme 2
        // since (epi)genotypes for deme 1 come first in freq_mat
        freq_mat_out(_,i) = (1-m)*freq_mat(_,i) + m*freq_mat(_,i+36);
        freq_mat_out(_,i+36) = (1-m)*freq_mat(_,i+36) + m*freq_mat(_,i);

    }

    return freq_mat_out;

}

// migration between demes where migrants from deme 1 to deme 2 are converted 
// from bb to BB at the neutral locus
// Note: this only works when both demes are fixed for bb
// [[Rcpp::export]]
NumericMatrix injection_migrationC(NumericMatrix freq_mat, 
                                   NumericVector m,
                                   NumericVector bb1_cols,
                                   NumericVector BB2_cols,
                                   NumericVector bb2_cols) {

    // looping over each column which contains a bb/BB genotype
    for(int i = 0; i < bb1_cols.length(); i++) {
        for(int j = 0; j < freq_mat.nrow(); j++) {
            // convert bb migrants from deme 1 to BB in deme 2
            freq_mat(j,BB2_cols[i]) = m[j] * freq_mat(j,bb1_cols[i]); 
            // normal migration from deme 2 to deme 1
            freq_mat(j,bb1_cols[i]) = (1-m[j]) * freq_mat(j,bb1_cols[i]) + m[j] * freq_mat(j,bb2_cols[i]);
            // loss of migrant bb individuals in deme 2
            freq_mat(j,bb2_cols[i]) = (1-m[j]) * freq_mat(j,bb2_cols[i]);
        }
    }

    return freq_mat;

}

// [[Rcpp::export]]
NumericMatrix selectionC(NumericMatrix freq_mat, NumericMatrix selection_mat) {

    // multiply genotype frequencies by their fitness
    for(int i = 0; i < 72; i++) {
        freq_mat(_,i) = freq_mat(_,i) * selection_mat(_,i);
    }

    double sum_pop1 = 0;
    double sum_pop2 = 0;
    
    for(int j = 0;j < freq_mat.nrow(); j++) {
        // calculating the total frequency (mean fitness) for each deme
        // for some reason you have to use Range(j,j) instead of j
        sum_pop1 =  std::accumulate(freq_mat(Range(j,j),Range(0,35)).begin(),
                                    freq_mat(Range(j,j),Range(0,35)).end(),
                                    0.0); // 0.0 to set the type to double
        sum_pop2 =  std::accumulate(freq_mat(Range(j,j),Range(36,71)).begin(),
                                    freq_mat(Range(j,j),Range(36,71)).end(),
                                    0.0);
        // devide each genotype by mean fitness.
        // to return the total genotype freqency to 1
        for(int k = 0; k < 36; k++) {
            freq_mat(j,k) /= sum_pop1;
            freq_mat(j,k+36) /= sum_pop2;
        }
    }

    return freq_mat;

}

// [[Rcpp::export]]
NumericMatrix recombinationC(NumericMatrix freq_mat,
                             NumericVector no_reco,
                             NumericVector reco,
                             NumericVector just_r1,
                             NumericVector just_r2,
                             NumericVector both_r,
                             NumericVector r1,
                             NumericVector r2) {

    // an output matrix is required since genotypes are being written to and used for calculations
    NumericMatrix freq_mat_out(freq_mat.nrow(),freq_mat.ncol());
    colnames(freq_mat_out) = colnames(freq_mat);

    // shift (epi)genotypes from their starting location i.e., reco to their end location e.g., just_r1, just_r2 or r2_r1
    // amount of transfer is based on recombination rate e.g., reco moves to just_r1 at rate r1 - r1*r2

    double r1_freqs = 0;
    double r2_freqs = 0;
    double both_r_freqs = 0;

    for(int k = 0; k < freq_mat.nrow(); k++) {
        // reco, just_r1, just_r2 and both_r have the same lengths
        // to we itterate through all of them in the same loop
        for(int i = 0; i < reco.length(); i++) {

            // calculate frequencies of each recombination event based on r1 and r2
            // for (epi)genotypes that are affected by recombination
            r1_freqs = r1[k] * freq_mat(k,reco[i]);
            r2_freqs = r2[k] * freq_mat(k,reco[i]);
            both_r_freqs = r1[k] * r2[k] * freq_mat(k,reco[i]);

            // add the frequencies of recombination events to the new (epi)genotypes after recombination
            freq_mat_out(k,just_r1[i]) = freq_mat_out(k,just_r1[i]) + r1_freqs - both_r_freqs;
            freq_mat_out(k,just_r2[i]) = freq_mat_out(k,just_r2[i]) + r2_freqs - both_r_freqs;
            freq_mat_out(k,both_r[i]) = freq_mat_out(k,both_r[i]) + both_r_freqs;
            // subtract the frequencies of recombination events from the (epi)genotypes prior to recombination
            freq_mat_out(k,reco[i]) = freq_mat_out(k,reco[i]) + freq_mat(k,reco[i]) - r1_freqs - r2_freqs + both_r_freqs;

        }

        // add (epi)genotypes that are not affected by recombination to the output matrix
        for(int j = 0; j < no_reco.length(); j++) {
            freq_mat_out(k,no_reco[j]) = freq_mat(k,no_reco[j]);
        }

    }

    return freq_mat_out;

}   

// [[Rcpp::export]]          
List calc_propsC (NumericMatrix freq_mat,
                  NumericVector strand1,
                  NumericVector strand2) {

    // create matricies to store haplotype freqencies                    
    NumericMatrix haplo_freq_pop1(freq_mat.nrow(),8);
    NumericMatrix haplo_freq_pop2(freq_mat.nrow(),8);

    // getting the cumulative frequency of haplotypes from genotypes e.g. aA_eE_bB has aeb on strand1 and AEB on strand2
    // there are 8 unique haplotypes from the 36 unique genotypes in each pop.

    for(int i = 0; i < 36; i++) {
        // Rcout << strand1[i] << "\n";
        // Rcout << strand2[i] << "\n";
        // Rcout << freq_cols[i] << "\n";
        haplo_freq_pop1(_,strand1[i]) = haplo_freq_pop1(_,strand1[i]) + freq_mat(_,i);
        haplo_freq_pop1(_,strand2[i]) = haplo_freq_pop1(_,strand2[i]) + freq_mat(_,i);
        haplo_freq_pop2(_,strand1[i]) = haplo_freq_pop2(_,strand1[i]) + freq_mat(_,i+36);
        haplo_freq_pop2(_,strand2[i]) = haplo_freq_pop2(_,strand2[i]) + freq_mat(_,i+36);
    }
    
    // we must divide total haplotype frequencies by two since we are suming genotype frequencies twice
    // i.e. once for the first strand and once for the second strand
    // This is done after the totals are calculated since when both strands are the same e.g., AA_EE_BB
    // this genotype contributes doubly to the AEB haplotype

    for(int j = 0; j < haplo_freq_pop1.ncol(); j++) {
        haplo_freq_pop1(_,j) = haplo_freq_pop1(_,j) / 2;
        haplo_freq_pop2(_,j) = haplo_freq_pop2(_,j) / 2;
    }
    
    // combine the haplotype matricies into a named list, so that it can be passed back to R
    List prop_list = List::create(Named("haplo_freq_pop1") = haplo_freq_pop1,
                            Named("haplo_freq_pop2") = haplo_freq_pop2);
                            
    
    return prop_list;
  
} 

// [[Rcpp::export]]
NumericMatrix matingC (List prop_list, NumericVector VectOf2s, CharacterVector freq_cols) {

    NumericMatrix pop1_props = prop_list[0];
    NumericMatrix pop2_props = prop_list[1];
    
    // here we are generating two matricies whos column correspond to each possible combination of haplotypes
    // e.g. the first matrix eggs1 will have the first haplotype aeb repeated 8 times for the first 8 columns
    // and the second matrix sprm1 will have haplotype 1 through 8 in the first 8 columns
    // then we fill the next 7 columns with the second haplotype in eggs1 and haplotypes 2 through 8 in sprm1, etc.

    int nrows = pop1_props.nrow();
    int ncols = 36;
    int cumu_i = 0;
    
    NumericMatrix eggs1(nrows,ncols);
    NumericMatrix sprm1(nrows,ncols);
    NumericMatrix eggs2(nrows,ncols);
    NumericMatrix sprm2(nrows,ncols);
    NumericMatrix freq_mat(nrows,2*ncols);
    colnames(freq_mat) = freq_cols; 

    for(int i = 0; i < 8; i++) {
        
        // Rcout << "cumu_i = " << cumu_i << "\n";
        // Rcout << "i = " << i << "\n";
        // Rcout << "i+j = ";
        // Rcout << "cumu_i + j = ";

        for(int j = 0;j < 8-i;j++) {

            eggs1(_,cumu_i+j) = pop1_props(_,i);
            sprm1(_,cumu_i+j) = pop1_props(_,j+i);
            eggs2(_,cumu_i+j) = pop2_props(_,i);
            sprm2(_,cumu_i+j) = pop2_props(_,j+i);

            // Rcout << i+j << " ";
            // Rcout << cumu_i + j << " ";

        }

        // Rcout << "\n";
        
        cumu_i = cumu_i + 8 - i;
        
    }

    // we then multiply out each unique haplotype combination to get our new genotype frequencies
    // the multiplication by VectOf2s is to account for the fact that when the haplotypes are different
    // they are form in two ways e.g. Aa_Ee_Bb_1 can form with eggs1 = aeb and sprm1 = AEB
    // or eggs1 = AEB and sprm1 = aeb

    for(int k = 0; k < nrows; k++) {
        for(int w = 0; w < 36; w++) {
            freq_mat(k,w) = eggs1(k,w) * sprm1(k,w) * VectOf2s[w];
            freq_mat(k,w+36) = eggs2(k,w) * sprm2(k,w) * VectOf2s[w];
        }
        
    } 

    return freq_mat;

}


// [[Rcpp::export]]
NumericMatrix epimutationC (NumericMatrix freq_mat,
                            NumericMatrix c_mat,
                            NumericMatrix i_mat,
                            NumericVector ee_cols,
                            NumericVector eE_cols,
                            NumericVector Ee_cols,
                            NumericVector EE_cols,
                            NumericVector ee_to_eE_cols,
                            NumericVector ee_to_Ee_cols,
                            NumericVector ee_to_EE_cols,
                            NumericVector eE_to_ee_cols,
                            NumericVector eE_to_Ee_cols,
                            NumericVector eE_to_EE_cols,
                            NumericVector Ee_to_ee_cols,
                            NumericVector Ee_to_eE_cols,
                            NumericVector Ee_to_EE_cols,
                            NumericVector EE_to_ee_cols,
                            NumericVector EE_to_eE_cols,
                            NumericVector EE_to_Ee_cols) {

    // creating output matrix since we are writing to columns that we are using for calculation
    NumericMatrix freq_mat_out(freq_mat.nrow(),freq_mat.ncol());
    colnames(freq_mat_out) = colnames(freq_mat);
    
    // XX_cols are all columns which have that epigenotype
    // XX_to_YY_cols are columns which have their frequency increase if an XX to YY transition happens

    // for some reason when you try and do this using the vectorized notation it does not work
    // i.e. freq_mat(_,i) instead of freq_mat(x,i)

    for(int x = 0; x < freq_mat.nrow(); x++) {

        // sub-loops are based on what input columns output columns are dependent on
        // e.g., all ee_to_YY_cols are the same length
        for(int i = 0; i < ee_cols.length(); i++) {        
            freq_mat_out(x,ee_cols[i])       += (1-i_mat(x,ee_cols[i])) * (1-i_mat(x,ee_cols[i])) * freq_mat(x,ee_cols[i]);
            freq_mat_out(x,ee_to_eE_cols[i]) += (1-i_mat(x,ee_cols[i])) * i_mat(x,ee_cols[i])     * freq_mat(x,ee_cols[i]);
            freq_mat_out(x,ee_to_Ee_cols[i]) += i_mat(x,ee_cols[i])     * (1-i_mat(x,ee_cols[i])) * freq_mat(x,ee_cols[i]);
            freq_mat_out(x,ee_to_EE_cols[i]) += i_mat(x,ee_cols[i])     * i_mat(x,ee_cols[i])     * freq_mat(x,ee_cols[i]);
        }

        for(int j = 0; j < eE_cols.length(); j++) {        
            freq_mat_out(x,eE_to_ee_cols[j]) += (1-i_mat(x,eE_cols[j])) * c_mat(x,eE_cols[j])     * freq_mat(x,eE_cols[j]);
            freq_mat_out(x,eE_cols[j])       += (1-i_mat(x,eE_cols[j])) * (1-c_mat(x,eE_cols[j])) * freq_mat(x,eE_cols[j]);
            freq_mat_out(x,eE_to_Ee_cols[j]) += i_mat(x,eE_cols[j])     * c_mat(x,eE_cols[j])     * freq_mat(x,eE_cols[j]);
            freq_mat_out(x,eE_to_EE_cols[j]) += i_mat(x,eE_cols[j])     * (1-c_mat(x,eE_cols[j])) * freq_mat(x,eE_cols[j]);
        }

        for(int k = 0; k < Ee_cols.length(); k++) {        
            freq_mat_out(x,Ee_to_ee_cols[k]) += c_mat(x,Ee_cols[k])     * (1-i_mat(x,Ee_cols[k])) * freq_mat(x,Ee_cols[k]);
            freq_mat_out(x,Ee_to_eE_cols[k]) += c_mat(x,Ee_cols[k])     * i_mat(x,Ee_cols[k])     * freq_mat(x,Ee_cols[k]);
            freq_mat_out(x,Ee_cols[k])       += (1-c_mat(x,Ee_cols[k])) * (1-i_mat(x,Ee_cols[k])) * freq_mat(x,Ee_cols[k]);
            freq_mat_out(x,Ee_to_EE_cols[k]) += (1-c_mat(x,Ee_cols[k])) * i_mat(x,Ee_cols[k])     * freq_mat(x,Ee_cols[k]);
        }

        for(int w = 0; w < EE_cols.length(); w++) {        
            freq_mat_out(x,EE_to_ee_cols[w]) += c_mat(x,EE_cols[w])     * c_mat(x,EE_cols[w])     * freq_mat(x,EE_cols[w]);
            freq_mat_out(x,EE_to_eE_cols[w]) += c_mat(x,EE_cols[w])     * (1-c_mat(x,EE_cols[w])) * freq_mat(x,EE_cols[w]);
            freq_mat_out(x,EE_to_Ee_cols[w]) += (1-c_mat(x,EE_cols[w])) * c_mat(x,EE_cols[w])     * freq_mat(x,EE_cols[w]);
            freq_mat_out(x,EE_cols[w])       += (1-c_mat(x,EE_cols[w])) * (1-c_mat(x,EE_cols[w])) * freq_mat(x,EE_cols[w]);
        }

    }

    // Rcout << "lastest epimutation" << "\n";
    
    return freq_mat_out;

}

// [[Rcpp::export]]
NumericVector RI_calcC (NumericMatrix out_mat,
                        NumericVector m,
                        NumericVector BB1_cols,
                        NumericVector BB2_cols,
                        NumericVector B_het1_cols,
                        NumericVector B_het2_cols) {

    // have to give m as a seperate vector since the position of m within
    // out_mat changes depending on which variable(s) you parallelize over

    // set up our variables used to calulate RI
    double prop_B_1 = 0;
    double prop_B_2 = 0;
    NumericMatrix RI_mat(out_mat.nrow(),2);
    colnames(RI_mat) = CharacterVector::create("RI","RI_2");

    for(int k = 0; k < out_mat.nrow(); k++) {

        // reset prop_B_X to zero for this row (parameter set)
        prop_B_1 = 0;
        prop_B_2 = 0;

        // add the frequency of BB individuals to the proportion of B in each deme
        for(int i = 0; i < BB1_cols.length(); i++) {
            prop_B_1 += out_mat(k,BB1_cols[i]);
            prop_B_2 += out_mat(k,BB2_cols[i]);
        }

        // add 0.5 times the frequency of Bb individuals to the proportion of B in each deme
        for(int j = 0; j < B_het1_cols.length(); j++) {
            prop_B_1 += 0.5*out_mat(k,B_het1_cols[j]);
            prop_B_2 += 0.5*out_mat(k,B_het2_cols[j]);
        }

        // calulate RI based on the average proportion of B between demes
        RI_mat(k,0) = 1 - ((0.5 * (prop_B_1 + prop_B_2)) / (0.5 * m[k]));
        // calculate RI based on the proportion of B in deme two
        RI_mat(k,1) = 1 - (prop_B_2 / (0.5 * m[k]));
        // if equilibrium is reached (and there is post-contact migration)
        // these values should be the same

    }

    return RI_mat;   

} 


// [[Rcpp::export]]
NumericMatrix calc_equilibriumC (NumericMatrix freq_mat,
                                 NumericVector homo_cols1,
                                 NumericVector homo_cols2,
                                 NumericVector hetero_cols1,
                                 NumericVector hetero_cols2,
                                 NumericMatrix old_prop_mat,
                                 double equilibrium_threshold) {

    // create a matrix to store frequencies of (epi)alleles
    NumericMatrix prop_mat(freq_mat.nrow(),3);
    colnames(prop_mat) = CharacterVector::create("prop1","prop2","equilibrium_reached");

    // create doubles to store proportions
    double prop1 = 0;
    double prop1_dif = 0;
    double prop2 = 0;
    double prop2_dif = 0;

    for(int k = 0; k < freq_mat.nrow(); k++) {

        // reset proportions to zero for this row (parameter set)
        prop1 = 0;
        prop2 = 0;

        // add the frequency of homozygous (epi)genotypes to the proportion
        for(int i = 0; i < homo_cols1.length(); i++) {
            prop1 += freq_mat(k,homo_cols1[i]);
            prop2 += freq_mat(k,homo_cols2[i]);
        }

        // add 0.5 times the frequency of heterozygous (epi)genotypes to the proportion
        for(int j = 0; j < hetero_cols1.length(); j++) {
            prop1 += 0.5*freq_mat(k,hetero_cols1[j]);
            prop2 += 0.5*freq_mat(k,hetero_cols2[j]);
        }

        // caluclate the change in proportion by subtracting the proportion
        // of the previous generation
        prop1_dif = old_prop_mat(k,0) - prop1;
        if(prop1_dif < 0) {
            prop1_dif = 0 - prop1_dif;
        }
        prop2_dif = old_prop_mat(k,1) - prop2;
        if(prop2_dif < 0) {
            prop2_dif = 0 - prop2_dif;
        }

        // set the equilibrium status to 1 if the change in frequency
        // is less than equilibrium_threshold
        if((prop1_dif < equilibrium_threshold) & (prop2_dif < equilibrium_threshold)) {
            prop_mat(k,0) = prop1;
            prop_mat(k,1) = prop2;
            prop_mat(k,2) = 1;
        } else {
            prop_mat(k,0) = prop1;
            prop_mat(k,1) = prop2;
        }

    }
    
    return prop_mat;   

} 
