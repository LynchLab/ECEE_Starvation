cd /Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_analysis
python calculate_well_mixed_hmm_wrapper.py > calculate_well_mixed_hmm_wrapper.o 2> calculate_well_mixed_hmm_wrapper.e
python calculate_clade_hmm_wrapper.py > calculate_clade_hmm_wrapper.o 2> calculate_clade_hmm_wrapper.e

cd /Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary
perl basal_fixation_probability_by_mutate_TP.pl
perl count_mutations_fixed_any.pl
perl count_mutations_fixed_basal_major.pl
perl list_basal_fixed_alleles.pl
perl list_clade_fixed_alleles.pl
perl list_basal_poly_alleles.pl
perl summary_by_mutate_TP.pl
perl summary_temporal.pl
perl summary.pl
R CMD BATCH calculate_longest_coexist_length.R
R CMD BATCH calculate_num_fixed_basal_byTP.R
R CMD BATCH select_gene_mutation_rate_related.R