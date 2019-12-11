# usr/bin/perl

open OUT1, ">/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_major_clade_fixed_alleles.txt";
open OUT2, ">/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_minor_clade_fixed_alleles.txt";
open OUT3, ">/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_major_minor_clade_fixed_alleles.txt";

$daf_cutoff_lower = 0.05;
$daf_cutoff_upper = 0.95;

my @pop = (101..124, 201..224, 301..324, 401..424, 501..524); #my @pop = (101);
my @wt = (1,12,13,24);
my @del_araBAD = (3,8,15,20);
my @del_mutL = (4,9,16,21);
my @del_mutL_araBAD = (6,11,18,23);
my @del_srlD = (5,10,17,22);
my @del_mutL_srlD = (2,7,14,19);

my %pop2pop_type;
foreach my $p(@wt){
	$pop2pop_type{$p} = "wt";
}
foreach my $p(@del_araBAD){
	$pop2pop_type{$p} = "del_araBAD";
}
foreach my $p(@del_mutL){
	$pop2pop_type{$p} = "del_mutL";
}
foreach my $p(@del_mutL_araBAD){
	$pop2pop_type{$p} = "del_mutL_araBAD";
}
foreach my $p(@del_srlD){
	$pop2pop_type{$p} = "del_srlD";
}
foreach my $p(@del_mutL_srlD){
	$pop2pop_type{$p} = "del_mutL_srlD";
}

# read the available time points
my %numT;
open IN, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/all_annotated_timecourse_available_timepoints.txt";

while (my $line = <IN>){
	chomp $line;
	$line =~ /(\d+)\s(\d)/;
	$numT{$1} = $2;
}
close IN;

# see the running results of clade-aware HHM 
my %cladeHMMResult;
open IN2, "/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_analysis/calculate_clade_hmm_wrapper.o";

while (my $line = <IN2>){
	chomp $line;
	if ($line =~ /(\d+): length of coexistence = (\d+)/){ #101: length of coexistence = 990
		if ($2 == 0){
			$cladeHMMResult{$1} = 0;  # "Insufficient evidence for frequency-dependence!"
		}else{
			$cladeHMMResult{$1} = 1;  # clade-aware HMM works well
		}
	}
	if ($line =~ /error:(\d+)/){ #error:105
		$cladeHMMResult{$1} = -1; # HMM was not running b/c limited number of time points or sum of good_majority_idxs=0
	}
}
close IN2;

print OUT1 "pop\tpop_type\tfixedTP\tmut_type\tpos\tgene_name\tmut_type2\tsel_coef\tsel_coef_TP1\tsel_coef_TP2\n";
print OUT2 "pop\tpop_type\tfixedTP\tmut_type\tpos\tgene_name\tmut_type2\tsel_coef\tsel_coef_TP1\tsel_coef_TP2\n";
print OUT3 "pop\tpop_type\tfixedTP\tmut_type\tpos\tgene_name\tmut_type2\tsel_coef\tsel_coef_TP1\tsel_coef_TP2\n";


foreach my $p(@pop){
	print "checking: ".$p."\n";

	if ($cladeHMMResult{$p} !=-1){

		open HMM2, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_haplotype_timecourse.txt";
		open DATA, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_annotated_timecourse.txt";
		open DAF, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/daf_time_matrix/".$p."_daf_time_matrix.txt";
	
		$varLine = <DATA>;
		$coverageLine = <DATA>;

		my $timeLine = <HMM2>;
		chomp($timeLine);
		my @timePoints = split /, /, $timeLine; 

		my $f_majorLine = <HMM2>;
		chomp($f_majorLine);
		my @f_major = split /, /, $f_majorLine; 
		my $f_minorLine = <HMM2>;
		chomp($f_minorLine);
		my @f_minor = split /, /, $f_minorLine; 


		foreach my $i(1..(8-2)){
			my $notUsedLine = <HMM2>;   #fmajors, fminors, hard_ns[0,:], hard_ns[5:,:].sum(axis=0), 
									#hard_ns[1,:], hard_ns[2,:], hard_ns[3,:], hard_ns[4,:]
									#clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3,'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
		}

		$varLine3 = <DAF>;

		while (my $line1 = <DATA>){

			my $line2 = <HMM2>;
			my $line3 = <DAF>;

			chomp($line1);
			chomp($line2);
			chomp($line3);

			my @temp1 = split /, /, $line1; #4174, thrC, T->G, synonymous,
			my $SNPType = $temp1[3];
			my $pos = $temp1[0];
			my $geneName = $temp1[1];
			$geneName =~ s/[\[\]]//g;
			my $type = $temp1[2];
			#print $SNPType."\n";

			my @temp2 = split /, /, $line2; 
			my $numTP = scalar(@temp2);
			#print $temp2[1]."\n";

			my @temp3 = split /, /, $line3;
			my @daf = @temp3[4..(scalar(@temp3)-1)]; 
			my @correctedDaf = (0);
			foreach my $i(1..($numTP-1)){
				if($temp2[$i] == 3 || $temp2[$i] == 6){
					push @correctedDaf, ($daf[$i] / $f_major[$i]);
				}elsif($temp2[$i] == 4 || $temp2[$i] == 7){
					push @correctedDaf, ($daf[$i] / $f_minor[$i]);
				}else{
					push @correctedDaf, $daf[$i];
				}
			}

			#print $daf[1]."\n";

			my $SNPstate = $temp2[$numTP-1];
			#print $SNPstate."\n";

			if ($SNPstate == 3){

				my $fixedTP = 0;
				foreach my $t(0..scalar(@temp2)-1){
					#print $temp2[$t]."\n";
					if($temp2[$t] == 3){
						$fixedTP = $timePoints[$t];
						#print $timePoints[$t]."\n";
						last;
					}
				}

				my $max_s = -1;
				my $max_d1 = -1;
				my $max_d2 = -1;

				foreach my $t(1..(scalar(@daf)-2)){
					if($daf[$t]>=$daf_cutoff_lower && $daf[$t+1]>=$daf_cutoff_lower && $daf[$t]<=$daf_cutoff_upper && $daf[$t+1]<=$daf_cutoff_upper){
						if($correctedDaf[$t]>=$daf_cutoff_lower && $correctedDaf[$t+1]>=$daf_cutoff_lower && $correctedDaf[$t]<=$daf_cutoff_upper && $correctedDaf[$t+1]<=$daf_cutoff_upper){
							my $s = log($correctedDaf[$t+1]*(1-$correctedDaf[$t])/$correctedDaf[$t]/(1-$correctedDaf[$t+1]))/($timePoints[$t+1]-$timePoints[$t]);
							if($s>$max_s && $s>0){
								$max_s = $s;
								$max_d1 = $timePoints[$t];
								$max_d2 = $timePoints[$t+1];
							}
						}
					}
				}

				if((!($type =~ /\./)) && ($SNPType ne "")){
					print OUT1 $p."\t".$pop2pop_type{($p%100)}."\t".$fixedTP."\t".$SNPType."\t".$pos."\t".$geneName."\t".$type."\t".$max_s."\t".$max_d1."\t".$max_d2."\n";
					print OUT3 $p."\t".$pop2pop_type{($p%100)}."\t".$fixedTP."\t".$SNPType."\t".$pos."\t".$geneName."\t".$type."\t".$max_s."\t".$max_d1."\t".$max_d2."\n";
				}
			}

			if ($SNPstate == 4){

				my $fixedTP = 0;
				foreach my $t(0..scalar(@temp2)-1){
					#print $temp2[$t]."\n";
					if($temp2[$t] == 4){
						$fixedTP = $timePoints[$t];
						#print $timePoints[$t]."\n";
						last;
					}
				}

				my $max_s = -1;
				my $max_d1 = -1;
				my $max_d2 = -1;

				foreach my $t(1..(scalar(@daf)-2)){
					if($daf[$t]>=$daf_cutoff_lower && $daf[$t+1]>=$daf_cutoff_lower && $daf[$t]<=$daf_cutoff_upper && $daf[$t+1]<=$daf_cutoff_upper){
						if($correctedDaf[$t]>=$daf_cutoff_lower && $correctedDaf[$t+1]>=$daf_cutoff_lower && $correctedDaf[$t]<=$daf_cutoff_upper && $correctedDaf[$t+1]<=$daf_cutoff_upper){
							my $s = log($correctedDaf[$t+1]*(1-$correctedDaf[$t])/$correctedDaf[$t]/(1-$correctedDaf[$t+1]))/($timePoints[$t+1]-$timePoints[$t]);
							if($s>$max_s && $s>0){
								$max_s = $s;
								$max_d1 = $timePoints[$t];
								$max_d2 = $timePoints[$t+1];
							}
						}
					}
				}
				
				if((!($type =~ /\./)) && ($SNPType ne "")){
					print OUT2 $p."\t".$pop2pop_type{($p%100)}."\t".$fixedTP."\t".$SNPType."\t".$pos."\t".$geneName."\t".$type."\t".$max_s."\t".$max_d1."\t".$max_d2."\n";
					print OUT3 $p."\t".$pop2pop_type{($p%100)}."\t".$fixedTP."\t".$SNPType."\t".$pos."\t".$geneName."\t".$type."\t".$max_s."\t".$max_d1."\t".$max_d2."\n";				
				}
			}
		}
		close DATA;
		close HMM2;
		close SELCOEF;
	}
}