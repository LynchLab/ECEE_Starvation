# usr/bin/perl

open OUT, ">/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_basal_poly_alleles.txt";

$daf_cutoff_lower = 0.05;
$daf_cutoff_upper = 0.95;

my @pop = (101..124, 201..224, 301..324, 401..424, 501..524); #
#my @pop = (101);
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

print OUT "pop\tpop_type\tfixedTP\tmut_type\tpos\tgene_name\tmut_type2\tsel_coef\tsel_coef_TP1\tsel_coef_TP2\n";

foreach my $p(@pop){
	print "checking: ".$p."\n";

	my @timePoints = ();

	open DATA, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_annotated_timecourse.txt";
	open DAF, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/daf_time_matrix/".$p."_daf_time_matrix.txt";

	if ($cladeHMMResult{$p} !=-1){
		open HMM, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_haplotype_timecourse.txt";
		my $timeLine = <HMM>;
		chomp($timeLine);
		@timePoints = split /, /, $timeLine;
		foreach my $i(1..8){
			my $notUsedLine = <HMM>;   #fmajors, fminors, hard_ns[0,:], hard_ns[5:,:].sum(axis=0), 
									#hard_ns[1,:], hard_ns[2,:], hard_ns[3,:], hard_ns[4,:]
									#clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3,'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
		}
	}else{
		open HMM, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_well_mixed_state_timecourse.txt";
		my $timeLine = <HMM>;
		chomp($timeLine);
		@timePoints = split /, /, $timeLine; 

		foreach my $i(1..4){
			my $notUsedLine = <HMM>;   #hard_ns[0,:], hard_ns[1,:], hard_ns[2,:], hard_ns[3,:]
									#well_mixed_hmm_states = {'A':0,'E':1,'F':2,'P':3}    
		}
	}	
	
	$varLine1 = <DATA>;
	$coverageLine = <DATA>;

	$varline3 = <DAF>;

	while (my $line1 = <DATA>){

		my $line2 = <HMM>;
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
		#print $type."\n";

		my @temp2 = split /, /, $line2; 
		my $numTP = scalar(@temp2);
		#print $temp2[1]."\n";

		my @temp3 = split /, /, $line3; #4174, thrC, T->G, synonymous,
		my @daf = @temp3[4..(scalar(@temp3)-1)]; 
		#print $daf[1]."\n";

		my $num_daf_50 = 0;
		foreach my $f(@daf){
			if($f >= 0.5){
				$num_daf_50 ++;
			}
		}

		my $SNPstate = $temp2[$numTP-1];
		#print $SNPstate."\n";

		if (($cladeHMMResult{$p}!=-1 && $SNPstate == 5 && $num_daf_50>=2) || ($cladeHMMResult{$p}==-1 && $SNPstate == 3 && $num_daf_50>=2) ){

			my $fixedTP = -1;
			
			my $max_s = -1;
			my $max_d1 = -1;
			my $max_d2 = -1;

			foreach my $t(1..(scalar(@daf)-2)){
				if($daf[$t]>=$daf_cutoff_lower && $daf[$t+1]>=$daf_cutoff_lower && $daf[$t]<=$daf_cutoff_upper && $daf[$t+1]<=$daf_cutoff_upper){
					my $s = log($daf[$t+1]*(1-$daf[$t])/$daf[$t]/(1-$daf[$t+1]))/($timePoints[$t+1]-$timePoints[$t]);
					if($s>$max_s && $s>0){
						$max_s = $s;
						$max_d1 = $timePoints[$t];
						$max_d2 = $timePoints[$t+1];
					}
				}
			}

			if(!($type =~ /\./)){
				print OUT $p."\t".$pop2pop_type{($p%100)}."\t".$fixedTP."\t".$SNPType."\t".$pos."\t".$geneName."\t".$type."\t".$max_s."\t".$max_d1."\t".$max_d2."\n";
			}
		}
	}
	close DATA;
	close HMM;
	close DAF;
}