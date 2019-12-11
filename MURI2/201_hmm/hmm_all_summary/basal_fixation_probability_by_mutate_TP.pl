# usr/bin/perl

open OUT, ">/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/basal_fixation_probability_by_mutate_TP.txt";

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
$numT{"208"} = 0;
$numT{"224"} = 0;

# see the running results of clade-aware HHM 
my %cladeHMMResult;
my %coexistLength;
open IN2, "/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_analysis/calculate_clade_hmm_wrapper.o";

while (my $line = <IN2>){
	chomp $line;
	if ($line =~ /(\d+): length of coexistence = (\d+)/){ #101: length of coexistence = 990
		if ($2 == 0){
			$cladeHMMResult{$1} = 0;  # "Insufficient evidence for frequency-dependence!"
			$coexistLength{$1} = 0;
		}else{
			$cladeHMMResult{$1} = 1;  # clade-aware HMM works well
			$coexistLength{$1} = $2;
		}
	}
	if ($line =~ /error:(\d+)/){ #error:105
		$cladeHMMResult{$1} = -1; # HMM was not running b/c limited number of time points or sum of good_majority_idxs=0
		$coexistLength{$1} = 0;
	}
}
close IN2;
$cladeHMMResult{"208"} = -1;
$coexistLength{"208"} = 0;
$cladeHMMResult{"224"} = -1;
$coexistLength{"224"} = 0;

# count the status of SNPs in the most recent available time point
my %count;
my @SNPTypeArray = ("S", "N", "I");
my @SNPStateArray = ("E", "F", "P");
my %pop2TimePointArray;

foreach my $p(@pop){
	print "checking: ".$p."\n";

	if($cladeHMMResult{$p}==-1){
	
		open DATA, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_annotated_timecourse.txt";
		open HMM1, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_well_mixed_state_timecourse.txt";
	
		$varLine = <DATA>;
		$CoverageLine = <DATA>;

		my $timePointsStr = <HMM1>;
		chomp($timePointsStr);
		my @timePointArray = split /, /, $timePointsStr;
		$pop2TimePointArray{$p} = [@timePointArray];

		foreach my $i(1..4){
			my $notUsedLine = <HMM1>;   #hard_ns[0,:], hard_ns[1,:], hard_ns[2,:], hard_ns[3,:]
									#well_mixed_hmm_states = {'A':0,'E':1,'F':2,'P':3}    
		}

		while (my $line1 = <DATA>){


			my $line2 = <HMM1>;

			chomp($line1);
			chomp($line2);

			my @temp1 = split /, /, $line1;
			my $SNPType = $temp1[3];
			#print $SNPType."\n";

			my @temp2 = split /, /, $line2; 
			my $numTP = scalar(@temp2);
			#print $temp2[1]."\n";		

			my $SNPstate = "";
			if($temp2[$numTP-1] == 1){
				$SNPstate = "E";
			}elsif($temp2[$numTP-1] == 2){
				$SNPstate = "F";
			}elsif($temp2[$numTP-1] == 3){
				$SNPstate = "P";
			}

			#find mutate time point
			my $mutateTP = 0;
			foreach my $t_i(0..($numTP-1)){
				#if($temp2[$t_i]>0 ){
				if($temp2[$t_i]>0 && $temp2[$t_i+1]>1){  #at least exist at two consecutive points
					$mutateTP = $timePointArray[$t_i];
					last;
				}
			}
		
			$count{$p}{$mutateTP}{"all"} ++;

			if($temp2[$numTP-1] == 2){
				$count{$p}{$mutateTP}{"fixed_basal"} ++;
			}
		}
		close DATA;
		close HMM1;

	}else{

		open HMM2, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_haplotype_timecourse.txt";
		open DATA, "/Users/wei-chinho/Documents/MURI2/00_data/server/weiho/annotated_timecourse/".$p."_annotated_timecourse.txt";
	
		$varLine = <DATA>;

		my $timePointsStr = <HMM2>;
		chomp($timePointsStr);
		my @timePointArray = split /, /, $timePointsStr;
		$pop2TimePointArray{$p} = [@timePointArray];
		
		foreach my $i(1..8){
			my $notUsedLine = <HMM2>;   #fmajors, fminors, hard_ns[0,:], hard_ns[5:,:].sum(axis=0), 
									#hard_ns[1,:], hard_ns[2,:], hard_ns[3,:], hard_ns[4,:]
									#clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3,'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
		}

		while (my $line1 = <DATA>){

			my $line2 = <HMM2>;
			chomp($line1);
			chomp($line2);

			my @temp1 = split /, /, $line1;
			my $SNPType = $temp1[3];
			#print $SNPType."\n";

			my @temp2 = split /, /, $line2; 
			my $numTP = scalar(@temp2);
			#print $temp2[1]."\n";

			#find mutate time point
			my $mutateTP = 0;
			foreach my $t_i(0..($numTP-1)){
				#if($temp2[$t_i]>0 ){
				if($temp2[$t_i]>0 && $temp2[$t_i+1]>1){  #at least exist at two consecutive points
					$mutateTP = $timePointArray[$t_i];
					last;
				}
			}

			$count{$p}{$mutateTP}{"all"} ++;
			if($temp2[$numTP-1] == 2){
				$count{$p}{$mutateTP}{"fixed_basal"} ++;
			}
		}
	}	
}

# print header
print OUT "pop\tpop_type\tnum_TP\tout\tcoexist_length\tmutate_TP\tall\tfixed_basal\n";

foreach my $p(@pop){
	foreach my $tp(@{$pop2TimePointArray{$p}}){
		if($tp > 0){
			print OUT $p."\t".$pop2pop_type{($p%100)}."\t".$numT{$p}."\t";
			print OUT $cladeHMMResult{$p}."\t".$coexistLength{$p}."\t".$tp."\t";
			print OUT (0+$count{$p}{$tp}{"all"})."\t".(0+$count{$p}{$tp}{"fixed_basal"})."\n";
		}
	}
}
