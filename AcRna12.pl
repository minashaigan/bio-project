#!d:/Applications/Strawberry/perl/bin/perl.exe
#Rna Secondary structure with ant colony
use warnings;
use strict;

## input : RNA sequence input
	open (my $in,"RnaSeq2.txt") or die "Can't read $!";  
	my $Rna = <$in>;
	$Rna = lc($Rna);
	$Rna=~ s/\R//g;
	my $Structure = <$in>;
	open(my $out,'>',"D:/result.txt") or die "Can't open file for writing: $!"; 
	close $out or die "Failed to close file: $!";
	my @seq1 = split('',$Rna);
	my @seq2 = split('',$Rna);
	my @stems = ();
	my @pheremones = ();
	my @heur_info = ();
	my $IL;
	@stems = Stems();
	@pheremones = Initial_pheremones();
	@heur_info = heuristic_information();
	$IL = max_length();
	my @allowed = ();
	my @solution = ();
	my $Sp;
	my $energy;
	my @ants_solution = ();
	my @subsolution = ();
	my $subsp;
	my $min_energy = 0;
	my $ii = 1;
	my $num_ant = scalar(@seq1);
	my @tempsolution = ();
	my $terminate = 0;
	while($terminate == 0){
		my $kk = 1;
		while($kk<=$num_ant){
			print "\n",$kk,"\n";
			@allowed = (0..scalar(@stems)-1);
			@solution = ();
			push @solution,Initial_stem();
			#push @solution,39;
			@allowed = update_allowed();
			#print "allowed : \n";
			#for(my $i=0;$i<scalar(@allowed);$i++){
			#	print $allowed[$i],":",$stems[$allowed[$i]]{u1},",",$stems[$allowed[$i]]{u2},",",$stems[$allowed[$i]]{u3},"\n";
			#}
			
			#print "solution : \n";
			#for(my $i=0;$i<scalar(@solution);$i++){
			#	print $solution[$i],":",$stems[$solution[$i]]{u1},",",$stems[$solution[$i]]{u2},",",$stems[$solution[$i]]{u3},"\n";
			#}
			#<STDIN>;
			SelectNextStem();
			$Sp = Structure();
			$energy = CompleteEnergy($Sp,$Rna);
			if($energy<$min_energy-$min_energy*0.1){
				$min_energy = $energy;
				@subsolution = @solution;
				$subsp = $Sp;
			}
			$ants_solution[$kk-1] = [@solution];
			$kk++;
			#print "solution : \n";
			#for(my $i=0;$i<scalar(@solution);$i++){
			#	print $solution[$i],":",$stems[$solution[$i]]{u1},",",$stems[$solution[$i]]{u2},",",$stems[$solution[$i]]{u3},"\n";
			#}
			#print "energy",$energy,"\n";
			#print $Sp,"\n";
			#print $Rna,"\n";
			#<STDIN>;
		}
		if(@tempsolution == @subsolution){
			$terminate = 1;
		}
		@tempsolution = @subsolution;
		print_substr();
		Updatepheromone();
		$ii++;
	}
	print_result();
	#print_Stems_inFile();
	#print_Pheremones_inFile();
	#print_heuristic_information_inFile();
	
#Functions 

# Create Dot Plot
sub DotPlot {	
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq1);$i++)
	{
		for(my $j=0;$j<scalar(@seq2);$j++)
		{
			if(($seq1[$i] eq "g" and $seq2[$j] eq "c") or ($seq1[$i]eq"a" and $seq2[$j]eq"u") or ($seq1[$i]eq"c" and $seq2[$j]eq"g") or ($seq1[$i]eq"u" and $seq2[$j]eq"a"))
			{
				$matrix[$i][$j]="\\";
			}
			else
			{
				$matrix[$i][$j]=" ";
			}
		}
		
	}
	return @matrix;
}
#End Create Dot Plot

# Declare Stems
sub Stems {
	my @matrix = DotPlot();
	# first declare stems with length 3
	# u1 : initial nucletide position
	# u2 : final nucletide position
	# u3 : length of the stem
	for(my $i=0;$i<scalar(@seq1);++$i){
		for(my $j=scalar(@seq2)-1;$j>=0;--$j){
			if($j<=$i){
				last;
			}
			if($matrix[$i][$j] eq "\\"){
				my $start=$i;
				my $end=$j;
				my $k;
				for($k=1;$k<scalar(@seq1);++$k){
					if($matrix[$i+$k][$j-$k] eq " "){
						if(abs($i+$k-1-($j-$k+1)+1) < 3){
							$k=0;
						}
						last;
					}
					if($j-$k<=$i+$k){
						$k=0;
						last;
					}
				}
				my $length=$k;
				if($length>=3){
					push @stems, {u1=>$start, u2=>$end, u3=>$length };
				}
			}
		}
	}
	my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	@stems = @sorted;
	return @stems;
}
#End Declare Stems

#print stems in file

#print_Stems_inFile();

sub print_Stems_inFile {
	open(my $out,'>',"D:/Stems.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++)
	{
		print $out $i," => ";
		print $out $stems[$i]{u1}," ",$stems[$i]{u2}," ",$stems[$i]{u3}," ";
		for(my $j=0;$j<$stems[$i]{u3};$j++){
			print $out $seq1[$stems[$i]{u1}+$j];
		}
		print $out " ";
		for(my $j=$stems[$i]{u3}-1;$j>=0;$j--){
			print $out $seq2[$stems[$i]{u2}-$j];
		}
		print $out "\n";
	}
	close $out or die "Failed to close file: $!";
}

#End print stems in file

#calculate initial pheremones

sub Initial_pheremones {
	for(my $i=0;$i<scalar(@stems);$i++){
		my $num = 0;
		for(my $j=0;$j<scalar(@stems);$j++){
			if( (($stems[$i]{u1}+$stems[$i]{u3}<=$stems[$j]{u1} and $stems[$j]{u2}<=$stems[$i]{u2}-$stems[$i]{u3} ) or ( $stems[$i]{u2}<$stems[$j]{u1} ) ) or 
				(($stems[$j]{u1}+$stems[$j]{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$stems[$j]{u2}-$stems[$j]{u3} ) or ( $stems[$j]{u2}<$stems[$i]{u1} )) ){
				$num++;
			}
		}
		#print $i,"= ",$num,"\n";
		for(my $j=0;$j<scalar(@stems);$j++){
			if( (($stems[$i]{u1}+$stems[$i]{u3}<=$stems[$j]{u1} and $stems[$j]{u2}<=$stems[$i]{u2}-$stems[$i]{u3} ) or ( $stems[$i]{u2}<$stems[$j]{u1} ) ) or 
				(($stems[$j]{u1}+$stems[$j]{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$stems[$j]{u2}-$stems[$j]{u3} ) or ( $stems[$j]{u2}<$stems[$i]{u1} )) ){
			
				$pheremones[$i][$j] = 1/$num;
			}
			else {
				$pheremones[$i][$j] = 0;
			}
		}
	}
	return @pheremones;
}

#End calculate initial pheremones

#print stems in file

#print_Pheremones_inFile();

sub print_Pheremones_inFile {
	open(my $out,'>',"D:/Pheremones.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			my $result = sprintf("%.6f", $pheremones[$i][$j]);
			print $out $result," ";
		}
		print $out "\n";
	}
	close $out or die "Failed to close file: $!";
}

#End print stems in file

#caculate heuristicinformation

sub heuristic_information {
	for(my $i=0;$i<scalar(@stems);$i++){
		my $num = 0;
		for(my $j=0;$j<scalar(@stems);$j++){
			if( (($stems[$i]{u1}+$stems[$i]{u3}<=$stems[$j]{u1} and $stems[$j]{u2}<=$stems[$i]{u2}-$stems[$i]{u3} ) or ( $stems[$i]{u2}<$stems[$j]{u1} )) or 
				(($stems[$j]{u1}+$stems[$j]{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$stems[$j]{u2}-$stems[$j]{u3} ) or ( $stems[$j]{u2}<$stems[$i]{u1} ))){
				$heur_info[$i][$j] = $stems[$j]{u3}*$stems[$j]{u3}/$stems[$i]{u3};
			}
			else {
				$heur_info[$i][$j] = 0;
			}
		}
	}
	return @heur_info;
}

#End caculate heuristic information

#print heuristic information in file

#print_heuristic_information_inFile();

sub print_heuristic_information_inFile {
	open(my $out,'>',"D:/heur_info.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			my $result = sprintf("%.2f", $heur_info[$i][$j]);
			print $out $result," ";
		}
		print $out "\n";
	}
	close $out or die "Failed to close file: $!";
}

#End print heuristic information in file

#max length

sub max_length {
	my $sum2 = scalar(@stems);
	my $landa = 0.13;
	my $i = scalar(@stems)-1;
	my $sum1 = 0;
	#my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	for($IL=$stems[scalar(@stems)-1]{u3};$IL>=3;$IL--){
		while($i>=0){
			if($stems[$i]{u3}>=$IL)
				{$sum1++;}
			else{last;}
			$i--;
		}
		#find maximum value of IL that satisfy in expression below
		if($sum1>=$landa*$sum2){
			last;
		}	
	}
	return $IL;
}

#End max length

#select initial stem whose length>=IL

sub Initial_stem {
	my @sel_stems = ();
	#my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	for(my $i=0;$i<scalar(@stems);$i++){
		if($stems[$i]{u3}>=$IL){
			push @sel_stems,$i;
		}
	}
	my $rnd_nm = $sel_stems[0] + int( rand( $sel_stems[scalar(@sel_stems)-1] - $sel_stems[0] ) );
	#my %first_stem = (i => $sep+1+$rnd_nm, u1 => $sel_stems[$rnd_nm]{u1}, u2 => $sel_stems[$rnd_nm]{u2}, u3 => $sel_stems[$rnd_nm]{u3});
	my $first_stem = $rnd_nm;
	print "\n",$stems[$first_stem]{u1},",",$stems[$first_stem]{u2},",",$stems[$first_stem]{u3},"\n";
	return $first_stem;
}

#End select initial stem whose length>=IL

#update allowed base on solution

sub update_allowed {
	my @update_allowed = ();
	for(my $i=0;$i<scalar(@allowed);$i++){
		my $flag;
		for(my $j=0;$j<scalar(@solution);$j++){
			#print "sol",$stems[$solution[$j]]{u1},",",$stems[$solution[$j]]{u2},",",$stems[$solution[$j]]{u3},"\n";
			#print "all",$stems[$allowed[$i]]{u1},",",$stems[$allowed[$i]]{u2},",",$stems[$allowed[$i]]{u3},"\n";
			if( (($stems[$solution[$j]]{u1}+$stems[$solution[$j]]{u3}<=$stems[$allowed[$i]]{u1}) and ($stems[$allowed[$i]]{u2}<=$stems[$solution[$j]]{u2}-$stems[$solution[$j]]{u3}) ) or ( $stems[$solution[$j]]{u2}<$stems[$allowed[$i]]{u1} ) 
				or (($stems[$allowed[$i]]{u1}+$stems[$allowed[$i]]{u3}<=$stems[$solution[$j]]{u1}) and ($stems[$solution[$j]]{u2}<=$stems[$allowed[$i]]{u2}-$stems[$allowed[$i]]{u3}) ) or ( $stems[$allowed[$i]]{u2}<$stems[$solution[$j]]{u1} ) ){
					$flag=1;
			}
			else {
				$flag=0;
				last;
			}
		}
		if($flag==1){
				push @update_allowed, $allowed[$i];
		}
	}		
	#<STDIN>;
	return @update_allowed;
}

#End update allowed base on solution

#Select next stem

sub SelectNextStem{
	while(@allowed){
		my $sum = 0;
		for(my $i=0;$i<scalar(@solution);$i++){
			for(my $j=0;$j<scalar(@allowed);$j++){
				#print "sol ",$pheremones[$solution[$i]][$allowed[$j]],"\n";
				#print "all ",$heur_info[$solution[$i]][$allowed[$j]],"\n";
				$sum += $pheremones[$solution[$i]][$allowed[$j]]*$pheremones[$solution[$i]][$allowed[$j]]*$heur_info[$solution[$i]][$allowed[$j]];
			}
		}
		my $temp = rand(1);
		my $rdm = 4*$temp*(1-$temp);
		my $flag = 0;
		for(my $j=0;$j<scalar(@allowed);$j++){
			my $cur_prob = 0;
			for(my $i=0;$i<scalar(@solution);$i++){
				$cur_prob += $pheremones[$solution[$i]][$allowed[$j]]*$pheremones[$solution[$i]][$allowed[$j]]*$heur_info[$solution[$i]][$allowed[$j]];
			}
			$cur_prob = $cur_prob/$sum;
			if($rdm<$cur_prob){
				push @solution, $allowed[$j];
				@allowed = update_allowed();
				last;
			}
		}
	}
}

#End Select next stem

#find structure of a solution

sub Structure {
	my @sp ;
	for(my $i=0;$i<scalar(@seq1);$i++){
		$sp[$i] = ".";
	}
	for(my $i=0;$i<scalar(@solution);$i++){
		for(my $j=0;$j<$stems[$solution[$i]]{u3};$j++){
			$sp[$stems[$solution[$i]]{u1}+$j]="(";
			$sp[$stems[$solution[$i]]{u2}-$j]=")";
		}
	}
	my $Sp = join("",@sp);
	return $Sp;
}

#End find structure of a solution

# Calculate Energy Solution
sub CompleteEnergy{
	my($SP,$RNA)=@_;
	my($prog);
	my($ES);	
	open(filehandelO3,'>',"solutionstem.txt");
	print(filehandelO3	$RNA,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <solutionstem.txt/;
	print $ES;
	#<STDIN>;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
	
}
# End Calculate Energy Solution

#update pheremones

sub Updatepheromone{
	my $rou = 0.2;
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			$pheremones[$i][$j] = $pheremones[$i][$j]*(0.8);
		}
	}
	for(my $i=0;$i<scalar(@ants_solution);$i++){
		for(my $j=0;$j<scalar(@{$ants_solution[$i]});$j++){
			for(my $k=$j+1;$k<scalar(@{$ants_solution[$i]});$k++){
					$pheremones[$ants_solution[$i][$j]][$ants_solution[$i][$k]] += $energy/($min_energy/0.6);
					$pheremones[$ants_solution[$i][$k]][$ants_solution[$i][$j]] += $energy/($min_energy/0.6);	
					#$pheremones[$ants_solution[$i][$j]][$ants_solution[$i][$k]] += $energy/scalar(@seq1);
					#$pheremones[$ants_solution[$i][$k]][$ants_solution[$i][$j]] += $energy/scalar(@seq1);	
			}
		}
	}
	return @pheremones;
}

#End update pheremones

#print substructure and result 

sub print_substr{
	open(my $out,'>>',"D:/result.txt") or die "Can't open file for writing: $!"; 
	print $out $subsp," -> ",$min_energy,"\n";
	#print $out "sub solution : \n";
	#for(my $i=0;$i<scalar(@subsolution);$i++){
	#	print $out $subsolution[$i],":",$stems[$subsolution[$i]]{u1},",",$stems[$subsolution[$i]]{u2},",",$stems[$subsolution[$i]]{u3},"\n";
	#}
	close $out or die "Failed to close file: $!";
}

sub print_result{
	my $j;
	open(my $out,'>>',"D:/result.txt") or die "Can't open file for writing: $!"; 
	print $out $Structure," : main structure\n";
	print $out "IL : ",$IL,"\n";
	print $out "number of iterations : ",$ii-1,"\n";
	my @Stack = ();
	my @MainStr = ();
	for(my $i=0;$i<length($Rna);++$i){
		if(substr($Structure,$i,1) eq "."){
			$MainStr[$i]=$i
		}
		elsif(substr($Structure,$i,1) eq "("){
			push(@Stack,$i);
		}
		else{
			$j=pop(@Stack);
			$MainStr[$i]=$j;
			$MainStr[$j]=$i;
		}
	}
	@Stack="";
	pop(@Stack);
	my @PredictStr = ();
	for(my $i=0;$i<length($Rna);++$i){
		if(substr($subsp,$i,1) eq "."){
			$PredictStr[$i]=$i
		}
		elsif(substr($subsp,$i,1) eq "("){
			push(@Stack,$i);
		}
		else{
			$j=pop(@Stack);
			$PredictStr[$i]=$j;
			$PredictStr[$j]=$i;
		}
	}
	my $TP=0;
	my $FP=0;
	for(my $i=0;$i<length($Rna);++$i){
		if($PredictStr[$i] eq $MainStr[$i] ){
			++$TP;
		}
		else{
			++$FP;
		}
	}
	print $out $TP,"	",$FP;
		
	close $out or die "Failed to close file: $!";
}

#End print substructure and result 

print "finish";
<>;