#!d:/Applications/Strawberry/perl/bin/perl.exe
#Rna Secondary structure with ant colony
use warnings;
use strict;

## input : RNA sequence input
	open (my $in,"D:/RnaSeq.txt") or die "Can't read $!";  
	my $string = do { local $/; <$in> };
	$string = lc($string);
	my @seq1 = split('',$string);
	my @seq2 = split('',$string);
## find all stems with dotplot
	create_DotPlot(@seq1,@seq2);
## calculate energy of each stems
	my %stack_energy;
	$stack_energy{"aa"}{"uu"} = -0.9;
	$stack_energy{"ac"}{"ug"} = -2.2;
	$stack_energy{"ag"}{"uc"} = -2.1;
	$stack_energy{"ag"}{"uu"} = -0.6;
	$stack_energy{"au"}{"ua"} = -1.1;
	$stack_energy{"au"}{"ug"} = -1.4;
	$stack_energy{"ca"}{"gu"} = -2.1;
	$stack_energy{"cc"}{"gg"} = -3.3;
	$stack_energy{"cg"}{"gc"} = -2.4;
	$stack_energy{"cg"}{"gu"} = -1.4;
	$stack_energy{"cu"}{"ga"} = -2.1;
	$stack_energy{"cu"}{"gg"} = -2.1;
	$stack_energy{"ga"}{"cu"} = -2.4;
	$stack_energy{"ga"}{"uu"} = -1.3;
	$stack_energy{"gc"}{"cg"} = -3.4;
	$stack_energy{"gc"}{"ug"} = -2.5;
	$stack_energy{"gg"}{"cc"} = -3.3;
	$stack_energy{"gg"}{"cu"} = -1.5;
	$stack_energy{"gg"}{"uc"} = -2.1;
	$stack_energy{"gg"}{"uu"} = -0.5;
	$stack_energy{"gu"}{"ca"} = -2.2;
	$stack_energy{"gu"}{"cg"} = -2.5;
	$stack_energy{"gu"}{"ua"} = -1.4;
	$stack_energy{"gu"}{"ug"} = 1.3;
	$stack_energy{"ua"}{"au"} = -1.3;
	$stack_energy{"ua"}{"gu"} = -1;
	$stack_energy{"uc"}{"ag"} = -2.4;
	$stack_energy{"uc"}{"gg"} = -1.5;
	$stack_energy{"ug"}{"ac"} = -2.1;
	$stack_energy{"ug"}{"au"} = -1;
	$stack_energy{"ug"}{"gc"} = -1.4;
	$stack_energy{"ug"}{"gu"} = 0.3;
	$stack_energy{"uu"}{"aa"} = -0.9;
	$stack_energy{"uu"}{"ag"} = -1.3;
	$stack_energy{"uu"}{"ga"} = -0.6;
	$stack_energy{"uu"}{"gg"} = -0.5;
	calculate_Stem_Energy();
## iteration i = 1
### find first stem : initialize IL
	#my %first_stem = find_First_Stem();
### initialize pheremone
	#initialize_first_pheremone_huristic(%first_stem) ;
	#print_first_iteration(%first_stem);
##
	my $Rna = "";
	for(my $i=0;$i<scalar(@seq1);$i++){
		$Rna .= $seq1[$i];
	}
	
	#$Rna .= "&";
	#for(my $i=0;$i<scalar(@seq2);$i++){
	#	$Rna .= $seq2[$i];
	#}
	my @pheremones = ();
	my @allowed = declare_Stems();
	for(my $i=0;$i<scalar(@allowed);$i++){
		my $num = 0;
		for(my $j=0;$j<scalar(@allowed);$j++){
			$pheremones[$i][$j] = 0;
		}
	}
#	print "hh";
#	<STDIN>;
#	my $i; 
#	$i=@allowed;
#	print($i);
#	<STDIN>;
	for(my $i=0;$i<scalar(@allowed);$i++){
		my $num = 0;
#		print $i,"\n";
		for(my $j=0;$j<scalar(@allowed);$j++){
			if( ($allowed[$i]{u1}+$allowed[$i]{u3}<=$allowed[$j]{u1} and $allowed[$j]{u2}<=$allowed[$i]{u2}-$allowed[$i]{u3} ) or ( $allowed[$i]{u2}<$allowed[$j]{u1} ) ){
				$num++;
			}
		}
		for(my $j=0;$j<scalar(@allowed);$j++){
			if( ($allowed[$i]{u1}+$allowed[$i]{u3}<=$allowed[$j]{u1} and $allowed[$j]{u2}<=$allowed[$i]{u2}-$allowed[$i]{u3} ) or ( $allowed[$i]{u2}<$allowed[$j]{u1} )){
				$pheremones[$i][$j] = 1/$num;
			}
		}
	}
	#print "girkardeh";
	#<STDIN>;
=comment
{
	open(my $outtt,'>>',"D:/pheremones.txt") or die;
	print $outtt "first pheremones:\n";
	for(my $i=0;$i<scalar(@allowed);$i++){
		for(my $j=0;$j<scalar(@allowed);$j++){
			my $r = sprintf("%.3f", $pheremones[$i][$j]);
			print $outtt $r," ";
		}
		print $outtt "\n";
	}
	close $outtt or die;
}
=cut
	#print_Stems_inFile();
	my @ants_solution_string =();
		my @ants_solution = ();
		my $count =0;
	my $ii = 1;
	my $best_index;
	my @best;
	while($ii<=15){
=comment
		open(my $outt,'>',"D:/pheremones.txt") or die;
		print $outt "iteration ",$ii,":\n";
		close $outt or die;
=cut
		my $rou = 0.2;
		for(my $i=0;$i<scalar(@allowed);$i++){
			for(my $j=0;$j<scalar(@allowed);$j++){
				$pheremones[$i][$j] = $rou * $pheremones[$i][$j];
			}
		}
		my $min=10000;
			
		for(my $k=1;$k<=10;$k++){
			my %first_stem;
			my @solution;
			if($k ==1 and $ii>1){
				%first_stem = (  u1 => $allowed[$ants_solution[$best_index][0]]{u1}, u2 => $allowed[$ants_solution[$best_index][0]]{u2}, u3 => $allowed[$ants_solution[$best_index][0]]{u3},);
				for(my $i=1;$i<scalar(@best);$i++){
					$solution[$i-1] = $best[$i];
				}
				
				print "best solution : ",@solution;
			}
			else{
				%first_stem = find_First_Stem();
				@solution = select_next_stem(%first_stem);
			}
			my $index=0;
			#find index a first stem
			for(my $i=0;$i<scalar(@allowed);$i++){
				if($first_stem{u1}==$allowed[$i]{u1} and $first_stem{u2}==$allowed[$i]{u2} and $first_stem{u3}==$allowed[$i]{u3}){
					$index = $i;
					last;
				}
			}
			##### print
			####
			print "\n",$k,":\n";
			####
			print "first stem ",$index," : ",$first_stem{u1},",",$first_stem{u2},",",$first_stem{u3},"\n";
			####
			for(my $j=0;$j<scalar(@solution);$j++){
				print $solution[$j],":",$allowed[$solution[$j]]{u1},",",$allowed[$solution[$j]]{u2},",",$allowed[$solution[$j]]{u3},"\n";
			}
			#####
			my @solution_string = ();
			for(my $i=0;$i<scalar(@seq1)+scalar(@seq2);$i++){
				$solution_string[$i]="n";
			}
			#put first stem in solution string
			for(my $j=0;$j<$allowed[$index]{u3};$j++){
					$solution_string[$allowed[$index]{u1}+$j] = $seq1[$allowed[$index]{u1}+$j];
					$solution_string[scalar(@seq1)+$allowed[$index]{u2}-$j] = $seq2[$allowed[$index]{u2}-$j];
				}
			#put soltions in solution string 
			for(my $i=0;$i<scalar(@solution);$i++){
				for(my $j=0;$j<$allowed[$solution[$i]]{u3};$j++){
					$solution_string[$allowed[$solution[$i]]{u1}+$j] = $seq1[$allowed[$solution[$i]]{u1}+$j];
					$solution_string[scalar(@seq1)+$allowed[$solution[$i]]{u2}-$j] = $seq2[$allowed[$solution[$i]]{u2}-$j];
				}
			}
			$ants_solution[$k-1] = [$index,@solution];
			$ants_solution_string[$k-1] = [@solution_string];
			my @sp ;
			for(my $i=0;$i<scalar(@seq1);$i++){
				$sp[$i]=".";
			}
			for(my $j=0;$j<$allowed[$index]{u3};$j++){
					$sp[$allowed[$index]{u1}+$j]="(";
					$sp[$allowed[$index]{u2}-$j]=")";
				}
			for(my $i=0;$i<scalar(@solution);$i++){
				for(my $j=0;$j<$allowed[$solution[$i]]{u3};$j++){
					$sp[$allowed[$solution[$i]]{u1}+$j]="(";
					$sp[$allowed[$solution[$i]]{u2}-$j]=")";
				}
			}
			my $Sp = join("",@sp);
=comment
			for(my $i=0;$i<scalar(@seq1);$i++){
				if($solution_string[$i]eq"n"){
					$Sp .= ".";
				}
				else{
					$Sp .= "(";
				}
			}
			$Sp .= "&";
			for(my $i=scalar(@seq1);$i<scalar(@seq1)+scalar(@seq2);$i++){
				if($solution_string[$i] eq "n"){
					$Sp .= ".";
				}
				else{
					$Sp .= ")";
				}
			}
=cut
		
			#my $solution_string;
			#$solution_string = join("", @solution_string[0..scalar(@seq1)-1]);
			#$solution_string .= "&" ;
			#$solution_string .= join("", @solution_string[scalar(@seq1)..scalar(@seq1)+scalar(@seq2)-1]);
			#print $solution_string;
=comment
			for(my $i=0;$i<scalar(@seq1)+scalar(@seq2);$i++){
				print $solution_string[$i];
			}
=cut
=comment
			open(my $out,'>',"D:/pheremones.txt") or die;
			print $out "pheremones ant",$k,":\n";
			for(my $i=0;$i<scalar(@allowed);$i++){
				for(my $j=0;$j<scalar(@allowed);$j++){
					my $result = sprintf("%.3f", $pheremones[$i][$j]);
					print $out $result," ";
				}
				print $out "\n";
			}
			close $out or die;
=cut
					#my $energy_result = CompleteEnergy($Sp,$Rna);
					#if($energy_result<$min){
					#$min = $energy_result;
					#$best_index = $index;
					#}
		}
		
		#### pheremones updates
			
			for(my $t=0;$t<scalar(@ants_solution);$t++){
				for(my $s=0;$s<scalar(@{$ants_solution[$t]});$s++){
					my @sp ;
					for(my $i=0;$i<scalar(@seq1);$i++){
						$sp[$i]=".";
					}
					for(my $j=0;$j<$allowed[$ants_solution[$t][0]]{u3};$j++){
							$sp[$allowed[$ants_solution[$t][0]]{u1}+$j]="(";
							$sp[$allowed[$ants_solution[$t][0]]{u2}-$j]=")";
						}
					print "\n";
					for(my $i=1;$i<scalar(@{$ants_solution[$t]});$i++){
						for(my $j=0;$j<$allowed[$ants_solution[$t][$i]]{u3};$j++){
							$sp[$allowed[$ants_solution[$t][$i]]{u1}+$j]="(";
							$sp[$allowed[$ants_solution[$t][$i]]{u2}-$j]=")";
						}
					}
					my $Sp = join("",@sp);
					my $energy_result = CompleteEnergy($Sp,$Rna);
					if($energy_result<$min){$min = $energy_result;$best_index = $t;@best = @{$ants_solution[$t]};}
					for(my $i=0;$i<scalar(@{$ants_solution[$t]});$i++){
						for(my $j=$i+1;$j<scalar(@{$ants_solution[$t]});$j++){
							$pheremones[$ants_solution[$t][$i]][$ants_solution[$t][$j]] += $energy_result/scalar(@seq1);
						}
					}
				}
			}
		$ii++;
		print $min;
		my @sp ;
		for(my $i=0;$i<scalar(@seq1);$i++){
			$sp[$i]=".";
		}
		for(my $j=0;$j<$allowed[$best[0]]{u3};$j++){
				$sp[$allowed[$best[0]]{u1}+$j]="(";
				$sp[$allowed[$best[0]]{u2}-$j]=")";
			}
		print "\n";
		for(my $i=1;$i<scalar(@best);$i++){
			for(my $j=0;$j<$allowed[$best[$i]]{u3};$j++){
				$sp[$allowed[$best[$i]]{u1}+$j]="(";
				$sp[$allowed[$best[$i]]{u2}-$j]=")";
			}
		}
		my $Sp = join("",@sp);
		print $Sp;
		<STDIN>;
		}
		
		#update pheremones
=comment
		print "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		for(my $i=0;$i<scalar(@ants_solution_string);$i++){
			my @temp = @{$ants_solution_string[$i]};
			my $Sp = "";
			for(my $i=0;$i<scalar(@seq1);$i++){
				if($temp[$i]eq"n"){
					$Sp .= ".";
				}
				else{
					$Sp .= "(";
				}
			}
			$Sp .= "&";
			for(my $i=scalar(@seq1);$i<scalar(@seq1)+scalar(@seq2);$i++){
				if($temp[$i] eq "n"){
					$Sp .= ".";
				}
				else{
					$Sp .= ")";
				}
			}
			my $result = CompleteEnergy($Sp,$Rna);
			my @fields = split /\s[(]\s*/, $result;
			my @result = split /[)]/, $fields[1];
			print "\n***\n",$result,"\n***\n";
			print "\n***\n",$result[0],"\n***\n";
		}

	#}
	#print_Solution_String_inFile();
=cut
#Functions
sub print_Seq1 {
	for my $nuc (@seq1){	
		print $nuc;}	
}
sub print_Seq2 {
	for my $nuc (@seq2) {
		print $nuc;}
}

sub create_DotPlot {	
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq1);$i++)
	{
		for(my $j=0;$j<scalar(@seq2);$j++)
		{
			if(($seq1[$i] eq "g" and $seq2[$j] eq "c") or ($seq1[$i]eq "g" and $seq2[$j]eq"u") or ($seq1[$i]eq"a" and $seq2[$j]eq"u") or ($seq1[$i]eq"c" and $seq2[$j]eq"g") or ($seq1[$i]eq"u" and $seq2[$j]eq"a") or ($seq1[$i]eq"u" and $seq2[$j]eq"g"))
			{
				if($i != $j){
					$matrix[$i][$j]="\\";
				}
			}
			else
			{
				$matrix[$i][$j]=" ";
			}
		}
		
	}
	return @matrix;
}
sub print_Matrix {
	my @matrix = create_DotPlot();
	for(my $i=0;$i<scalar(@seq1);$i++){
		for(my $j=0;$j<scalar(@seq2);$j++){
			print $matrix[$i][$j];
		}
		print "\n";
	}
}
sub declare_Stems {
	my @matrix = create_DotPlot();
	my @stems = ();
	# first declare stems with length 3
	# u1 : initial rebonucletide position
	# u2 : final rebonucletide position
	# u3 : length of the stem
	my $length = 3;
	my $count = 0;
	
	for(my $i = 0; $i<scalar(@seq1); ++$i)
	{
		for(my $j = $i ; $j<scalar(@seq2); ++$j)
		{
			for(my $w = 0 ; $w<$length; ++$w){
				if(($i+$w)<scalar(@seq1) and ($j-$w)<scalar(@seq2))
				{
					if($matrix[$i+$w][$j-$w] eq "\\" and $i != $j)
						{
							++$count;
						}
				}
			}
			if($count == 3)
			{
				if($i+$length<$j-$length){
					push @stems, {u1=>$i, u2=>$j, u3=>$length };
					
				}
			}
			$count=0;	
		}
		
	}
	my @temp_array = ();
	for(my $i=0;$i<scalar(@stems);++$i)
		{	
		my $pos1 = $stems[$i]{u1}+2;
		my $pos2 = $stems[$i]{u2}-2;
		#my $num=0;
		for(my $length = 4;$length<scalar(@seq1);++$length)
			{
				if(($pos1)<scalar(@seq1)-1 and ($pos2)<scalar(@seq2)-1)
				{
					if($matrix[++$pos1][--$pos2] eq "\\" and $pos1 != $pos2)
					{
						if($stems[$i]{u1}+$length<$stems[$i]{u2}-$length){
							push @temp_array, {u1=>$stems[$i]{u1}, u2=>$stems[$i]{u2}, u3=>$length };
							#print $stems[$i]{u1};
							#$num++;
						}
					}
					else 
					{
					last;
					}
				}
			}
		}

	my @sorted_temp_array =  sort { $a->{u3} <=> $b->{u3} } @temp_array;
	push(@stems,@sorted_temp_array);
	
	return @stems;
}
sub print_Stems_inFile {
	my @stems = declare_Stems();
	open(my $out,'>>',"D:/output_RnaSeq.txt") or die "Can't open file for writing: $!"; 
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

sub print_Energy_inFile {
	open (my $out,'>>',"D:/stack_energy.txt") or die;
	foreach my $top (sort { $a cmp $b} keys %stack_energy)
	{
		foreach my $bottom ( keys %{$stack_energy{$top}})
		{print $out "$top, $bottom: $stack_energy{$top}{$bottom}\n";}
	}
	close $out or die;
}
sub calculate_Stem_Energy {
	my @stems = declare_Stems();
	my @energys = ();
	for(my $i=0;$i<scalar(@stems);++$i){
		my $sum_energy = 0;
		for(my $j=0;$j<$stems[$i]{u3}-1;++$j){
			foreach my $top (sort { $a cmp $b} keys %stack_energy){
				foreach my $bottom ( keys %{$stack_energy{$top}}){
					my $pair1 = join('',$seq1[$stems[$i]{u1}+$j],$seq1[$stems[$i]{u1}+$j+1]);
					my $pair2 = join('',$seq2[$stems[$i]{u2}-$stems[$i]{u3}+$j+1],$seq2[$stems[$i]{u2}-$stems[$i]{u3}+$j+2]);
					if($pair1 eq $top and $pair2 eq $bottom)
					{
						$sum_energy = $sum_energy + $stack_energy{$top}{$bottom};
					}
				}
			}	
		}
		push @energys,$sum_energy;
	}
	return @energys;
}
sub print_Stems_Energy_inFile {
	my @stems = declare_Stems();
	my @energys = calculate_Stem_Energy();
	my $i=0;
	open(my $out,'>>',"D:/output_RnaSeq_Energy.txt") or die;
	foreach my $energy(@energys)
	{
		print $out $i+1," => ";
		print $out $stems[$i]{u1}," ",$stems[$i]{u2}," ",$stems[$i]{u3};
		print $out " : $energy\n";
		++$i;
	}
	close $out or die;
}
sub find_First_Stem {
	my $IL;
	my @lengths = ();
	my @stems = declare_Stems();
	my $sum2 = scalar(@stems);
	my $landa = 0.13;
	my $i = scalar(@stems)-1;
	my $sum1 = 0;
	for($IL=$stems[scalar(@stems)-1]{u3};$IL>=3;$IL--){
		for(;$i>=0;$i--){
			if($stems[$i]{u3}>=$IL)
				{$sum1++;}
			else{last;}
		}
		#find maximum value of IL that satisfy in expression below
		if($sum1>=$landa*$sum2){
			last;}	
	}
	#print $IL; #4
	my @temp_stems = ();
	for(my $j = scalar(@stems)-1;$j>=0;$j--){
		if($stems[$j]{u3}>=$IL){
		push @temp_stems,$stems[$j];
		}
	}
	#select randomly from stems whose length is equal to or greater than IL
	my $rnd_nm = 0 + int(rand( scalar(@temp_stems) - 0));
	my %first_stem = (  u1 => $temp_stems[$rnd_nm]{u1}, u2 => $temp_stems[$rnd_nm]{u2}, u3 => $temp_stems[$rnd_nm]{u3},);
	#print $first_stem{u1},",",$first_stem{u2},",",$first_stem{u3};
	return %first_stem;
}
sub initialize_first_pheremone_huristic {
	my @stems = declare_Stems();
	my @first_pheremones = ();
	my %first_stem = @_;
	my $index;
	for(my $i=0;$i<scalar(@stems);$i++){
		if($first_stem{u1}==$stems[$i]{u1} and $first_stem{u2}==$stems[$i]{u2} and $first_stem{u3}==$stems[$i]{u3}){
			$index = $i;
			last;
		}
	}
	#print $first_stem{u1},",",$first_stem{u2},",",$first_stem{u3};
	my $num = 0;
	for(my $i=0;$i<scalar(@stems);$i++){
		if( ($first_stem{u1}+$first_stem{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$first_stem{u2}-$first_stem{u3} ) or ( $first_stem{u2}<$stems[$i]{u1} )){
			$num++;
		}
	}
	#print "\n",$num,"\n";
	for(my $i=0;$i<scalar(@stems);$i++){
		if( ($first_stem{u1}+$first_stem{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$first_stem{u2}-$first_stem{u3}) or ( $first_stem{u2}<$stems[$i]{u1} ) ){
			my $heuristic = $stems[$i]{u3}*$stems[$i]{u3}/$first_stem{u3};
			push @first_pheremones, {i=>$i, p=>1/$num, h=>$heuristic };			
		}
	}
	return @first_pheremones;
}
sub print_Pheremones_inFile {
	my @stems = declare_Stems();
	open(my $out,'>>',"D:/pheremones.txt") or die;
	for(my $i=0;$i<scalar(@stems);$i++)
	{
		for(my $j=0;$j<scalar(@stems);$j++){
			my $result = sprintf("%.2f", $pheremones[$i][$j]);
			print $out $result," ";
		}
		print $out "\n";
	}
	close $out or die;
}
sub print_first_iteration {
	my %first_stem = @_;
	my @first_pheremones = initialize_first_pheremone_huristic(%first_stem);
	print $first_stem{u1},",",$first_stem{u2},",",$first_stem{u3},"\n";
	for(my $i=0;$i<scalar(@first_pheremones);$i++){
		print $first_pheremones[$i]{i},",";
	}
	print "\n";
}
sub initialize_pheremone_huristic {
	my @stems = declare_Stems();
	my %first_stem = @_;
	#my @first_pheremones = initialize_first_pheremone_huristic(%first_stem);
	################################################################################################
	my @first_pheremones = ();
	my $index;
	for(my $i=0;$i<scalar(@stems);$i++){
		if($first_stem{u1}==$stems[$i]{u1} and $first_stem{u2}==$stems[$i]{u2} and $first_stem{u3}==$stems[$i]{u3}){
			$index = $i;
			last;
		}
	}
	#print "\n",$num,"\n";
	for(my $i=0;$i<scalar(@stems);$i++){
		if( ($first_stem{u1}+$first_stem{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$first_stem{u2}-$first_stem{u3}) or ( $first_stem{u2}<$stems[$i]{u1}) ){
			my $heuristic = $stems[$i]{u3}*$stems[$i]{u3}/$first_stem{u3};
			push @first_pheremones, {i=>$i, p=>$pheremones[$index][$i], h=>$heuristic };			
		}
	}
	#######################################################################################
	my @probabilities = ();
	#print $first_stem{u1},",",$first_stem{u2},",",$first_stem{u3};
	my $sum = 0;
	for(my $i=0;$i<scalar(@first_pheremones);$i++){
		$sum += $first_pheremones[$i]{p}*$first_pheremones[$i]{h};
	}
	for(my $i=0;$i<scalar(@first_pheremones);$i++){
		push @probabilities, {index=>$first_pheremones[$i]{i}, p=>($first_pheremones[$i]{p}*$first_pheremones[$i]{h})/$sum };
	}
	print "\n";
	return @probabilities;
}
sub select_next_stem {
	my %first_stem = @_;
	my @probabilities = initialize_pheremone_huristic(%first_stem);
	my @solution = ();
	my $count = 0;
	my @stems = declare_Stems();
	while(@probabilities){
		my $sum = 0;
		for(my $i=0;$i<scalar(@probabilities);$i++){
			$sum += $probabilities[$i]{p};
		}
		my $rdm = rand(1);
		my $pheremone_first_stem = 0;
		for(my $i=0;$i<scalar(@probabilities);$i++){
			$pheremone_first_stem += $probabilities[$i]{p}/$sum;
			if($rdm<$pheremone_first_stem){
				push @solution,($probabilities[$i]{index});
				print $probabilities[$i]{index};
#				<STDIN>;
				for(my $t=0;$t<scalar(@probabilities);$t++){
					my $j = $solution[$count];
					my $w = $probabilities[$t]{index};
					if( ( $stems[$j]{u1}+$stems[$j]{u3}<=$stems[$w]{u1} and $stems[$w]{u2}<=$stems[$j]{u2}-$stems[$j]{u3} ) or ( $stems[$j]{u2}<$stems[$w]{u1} ) ){}
					else{
						splice @probabilities, $t, 1;
						$t--;
					}
				}
				$count++;
				last;
			}		
		}
	}
	return @solution;
}
=comment
sub print_Solution_String_inFile {
	my @solution_string = @_;
	my $Rna = "";
	my $Sp = "";
	open(my $out,'>',"D:/solutionstem.txt") or die;
	for(my $i=0;$i<scalar(@seq1);$i++){
		print $out $seq1[$i];
		$Rna .= $seq1[$i];
	}
	print $out "&";
	$Rna .= "&";
	for(my $i=0;$i<scalar(@seq2);$i++){
		print $out $seq2[$i];
		$Rna .= $seq2[$i];
	}
	print $out "\n";
	for(my $i=0;$i<scalar(@seq1);$i++){
		if($solution_string[$i]eq"n"){
			print $out ".";
			$Sp .= ".";
		}
		else{
			print $out "(";
			$Sp .= "(";
		}
	}
	print $out "&";
	$Sp .= "&";
	for(my $i=scalar(@seq1);$i<scalar(@seq1)+scalar(@seq2);$i++){
		if($solution_string[$i] eq "n"){
			print $out ".";
			$Sp .= ".";
		}
		else{
			print $out ")";
			$Sp .= ")";
		}
	}
	print $out "\n@";
	close $out or die;
	return ($Sp,$Rna);
}
=cut
sub CompleteEnergy{
	my($SP,$RNA)=@_;
	my($prog);
	my($ES);	
	open(filehandelO3,'>',"solutionstem.txt");
	print(filehandelO3	$RNA,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <solutionstem.txt/;
	print $ES;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
}

<>;
