use warnings;
use strict;
main();


sub main()
{
	my($seqRNA);
	my(@ListStem);
	my(@Eta);
	my(@Teta);
	my(@Ant);
	my($FirstStem);
	my(@Allowed);
	my(@Pop);
	my($k);
	my($i);
	my($j);
	my(@Delta);
	my(@T);
	my($Min)=0;
	my($MinTotal)=0;
	my($GTotal);
	my($General);
    unless(open(fIN,"RnaSeq2.txt"))
    {
        print("RnaSeq2.txt");
        exit;
    }
    unless(open(fO,">1.txt"))
    {
        print("1.txt");
        exit;
    }
	$seqRNA=<fIN>;
	$seqRNA=uc($seqRNA);
	chomp($seqRNA);
	ExtractStems(\@ListStem,\@Eta,$seqRNA);
	IniPhe(\@ListStem,\@Teta);
	for($k=0;$k<10;++$k)
	{
		print($k,"\n");
		$Min=0;
		@Pop="";
		pop(@Pop);
		for($i=0;$i<@ListStem;++$i)
		{
			for($j=0;$j<@ListStem;++$j)
			{
				$Delta[$i][$j]=0;
			}
		}
		for($i=0;$i<20;++$i)
		{
			@Ant="";
			pop(@Ant);
			$FirstStem=FindFirstStem(\@ListStem);
			push(@Ant,$FirstStem);
			#push(@Ant,31);
			print "first : ",$ListStem[$Ant[0]],"\n";
			@Allowed=MakeAllowed(\@ListStem,\@Ant);
			for(my $rt=0;$rt<scalar(@Allowed);$rt++){
				if($Allowed[$rt]==1){
					print $ListStem[$rt],"\n";
				}
			}
			#<STDIN>;
			SelectNextStem(\@Ant,\@Allowed,\@ListStem,\@Eta,\@Teta);
			push(@Pop,MakeStr(\@Ant,$seqRNA,\@ListStem));
			@T=split(",",$Pop[@Pop-1]);			
			if($T[1]<$Min)
			{
				$General=$Pop[@Pop-1];
				$Min=$T[1];
				if($MinTotal==0)
				{
					$MinTotal=$Min;
					$GTotal=$General;
				}
			}
#			MakeDelta(\@Ant,\@Delta,$T[1]/length($seqRNA));
			MakeDelta(\@Ant,\@Delta,$T[1]/($MinTotal/0.6));
			for(my $rt=0;$rt<scalar(@Ant);$rt++){
				print $ListStem[$Ant[$rt]],"\n";
			}
			print "energy",$T[1],"\n";
			#<STDIN>;
		}
#		print("dddd=",$General,"\n");
#		<STDIN>;
		if($MinTotal >$Min)
		{
			$GTotal=$General;
			$MinTotal=$Min;
			$k=0
		}
#			print("fff=",$GTotal,"\n");
#			<STDIN>;
		Update(\@Teta,\@Delta);
	}
	for($i=0;$i<@Pop;++$i)
	{
		print(fO $Pop[$i],"\n");
	}

	print(fO "==================\n");
	print(fO $GTotal,"\n");
	$seqRNA=<fIN>;
	chomp($seqRNA);
	my(@Stack);
	my(@MainStr);
	for($i=0;$i<length($seqRNA);++$i)
	{
		if(substr($seqRNA,$i,1) eq ".")
		{
			$MainStr[$i]=$i
		}
		elsif(substr($seqRNA,$i,1) eq "(")
		{
			push(@Stack,$i);
		}
		else
		{
			$j=pop(@Stack);
			$MainStr[$i]=$j;
			$MainStr[$j]=$i;
		}
	}
	@Stack="";
	pop(@Stack);
	my(@PredictStr);
	for($i=0;$i<length($seqRNA);++$i)
	{
		if(substr($GTotal,$i,1) eq ".")
		{
			$PredictStr[$i]=$i
		}
		elsif(substr($GTotal,$i,1) eq "(")
		{
			push(@Stack,$i);
		}
		else
		{
			$j=pop(@Stack);
			$PredictStr[$i]=$j;
			$PredictStr[$j]=$i;
		}
	}
	my($TP)=0;
	my($FP)=0;
	for($i=0;$i<length($seqRNA);++$i)
	{
		if($PredictStr[$i] eq $MainStr[$i] )
		{
			++$TP;
		}
		else
		{
			++$FP;
		}
	}
	print(fO $TP,"	",$FP);
}
sub Update()
{
	my($Teta,$Delta)=@_;
	my($i);
	my($j);
	for($i=0;$i<@$Teta;++$i)
	{
		for($j=0;$j<@$Teta;++$j)
		{
			$$Teta[$i][$j]=0.8*$$Teta[$i][$j]+$$Delta[$i][$j];
		}
	}
}

sub MakeDelta()
{
	my($Ant,$Delta,$Val)=@_;
	my($i);
	my($j);
	for($i=0;$i<@$Ant;++$i)
	{
		for($j=$i;$j<@$Ant;++$j)
		{
			$$Delta[$$Ant[$i]][$$Ant[$j]]+=$Val;
			$$Delta[$$Ant[$j]][$$Ant[$i]]+=$Val;
		}
	}
}
sub MakeStr()
{
	my($Ant,$seqRNA,$ListStem)=@_;
	my($i);
	my($j);
	my(@Str);
	my(@T);
	my($E);
	for($i=0;$i<length($seqRNA);++$i)
	{
		$Str[$i]=".";
	}

	for($i=0;$i<@$Ant;++$i)
	{
		@T=split(",",$$ListStem[$$Ant[$i]]);
		for($j=0;$j<$T[2];++$j)
		{
			if($Str[$T[0]+$j] eq ".")
			{
				$Str[$T[0]+$j]="(";
			}
			else
			{
				print("error1");
				<STDIN>;
			}
			if($Str[$T[1]-$j] eq ".")
			{
				$Str[$T[1]-$j]=")";
			}
			else
			{
				print("error2");
				<STDIN>;
			}
		}
	}
	$E=CompleteEnergy(join("",@Str),$seqRNA);
	print join("",@Str),"\n";
	return(join("",@Str).",".$E);
}
sub CompleteEnergy{
	my($SP,$RNA)=@_;
	my($prog);
	my($ES);	
	open(filehandelO3,'>',"solutionstem.txt");
	print(filehandelO3	$RNA,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <solutionstem.txt/;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
}


sub SelectNextStem()
{
	my($Ant,$Allowed,$ListStem,$Eta,$Teta)=@_;
	my($Sum1);
	my($Sum);
	my($i);
	my($j);
	my($t);
	my($flag);
	my(@Prob);
	my($k);
	while(NumberofOne($Allowed)!=0)
	{
		@Prob=P($ListStem,$Ant,$Allowed,$Eta,$Teta);
		$Sum1=0;
		for($k=0;$k<@$Ant;++$k)
		{
			for($j=0;$j<@$ListStem;++$j)
			{
				if($$Allowed[$j]==1)
				{
					$Sum1+=$Prob[$$Ant[$k]][$j];
				}
			}
		}
	#	$Sum=0;
	#	$flag=0;
	#	$t=rand(1);
	#	for($k=0;$k<@$Ant;++$k)
	#	{
	#		for($j=0;$j<@$ListStem;++$j)
	#		{
	#			if($$Allowed[$j]==1)
	#			{	
	#				$Sum+=$Prob[$$Ant[$k]][$j]/$Sum1;
	#				if($t<$Sum)
	#				{
	#					push(@$Ant,$j);
	#					@$Allowed=MakeAllowed($ListStem,$Ant);
	#					$flag=1;
	#					last;
	#				}
	#			}
	#		}
	#		if($flag==1)
	#		{
	#			last;
	#		}
	#	}
		$Sum=0;
		$flag=0;
		$t=rand(1);
		$t = 4*$t*(1-$t);
		for($j=0;$j<@$ListStem;++$j)
		{
			$Sum=0;
			if($$Allowed[$j]==1)
			{
				for($k=0;$k<@$Ant;++$k)
				{
					$Sum+=$Prob[$$Ant[$k]][$j];
				}
				$Sum=$Sum/$Sum1;
				if($t<$Sum)
				{
					push(@$Ant,$j);
					@$Allowed=MakeAllowed($ListStem,$Ant);
					last;
				}
			}
		}
	}
	
}
sub NumberofOne()
{
	my($Allowed)=@_;
	my($i);
	for($i=0;$i<@$Allowed;++$i)
	{
		if($$Allowed[$i]==1)
		{
			return(1)
		}
	}
	return(0);
}
sub P()
{
	my($ListStem,$Ant,$Allowed,$Eta,$Teta)=@_;
	my(@Prob);
	my($Sum);
	my($i);
	my($j);
	my($k);
	for($k=0;$k<@$Ant;++$k)
	{
		$i=$$Ant[$k];
		$Sum=0;
		for($j=0;$j<@$ListStem;++$j	)
		{
			if($$Allowed[$j]==1)
			{
				$Sum+=$$Teta[$i][$j]*$$Teta[$i][$j]*$$Eta[$i][$j];
				if($Sum==0)
				{
					print("error");
					<STDIN>;
				}
			}
		}

		for($j=0;$j<@$ListStem;++$j	)
		{
			$Prob[$i][$j]=0;
			if($$Allowed[$j]==1)
			{
				$Prob[$i][$j]=$$Teta[$i][$j]*$$Teta[$i][$j]*$$Eta[$i][$j];
			}
		}
	}
	return(@Prob);
}



sub MakeAllowed()
{
	my($ListStem,$Ant)=@_;
	my(@Allowed);
	my($j);
	my($i);
	for($i=0;$i<@$ListStem;++$i)
	{
		for($j=0;$j<@$Ant;++$j)
		{
			if(Consistance($$ListStem[$i],$$ListStem[$$Ant[$j]])==0)
			{
				last;
			}
		}
		if($j==@$Ant)
		{
			$Allowed[$i]=1;
		}
		else
		{
			$Allowed[$i]=0;
		}
	}
	return(@Allowed);
}




sub FindFirstStem()
{
	my($ListStem)=@_;
	my($i);
	my($Max)=-1;
	my($Landa)=0.13;
	my($IL);
	my(@T);
	my(@List);

	my(%Hash);
	my($Sum);
	my($Sum1);
	
	for($i=0;$i<@$ListStem;++$i)
	{
		@T=split(",",$$ListStem[$i]);
		$Hash{$T[2]}++;
		if($T[2]>$Max)
		{
			$Max=$T[2];
		}

	}
	$Sum=@$ListStem;
	$Sum=$Landa*$Sum;
	$Sum1=0;
	for($i=$Max;$i>=0;--$i)
	{
		$Sum1+=$Hash{$i};
		if($Sum1>=$Sum)
		{
			last;
		}
	}
	$IL=$i;
	@List=@$ListStem;
	MergeSort(\@List,0,@List-1);
	for($i=@List-1;$i>=0;--$i)
	{
		@T=split(",",$List[$i]);
		if($T[2]>=$IL)
		{
			last;
		}
	}
	my $temp = int(rand($i+1));
	my $ui;
	for($ui=0;$ui<scalar(@List);$ui++){
		if($List[$temp] eq $$ListStem[$ui]){
			last
		}
	}
	return($ui);
}

sub IniPhe()
{
	my($ListStem,$Teta)=@_;
	my($i);
	my($num);
	my($j);
	for($i=0;$i<@$ListStem;++$i)
	{
		$$Teta[$i][$i]=0;
		$num=0;
		for($j=0;$j<@$ListStem;++$j)
		{
			if(Consistance($$ListStem[$i],$$ListStem[$j])==1)
			{
				++$num;
			}
		}
		for($j=0;$j<@$ListStem;++$j)
		{
			if(Consistance($$ListStem[$i],$$ListStem[$j])==1)   # fagat sazegar barmidarad
#			if($i!=$j)
			{
				$$Teta[$i][$j]=2/$num;
			}
			else
			{
				$$Teta[$i][$j]=0;
			}
		}
	}
}


sub ExtractStems()
{
	my($ListStem,$Eta,$seqRNA)=@_;
	my($i);
	my($j);
	my($start);
	my($end);
	my($k);
	my($lent);
	my(@U);
	my(@W);
	for($i=0;$i<length($seqRNA);++$i)
	{
		for($j=length($seqRNA)-1;$j>=0;--$j)
		{

			if($j<=$i)
			{
				last;
			}
			if(complement(substr($seqRNA,$i,1),substr($seqRNA,$j,1))==1)
			{
				$start=$i;
				$end=$j;
				for($k=1;$k<length($seqRNA);++$k)
				{
					if(complement(substr($seqRNA,$i+$k,1),substr($seqRNA,$j-$k,1))==0)
					{
						if(abs($i+$k-1-($j-$k+1)+1)<3)
						{
							$k=0;
						}
						last;
					}
					if($j-$k<=$i+$k)
					{
						$k=0;
						last;
					}
				}
				$lent=$k;
				if($lent>=3)
				{
						push(@$ListStem,$start.",".$end.",".$lent);
				}
			}
		}
	}
	my($T);
	for($i=0;$i<@$ListStem;++$i)
	{
		$$Eta[$i][$i]=0;
		for($j=$i+1;$j<@$ListStem;++$j)
		{
			$$Eta[$i][$j]=0;
			$$Eta[$j][$i]=0;
			$T=Consistance($$ListStem[$i],$$ListStem[$j]);
			if($T==1)
			{
				@U=split(",",$$ListStem[$i]);
				@W=split(",",$$ListStem[$j]);
				$$Eta[$i][$j]=$W[2]/$U[2]*$W[2];
				$$Eta[$j][$i]=$U[2]/$W[2]*$U[2];
			}
		}
	}	
	return();
}

sub complement()
{
	my($A,$B)=@_;
	if(Com($A) eq $B)
	{
		return(1);
	}
	return(0);
}

sub Com()
{
	my($A)=@_;
	if($A eq "A")
	{
		return("U");
	}
	elsif($A eq "C")
	{
		return("G");
	}	
	elsif($A eq "G")
	{
		return("C");
	}
	elsif($A eq "U")
	{
		return("A");
	}
}

sub Consistance()
{
	my($L1,$L2)=@_;
	my(@A);
	my(@B);
	if($L1 eq $L2)
	{
		return(0);
	}
	@A=split(",",$L1);
	@B=split(",",$L2);
	if($A[0]+$A[2]<=$B[0] && $B[1]<=$A[1]-$A[2] )
	{
		return(1);
	}
	elsif($B[0]+$B[2]<=$A[0] && $A[1]<=$B[1]-$B[2] )
	{
		return(1);
	}
	elsif( $A[1]<$B[0] || $B[1]<$A[0])
	{
		return(1);
	}
	return(0);
}

sub MergeSort()
{
        my($A,$s,$e)=@_;
        if($s<$e)
        {
                MergeSort($A,$s,int(($s+$e)/2));
                MergeSort($A,int((($s+$e)/2)+1),$e);
                Merge($A,$s,int(($s+$e)/2),$e);
        }
}
sub Merge()
{
        my($A,$p,$q,$r)=@_;
        my($n1);
        my($n2);
        my(@L);
        my(@R);
        my($i);
        my($j);
        my($k);
        my(@T1);
        my(@T2);
        $n1=$q-$p+1;
        $n2=$r-$q;
        for($i=0;$i<$n1;++$i)
        {
                $L[$i]=$$A[$p+$i]

        }
        for($i=0;$i<$n2;++$i)
        {
                $R[$i]=$$A[$q+$i+1]
        }

        $i=0;
        $j=0;
        for($k=$p;$k<=$r;++$k)
        {
				if($i==$n1 || $j==$n2)
				{
					last;
				}  
				@T1=split(",",$L[$i]);
                @T2=split(",",$R[$j]);
                if($T1[@T1-1]>$T2[@T1-1])
                {
                        $$A[$k]=$L[$i];
                        ++$i;
                }
                elsif($T1[@T1-1]==$T2[@T1-1])
                {

                        if($T1[@T1-2]>$T2[@T1-2])
                        {
                                $$A[$k]=$L[$i];
                                ++$i;
                        }
                        else
                        {
                                $$A[$k]=$R[$j];
                                ++$j;
                        }
                }
                else
                {
                        $$A[$k]=$R[$j];
                        ++$j;
                }
        }
        for($p=$i;$p<$n1;++$p)
        {
			$$A[$k++]=$L[$p];
		}
        for($p=$j;$p<$n2;++$p)
        {
			$$A[$k++]=$R[$p];
		}

	
}

