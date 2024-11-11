# USE $Int_low=-1 FOR Z-SCORE FILES AND $Int_low=0 FOR NZ FILES

use strict;
use warnings;

my $filname_inp=$ARGV[0];
my $motif=$ARGV[1];
my $len_seq=$ARGV[2];
my $colm=$ARGV[3];

my $filname_res='perl_'.$motif.'_'.$filname_inp;
my $Int_low=0;

my $len_motif=length($motif);
my $mot_rev=reverse $motif;
$mot_rev=~tr/ACGT/TGCA/;

my %Int=();
open(fl_inp,"$filname_inp");
while (<fl_inp>)
{
chomp;
my @elements=split(/\s+/,$_);
$Int{$elements[0]}=$elements[1+$colm];
$Int{$elements[1]}=$elements[1+$colm];
}
close(fl_inp);


open(fl_op,">$filname_res");

my %rev=();
my %no_mis=();
my %mis_pos=();
my %mot_pos=();

my $mism_now=0;
my $Int_now=$Int_low;

(my $t1,my $mot_posx)=permutate_flanks2($motif,$len_seq);
my @t2=@{$t1};
my %mot_posp=%{$mot_posx};


foreach my $k(@t2)
{
	my $rev_v=reverse $k;
	$rev_v=~tr/ACGT/TGCA/;
	unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
	{
		$rev{$k}=$rev_v;
		$no_mis{$k}=$mism_now;
		@{$mis_pos{$k}}=(0) x $len_motif;
		$mot_pos{$k}=$mot_posp{$k}+2;
		if (exists $Int{$k})
		{$Int_now=$Int{$k};}
		else {$Int_now=$Int_low;}
		print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$k}})."\t".$mot_pos{$k}."\t".$no_mis{$k}."\n";	}
}
#print join("\n",keys(%rev));


$mism_now=1;
# MISMATCH 1
# including case for 1 nucleotide lost in ends

	my $mism_motif=substr($motif,1,$len_motif-1);
	#print "$mism_motif\n";
	
	( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
	 @t2=@{$t1};
	 %mot_posp=%{$mot_posx};

	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=(1, (0) x ($len_motif-1));
			my $x=$k; my $y=$rev_v;my $mot_pos_temp=$mot_posp{$k}-1;
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx}; $mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
			}
	}
	
	$mism_motif=substr($motif,0,$len_motif-1);
	#print "$mism_motif\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
	@t2=@{$t1};
	%mot_posp=%{$mot_posx};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=((0) x ($len_motif-1),1);
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k};
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx}; $mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
			}
	}
	
#print fl_op "abc\n";

my @acgt=('A','C','G','T');
for(my $i=0;$i<$len_motif;$i++)
{
	foreach my $mism(@acgt)
	{
	$mism_motif=substr($motif,0,$i).$mism.substr($motif,$i+1,$len_motif-$i-1);
	#print "$mism_motif\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
	@t2=@{$t1};
	%mot_posp=%{$mot_posx};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=((0) x $i,1,(0) x ($len_motif-$i-1));
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k};
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
		}
	}
	
	}

}


$mism_now=2;
#print "\n\n";
# MISMATCH 2
# including case for 1mismatch in motif and 1 nucleotide lost in ends
$mism_motif=substr($motif,1,$len_motif-1);
for(my $i=0;$i<$len_motif-1;$i++)
{
	foreach my $mism(@acgt)
	{
	my $mism_motifx=substr($mism_motif,0,$i).$mism.substr($mism_motif,$i+1,$len_motif-$i-1-1);
	#print "$mism_motifx\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motifx,$len_seq);
	@t2=@{$t1};
	%mot_posp=%{$mot_posx};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=(1,(0) x $i,1,(0) x ($len_motif-$i-1-1));
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k}-1;
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
		}
	}
	
	}

}
	#print fl_op "abc2\n";

$mism_motif=substr($motif,0,$len_motif-1);
for(my $i=0;$i<$len_motif-1;$i++)
{
	foreach my $mism(@acgt)
	{
	my $mism_motifx=substr($mism_motif,0,$i).$mism.substr($mism_motif,$i+1,$len_motif-$i-1-1);
	#print "$mism_motifx\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motifx,$len_seq);
	%mot_posp=%{$mot_posx};
	@t2=@{$t1};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=((0) x $i,1,(0) x ($len_motif-$i-1-1),1);
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k};
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
		}
	}
	
	}

}
	#print fl_op "abc3\n";

# including case for 2 nucleotide lost in ends

	$mism_motif=substr($motif,2,$len_motif-2);
	#print "$mism_motif\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
	@t2=@{$t1};
	%mot_posp=%{$mot_posx};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=(1,1, (0) x ($len_motif-2));
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k}-2;
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
		}
	}
	
	$mism_motif=substr($motif,0,$len_motif-2);
	#print "$mism_motif\n";
	( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
	@t2=@{$t1};
	%mot_posp=%{$mot_posx};
	foreach my $k(@t2)
	{
		my $rev_v=reverse $k;
		$rev_v=~tr/ACGT/TGCA/;
		unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
		{
			my @tem1=( (0) x ($len_motif-2),1,1);
			my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k};
			(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
			if ($frwd_rev==1)
			{
				$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
			}
			$rev{$x}=$y;
			$no_mis{$x}=$mism_now;
			@{$mis_pos{$x}}=@tem1;
			$mot_pos{$x}=$mot_pos_temp+2;
			if (exists ($Int{$x}))
			{$Int_now=$Int{$x};}
			else {$Int_now=$Int_low;}
			print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
		}
	}
	
	#print fl_op "abc4\n";
my $two_nucleotidex=all_nmers_rna(2);
my @two_nucleotide=@{$two_nucleotidex};

for(my $i=0;$i<$len_motif;$i++)
{
	for(my $j=$i+1;$j<$len_motif;$j++)
	{
		foreach my $mism(@two_nucleotide)
		{
		$mism_motif=substr($motif,0,$i).substr($mism,0,1).substr($motif,$i+1,$j-$i-1).substr($mism,1,1).substr($motif,$j+1,$len_motif-$j-1);
		#print "$i\t$j\t$mism_motif\n";
		( $t1, $mot_posx)=permutate_flanks2($mism_motif,$len_seq);
		@t2=@{$t1};
		%mot_posp=%{$mot_posx};
		foreach my $k(@t2)
		{
			my $rev_v=reverse $k;
			$rev_v=~tr/ACGT/TGCA/;
			unless ( defined($rev{$k}) || defined($rev{$rev_v}) ) 
			{
				my @tem1=((0) x $i,1,(0) x ($j-$i-1),1,(0) x ($len_motif-$j-1));
				my $x=$k; my $y=$rev_v; my $mot_pos_temp=$mot_posp{$k};
				(my $mis_pos_revx, my $frwd_rev)=mism($k,$mot_pos_temp,$mot_rev,$len_motif,$motif,$mism_now,\@tem1);
				if ($frwd_rev==1)
				{
					$y=$k; $x=$rev_v; @tem1=@{$mis_pos_revx};$mot_pos_temp=$len_seq-$mot_pos_temp-$len_motif+2;
				}
				$rev{$x}=$y;
				$no_mis{$x}=$mism_now;
				@{$mis_pos{$x}}=@tem1;
				$mot_pos{$x}=$mot_pos_temp+2;
				if (exists ($Int{$x}))
				{$Int_now=$Int{$x};}
				else {$Int_now=$Int_low;}
				print fl_op $k."\t".$Int_now."\t".join("\t",@{$mis_pos{$x}})."\t".$mot_pos{$x}."\t".$no_mis{$x}."\n";
			}
		}
		
		}

	}
}


#print fl_op join("\n",keys(%rev));
close(fl_op);










sub mism
{
my $k=$_[0];
my $mot_posp1=$_[1];
my $mot_rev=$_[2];
my $len_motif=$_[3];
my $motif=$_[4];
my $mism_now=$_[5];
my @mis_pos=@{$_[6]};

my $no_mis_rev=0;
my @mis_pos_rev=();

my $frwd_rev=0; # 0==frwd 1==rev
my $mis_pos_sum=0;
my $mis_pos_rev_sum=0;


$k='XX'.$k.'XX';
$k=substr($k,$mot_posp1+1,$len_motif);

for (my $i=0;$i<$len_motif;$i++)
{
	if (substr($mot_rev,$i,1) eq substr($k,$i,1))
	{
		push @mis_pos_rev,0;
		$mis_pos_rev_sum=2*$mis_pos_rev_sum;
	}
	else
	{
		push @mis_pos_rev,1; 
		$no_mis_rev++;
		$mis_pos_rev_sum=2*$mis_pos_rev_sum+1;
	}
	$mis_pos_sum=2*$mis_pos_sum+$mis_pos[$i];
}

if (($mis_pos_sum > $mis_pos_rev_sum) && ($no_mis_rev==$mism_now))
{
	$frwd_rev=1;
}
return (\@mis_pos_rev, $frwd_rev);
}



sub permutate_flanks2
{
my $motif=$_[0];
my $tot_len=$_[1];
my $len_mot=length($motif);
my @op=();

my $all_nm=all_nmers_rna($tot_len-$len_mot);
my @all_nm_r=@{$all_nm};
my %mot_posx;

for(my $i=0;$i<=$tot_len-$len_mot;$i++)
{
	foreach my $v(@all_nm_r)
	{
		my $x=substr($v,0,$i).$motif.substr($v,$i,$tot_len-$len_mot-$i);
		push @op,$x;
		$mot_posx{$x}=$i+1;
	}
}

return (\@op,\%mot_posx);

}

sub all_nmers_rna
{
my $len=$_[0];
my $mut=0;
for(my $i=0;$i<$len;$i++)
{
$mut=$mut*10+1;
}
my @arr=();
my %x=();
for(my $i=0;$i<4**$len;$i++)
{
my $a=$i;
my $v=0;
my $m=1;
while ($a>0)
{
$v=$v+$m*($a%4);
$m=$m*10;
$a=($a-($a%4))/4;
}
$v=$v+$mut;
$v=~tr/1234/ACGT/;
push @arr,$v;
}
return (\@arr);;
}
