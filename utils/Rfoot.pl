use Getopt::Std;

### get arguments ###
my %args; 
getopt("itowcrs",\%args);
my $samfile=$args{i}; 
my $transfile=$args{t}; 
my $outfile=$args{o}; 
my $win=$args{w};
my $readcoff=$args{r};

($win=60) if (!$win);
($readcoff=10) if (!$readcoff);

############
my %dist;
for (my $i=15; $i<=50; $i++) {
	$dist{$i}=int($i/2);
}

my %read;
open (IN, "$samfile");
while (<IN>) {
	chomp;
	if ($_ !~ /^@/) {
		my @s1=split /\t/, $_;
		my @s2=split /\D+/, $s1[5];
		my @s3=split /\d+/, $s1[5];
		my $len=0;
		for (my $i=0; $i<$#s3; $i++) {
			if ($s3[$i+1] eq 'M') {
				$len+=$s2[$i];
			}
		}
		if (exists ($dist{$len})) {
			my $asite=$dist{$len};
			my $loc=$s1[3];
			my $ra=0;
			my $ind=1;
			if ($s1[1]==0) {
				for (my $i=0; $i<$#s3 && $ind==1; $i++) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc+=$asite;
							$ind=0;
						} else {
							$loc+=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc+=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc-=$s2[$i];
					}
				}
				my $k=$s1[2].":"."+".":".($loc-1);
				$read{$k}++;
			} else {
				for (my $i=0; $i<$#s3; $i++) {
					$loc+=$s2[$i];
				}
				$loc--;
				for (my $i=($#s3-1); $i>=0 && $ind==1; $i--) {
					if ($s3[$i+1] eq 'M') {
						if ($s2[$i] >= ($asite+1)) {
							$loc-=$asite;
							$ind=0;
						} else {
							$loc-=$s2[$i];
							$asite-=$s2[$i];
						}
					} elsif ($s3[$i+1] eq 'N') {
						$loc-=$s2[$i];
					} elsif ($s3[$i+1] eq 'D') {
						$loc+=$s2[$i];
					}
				}
				my $k=$s1[2].":"."-".":".($loc-1);
				$read{$k}++;
			}
		}
	}
}
close IN;

#################
my $ds=3;
open (AN, "$transfile");
open (OUT, ">$outfile");
print OUT "geneID"."\t"."transcriptID"."\t"."chrom"."\t"."strand"."\t"."start"."\t"."end"."\t"."length"."\t"."len.map"."\t"."read.num"."\t"."entropy"."\t"."MAXentropy"."\t"."pme"."\n";
while (<AN>) {
	chomp;
	my @s1=split /\t/, $_;
	my @s2=split /,/, $s1[9];
	my @s3=split /,/, $s1[10];
	my $k=$s1[0].":".$s1[1].":".$s1[2].":".$s1[3];
	my @tog;
	for (my $i=0; $i<=$#s2; $i++) {
		for (my $j=$s2[$i]; $j<$s3[$i]; $j++) {
			push @tog, $j;
		}
	}
	for (my $sam=0; $sam<=($#tog-$win+1); $sam++) {
	my @post;
	for (my $j=0; $j<$win; $j++) {
		push @post, $tog[$sam+$j];
	}
	my $len1=0;
	my $len2=0;
	my $tot=0;
	my @val;
	if ($s1[3] eq '-') {
		@post=reverse(@post);
	}
	for (my $m=0; $m<=$#post; $m++) {
		$len1++;
		my $k=$s1[2].":".$s1[3].":".$post[$m]; 
		if (! exists ($map{$k})) {
			$len2++;
		}
		if (exists ($read{$k})) {
			$tot+=$read{$k};
			$val[int(($len1-1)/$ds)]+=$read{$k};
		}
	}
	if ($tot >= $readcoff) {
		my $ent=0;
		my $ten=0;
		my $a=int(($len2+2)/$ds);
		my $b=$tot;
		my $t1=int(($a+$b-1)/$b);
		my @val2;
		for ($i=0; $i<=$#val; $i++) {
			if ($val[$i] > 0) {
				$val2[int($i/$t1)]+=$val[$i];
			}
		}
		for ($i=0; $i<=$#val2; $i++) {
			if ($val2[$i] > 0) {
				my $p=$val2[$i]/($tot);
				$ent+=$p*log(1/$p);
			}
		}
		my $t2=int($a/$t1);
		my $d1=int($b/$t2);
		my $d2=$b%$t2;
		my @va;
		for (my $i=0; $i<$t2; $i++) {
			$va[$i]=$d1;
		}
		for (my $j=0; $j<$d2; $j++) {
			$va[$j]++;
		}
		for (my $i=0; $i<=$#va; $i++) {
			if ($va[$i] > 0) {
				$p=$va[$i]/($tot);
				$ten+=$p*log(1/$p);
			}
		}
		my $per;
		if ($ten == 0) {
			$per=1;
		} else {
			$per=$ent/$ten;
		}
		if ($tot > 0) {
			my $out=$s1[0]."\t".$s1[1]."\t".$s1[2]."\t".$s1[3]."\t".$post[0]."\t".$post[$#post];
			$out.="\t".$len1."\t".$len2."\t".$tot;
			$out.="\t".sprintf("%.3f", $ent);
			$out.="\t".sprintf("%.3f", $ten);
			$out.="\t".sprintf("%.3f", $per);
			print OUT $out."\n";
		}
	}
	}
}
close AN;
close OUT;


