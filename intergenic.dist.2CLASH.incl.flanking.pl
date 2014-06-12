#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;



	
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
R.D. Emes  Keele University 2008
USAGE : calculate intergenic distances and potential for clashes between ORF of all genes
Takes into account that both upstream and downstream genes need to fit their UTRs

-d path to database of PlasmoDB file
plasmoDB_data.RDE file of tab delimited [name chr start stop strand]
-t transcript - ORF distance to assume

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
";
###
# works on the assumpion that we have three genes left-current-right
# determines up/down stream by orinetation.
# ie if current + then upstream = current start - left end
# if current - then upstream = right start - current end 


my $database;
my $transcript_length_minus_ORF;


GetOptions(
	'd|database:s' => \$database,
	't|transcript:s' => \$transcript_length_minus_ORF,
          );

if( ! defined $database) {
print "\n\nWARNING: Cannot proceed without a valid database \n$usage\n\n"; exit;
}


open OUT, ">$database\.clashes.with.flanking.$transcript_length_minus_ORF";

my @file_data = get_file_data($database);#retrieve file data

my $nextline = 0;
my $linecounter = 0; 
my $line = "";

my $leftgene = "";
my $leftchr = "";
my $leftstart = "";
my $leftend= "";
my $leftstrand= "";

my $currgene = "";
my $currchr = "";
my $currstart = "";
my $currend= "";
my $currstrand= "";
my $currpeplength = "";
my $currcodinglegth = "";

my $rightgene = "";
my $rightchr = "";
my $rightstart = "";
my $rightend= "";
my $rightstrand= "";


my $check = 0;

my $inter_type_left = ""; # T:H, H:H etc
my $inter_type_right = ""; # T:H, H:H etc
#~ my $check = 0;

print OUT "Left Gene\tLeft Gene \(Start_End_strand\)\tInter dist left\tInter type left\tCurrent Gene\tCurrent Gene \(Start_End_strand\)\tInter dist right\tInter type right\tRight Gene\tRight Gene \(Start_End_strand\)\t";
my $i = 100;
while ($i > -1)
{
print OUT "$i\t";
$i--;
}
print OUT "\n";



foreach $line(@file_data)
{
	chomp $line;
	my $previousline = $file_data[$linecounter-1];
	my $currentline = $file_data[$linecounter];
	my $nextline = $file_data[$linecounter+1];
	
	chomp $currentline;
	if (defined $nextline && $previousline)
		{
		chomp $nextline;
		my @nextdata = split  '\t', $nextline;
		$rightgene = $nextdata[0];
		$rightchr = $nextdata[1];
		$rightstart = $nextdata[2];
		$rightend = $nextdata[3];
		$rightstrand = $nextdata[4];	
		
		chomp $previousline;
		my @previousdata = split  '\t', $previousline;
		$leftgene = $previousdata[0];
		$leftchr = $previousdata[1];
		$leftstart = $previousdata[2];
		$leftend = $previousdata[3];
		$leftstrand = $previousdata[4];
		
	
		my @currentdata = split  '\t', $currentline;
		$currgene = $currentdata[0];
		$currchr = $currentdata[1];
		$currstart = $currentdata[2];
		$currend = $currentdata[3];
		$currstrand = $currentdata[4];
	

		
		
		if ($currchr eq $rightchr && $currchr eq $leftchr)
			{
			my $interdistA = ($currstart - $leftend);
			my $interdistB = ($rightstart - $currend);
			
			if ($leftstrand eq "+" && $currstrand eq "-")
				{
				$inter_type_left = "T2T";
				}
			elsif ($leftstrand eq "+" && $currstrand eq "+")
				{
				$inter_type_left = "T2H";
				}
			elsif ($leftstrand eq "-" && $currstrand eq "-")
				{
				$inter_type_left = "H2T";
				}
			elsif ($leftstrand eq "-" && $currstrand eq "+")
				{
				$inter_type_left = "H2H";
				}
			
			if ($currstrand eq "-" && $rightstrand eq "+")
				{
				$inter_type_right = "H2H";
				}
			elsif ($currstrand eq "+" && $rightstrand eq "+")
				{
				$inter_type_right = "T2H";
				}
			elsif ($currstrand eq "-" && $rightstrand eq "-")
				{
				$inter_type_right = "H2T";
				}
			elsif ($currstrand eq "+" && $rightstrand eq "-")
				{
				$inter_type_right = "T2T";
				}		
			
				my $upstream_dist = 0;
				my $downstream_dist = 0;			
			
				my $percent_count = 100;
				print OUT "$leftgene\t$leftstart $leftend $leftstrand\t$interdistA\t$inter_type_left\t$currgene\t$currstart $currend $currstrand\t$interdistB\t$inter_type_right\t$rightgene\t$rightstart $rightend $rightstrand\t";
				while ($percent_count > -1)
				{
				my $up_percent = $percent_count;
				my $down_percent = (100 - $up_percent);
				my $up_distance = sprintf("%.0f", (($transcript_length_minus_ORF/100)*$up_percent));
				my $down_distance = sprintf("%.0f", (($transcript_length_minus_ORF/100)*$down_percent));

			
			
			
			if ($currstrand eq "+") 
				{
				$upstream_dist = ($currstart - $leftend); # distance 
				$downstream_dist = ($rightstart- $currend);
				
				if ($inter_type_left eq "T2H") 
					{
					if ($upstream_dist >= ($down_distance + $up_distance))
						{
						if ($inter_type_right eq "T2H")
							{
							if ($downstream_dist >= ($down_distance + $up_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						if ($inter_type_right eq "T2T")
							{
							if ($downstream_dist >= ($down_distance + $down_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						}
					else {
						print OUT "1\t";
						}
					}
				if ($inter_type_left eq "H2H") 
					{
					if ($upstream_dist >= ($up_distance + $up_distance))
						{
						if ($inter_type_right eq "T2H")
							{
							if ($downstream_dist >= ($down_distance + $up_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						if ($inter_type_right eq "T2T")
							{
							if ($downstream_dist >= ($down_distance + $down_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}	
							}
						}
					else {
						print OUT "1\t";
						}
					}	
				}
			
			elsif ($currstrand eq "-") 
				{
				$upstream_dist = ($rightstart - $currend);
				$downstream_dist = ($currstart- $leftend);
				
				if ($inter_type_left eq "T2T")
					{
					if ($downstream_dist >= ($down_distance + $down_distance))
						{
						if ($inter_type_right eq "H2H")
							{
							if ($upstream_dist >= ($up_distance + $up_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						if ($inter_type_right eq "H2T")
							{
							if ($upstream_dist >= ($up_distance + $down_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}	
							}
						}
					else {
						print OUT "1\t";
						}
					}
				if ($inter_type_left eq "H2T")
					{
					if ($downstream_dist >= ($up_distance + $down_distance))
						{
						if ($inter_type_right eq "H2H")
							{
							if ($upstream_dist >= ($up_distance + $up_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						if ($inter_type_right eq "H2T")
							{
							if ($upstream_dist >= ($up_distance + $down_distance)) {print OUT "0\t";}
							else {print OUT "1\t";}
							}
						}
					else {
						print OUT "1\t";
						}
					}
				
				}
			
			

			$percent_count --;
			}
			
			
			
	
	#~ $check = 0;

print OUT "\n";
	}

}
$linecounter++;
#~ print "$leftgene\t$leftstart $leftend $leftstrand\t$inter_type_left\t$currgene\t$currstart $currend $currstrand\t$inter_type_right\t$rightgene\t$rightstart $rightend $rightstrand\t$upstream_dist\t$downstream_dist\t$max_percent_up\t$max_percent_down\t$clashes\n";

}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_file_data
# A subroutine to get data from a file given its filename
sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}
