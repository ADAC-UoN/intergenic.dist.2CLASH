#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;



	
my $usage = "
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
R.D. Emes  Keele University 2008
USAGE : calculate intergenic distances and potential for clashes between ORF of all genes

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

open OUT, ">$database\.clash\.output\.$transcript_length_minus_ORF";


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

print OUT "Left Gene\tLeft Gene \(Start_End_strand\)\tinter_type\(Left-Current\)\tCurrent Gene\tCurrent Gene \(Start_End_strand\)\tinter_type\(Current-Right\)\tRight Gene\tRight Gene \(Start_End_strand\)\tUpstream dist(bp)\tDownstream dist(bp)\tMax Percent available (upstream)\tMax Percent available (downstream)\t";
my $i = 100;
while ($i > -1)
{
print  OUT "$i\t";
$i--;
}
print  OUT "\n";

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
	#~ $currpeplength = $currentdata[10];
	#~ $currcodinglegth = ($currpeplength*3);
	
	

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
		
		
		
		
		if ($currstrand eq "+") 
			{
			$upstream_dist = ($currstart - $leftend); # distance 
			$downstream_dist = ($rightstart- $currend);
			}
		elsif ($currstrand eq "-") 
			{
			$upstream_dist = ($rightstart - $currend);
			$downstream_dist = ($currstart- $leftend);
			}
		
		
		my $max_percent_left = sprintf("%.0f", (($interdistA/$transcript_length_minus_ORF)*100));
		my $max_percent_right = sprintf("%.0f", (($interdistB/$transcript_length_minus_ORF)*100));
		
		my $max_percent_up = sprintf("%.0f", (($upstream_dist/$transcript_length_minus_ORF)*100));
		my $max_percent_down = sprintf("%.0f", (($downstream_dist/$transcript_length_minus_ORF)*100));
		
		
		if ($max_percent_left > 100){$max_percent_left = 100;}
		if ($max_percent_right > 100){$max_percent_right = 100;}
		
		if ($max_percent_up > 100){$max_percent_up = 100;}
		if ($max_percent_down > 100){$max_percent_down =100;}		
		
		my $clashes = "";
		my $print_1 = "1\t";
		my $print_0 = "0\t";
		my $left_1 = "";
		my $right_1 = "";
		my $mid_0 = "";
		
		
		if ($max_percent_up == 100 && $max_percent_down == 100) # no clash therefore print 101 * 0.
			{
			$clashes = $print_0 x 101;
			}
		elsif (($max_percent_up + $max_percent_down) < 100) # will always clash therefore print 101 * 1.
			{
			$clashes = $print_1 x 101;
			}
		elsif ($max_percent_up == 50 && $max_percent_down == 50) # oddity 
			{
			$clashes = ($print_1 x 50).($print_0 x 1).($print_1 x 50);
			}
			
		elsif ($max_percent_up < $max_percent_down)  # checked 
			{
			$left_1 = (100 - $max_percent_up);					# calculate how much fits upstream  
			$right_1 = (100 - $max_percent_down); 				# calculate what fits downstream 
			$mid_0 = (101 - ($left_1 + $right_1)); 				# calculate how much has to fit in the middle 
      			$clashes = ($print_1 x ($left_1)).($print_0 x ($mid_0)).($print_1 x ($right_1));
			} 
		 
		elsif ($max_percent_up > $max_percent_down) # checked                    
			{
			$left_1 = (100 - $max_percent_up); 
			$right_1 = (100 - $max_percent_down);
			$mid_0 = (101 - ($left_1 + $right_1)); 
			$clashes = ($print_1 x ($left_1)).($print_0 x ($mid_0)).($print_1 x ($right_1));
			}
			
		elsif ($max_percent_up == $max_percent_down) # 
			{
			$left_1 = (100 - $max_percent_up); 
			$right_1 = (100 - $max_percent_down);
			$mid_0 = (101 - ($left_1 + $right_1)); 
			$clashes = ($print_1 x ($left_1)).($print_0 x ($mid_0)).($print_1 x ($right_1));
			}
		
		
		print  OUT "$leftgene\t$leftstart $leftend $leftstrand\t$inter_type_left\t$currgene\t$currstart $currend $currstrand\t$inter_type_right\t$rightgene\t$rightstart $rightend $rightstrand\t$upstream_dist\t$downstream_dist\t$max_percent_up\t$max_percent_down\t$clashes\n";
		}
		
		
		
$linecounter++;
#~ $check = 0;
}
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
