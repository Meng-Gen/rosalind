use strict;
use warnings;

while (my $dna = <>) {
    chomp $dna;
    my @symbols = split //, $dna;
    my %count = (
        'A' => 0,
        'C' => 0,
        'G' => 0,
        'T' => 0,
    );
    for my $symbol (@symbols) {
        $count{$symbol}++;
    }
    print $count{'A'} . ' ' , 
          $count{'C'} . ' ' , 
          $count{'G'} . ' ' ,
          $count{'T'} . "\n";
}