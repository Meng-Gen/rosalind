dna = STDIN.read.chomp
puts [dna.count('A'), dna.count('C'), dna.count('G'), dna.count('T')].join ' '
