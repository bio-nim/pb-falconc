Test two: No sequences should be rotated with an empty whitelist

  $ touch empty.whitelist
  $ falconc circ-orient --window 20 -s 2 -i small.fasta -w empty.whitelist -o small.out2.fasta > junk
  $ cat small.out2.fasta
  >first shifted_by_bp:-0/180
  AGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAACACAAGAAGCTATACAATTGTCTTTTCAAGGCTCGTCTGTG
  CTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGACCCATTGGGAAGACTCTTCT
  ATAAAAACAGATGAAATTTT
  >second shifted_by_bp:-0/180
  AGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAACACAAGAAGCTATACAATTGTCTTTTCAAGGCTCGTCTGTG
  CTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGACCCATTGGGAAGACTCTTCT
  ATAAAAACAGATGAAATTTT
  >third shifted_by_bp:-0/180
  AGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAACACAAGAAGCTATACAATTGTCTTTTCAAGGCTCGTCTGTG
  CTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGACCCATTGGGAAGACTCTTCT
  ATAAAAACAGATGAAATTTT
  >fourth shifted_by_bp:-0/180
  AGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAACACAAGAAGCTATACAATTGTCTTTTCAAGGCTCGTCTGTG
  CTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGACCCATTGGGAAGACTCTTCT
  ATAAAAACAGATGAAATTTT