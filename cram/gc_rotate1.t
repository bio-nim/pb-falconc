Test one: All sequences should be rotated without a whitelist

  $ falconc circ-orient --window 20 -s 2 -i small.fasta -o small.out.fasta > junk
  $ cat small.out.fasta
  >first shifted_by_bp:-60/180
  CTTTTCAAGGCTCGTCTGTGCTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGA
  CCCATTGGGAAGACTCTTCTATAAAAACAGATGAAATTTTAGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAAC
  ACAAGAAGCTATACAATTGT
  >second shifted_by_bp:-60/180
  CTTTTCAAGGCTCGTCTGTGCTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGA
  CCCATTGGGAAGACTCTTCTATAAAAACAGATGAAATTTTAGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAAC
  ACAAGAAGCTATACAATTGT
  >third shifted_by_bp:-60/180
  CTTTTCAAGGCTCGTCTGTGCTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGA
  CCCATTGGGAAGACTCTTCTATAAAAACAGATGAAATTTTAGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAAC
  ACAAGAAGCTATACAATTGT
  >fourth shifted_by_bp:-60/180
  CTTTTCAAGGCTCGTCTGTGCTTTAGATAATTATGAAAAGGAAACAGTTGGGGTGTGCTTGCTTATAGGAGTCTATATGA
  CCCATTGGGAAGACTCTTCTATAAAAACAGATGAAATTTTAGGAATCCTCTTTGTGGGCACTTAAAGCCCATGAAGCAAC
  ACAAGAAGCTATACAATTGT
