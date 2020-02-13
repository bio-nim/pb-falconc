Test circ-orient

  $ rm -f circ-rotate.fq
  $ falconc circ-orient --window 500 -o circ-rotate.fq -i /pbi/dept/secondary/siv/testdata/hgap/rotate/hg.fq > /dev/null 2>&1
  $ diff circ-rotate.fq /pbi/dept/secondary/siv/testdata/hgap/rotate/hg-win500.fq
  $ rm -f circ-rotate.fq
