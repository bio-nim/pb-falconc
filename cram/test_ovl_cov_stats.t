Test on a fofn of overlap piles.
  $ ls -1 ${PROJECT_DIR}/tests/data/chimera/*.ovlp > in.fofn
  > ${PROJECT_DIR}/src/falconc ovl-cov-stats --in-fn in.fofn 2>/dev/null
  {
      "mean": 18.125,
      "std": 9.873417594733853,
      "median": 20.5,
      "min": 3,
      "max": 34,
      "n": 16
  }

Test on a non-fofn input.
  $ ${PROJECT_DIR}/src/falconc ovl-cov-stats --in-fn ${PROJECT_DIR}/tests/data/chimera/tp_15210.ovlp 2>/dev/null
  {
      "mean": 19.0,
      "std": 11.0,
      "median": 19.0,
      "min": 8,
      "max": 30,
      "n": 2
  }

