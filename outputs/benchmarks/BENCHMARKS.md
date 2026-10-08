# Benchmark results: analysistools (old) vs fastmsna (new)

Generated 2026-10-08 01:52 on Windows, R 4.5.1, 8 cores. Each case in a fresh R session.

`row-indicators/s` = dataset rows x LOA rows / seconds. Memory = peak R heap during the analysis.

| part | rows | analysis | LOA rows | old (s) | new (s) | speed-up | old peak MB | new peak MB | old row-ind/s | new row-ind/s | results equal | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A_full_loa |   1,000 | all types | 154 | 157.82 | 0.92 | 171.5x | 325.6 | 93 |    976 |   167,391 | TRUE |  |
| A_full_loa |   5,000 | all types | 154 | 152.41 | 1.05 | 145.2x | 253 | 92 |  5,052 |   733,333 | TRUE |  |
| A_full_loa |  20,000 | all types | 154 | 211.76 | 1.97 | 107.5x | 366.2 | 179 | 14,545 | 1,563,452 | TRUE |  |
| A_full_loa | 100,000 | all types | 154 | 644.87 | 4.84 | 133.2x | 1176.5 | 531 | 23,881 | 3,181,818 | TRUE |  |
| B_per_function |  20,000 | prop_select_one | 13 |  56.60 | 0.33 | 171.5x | 357.6 | 85 |  4,594 |   787,879 | TRUE |  |
| B_per_function |  20,000 | prop_select_multiple | 4 |  16.03 | 0.11 | 145.7x | 300.8 | 68 |  4,991 |   727,273 | TRUE |  |
| B_per_function |  20,000 | mean | 6 |   6.71 | 0.28 | 24x | 366.1 | 40 | 17,884 |   428,571 | TRUE |  |
| B_per_function |  20,000 | median | 6 |   8.85 | 0.29 | 30.5x | 295.8 | 101 | 13,559 |   413,793 | TRUE |  |
| B_per_function |  20,000 | ratio | 2 |   2.26 | 0.05 | 45.2x | 303.8 | 37 | 17,699 |   800,000 | TRUE |  |
| B_per_function | 100,000 | prop_select_one | 13 | 121.87 | 0.81 | 150.5x | 1143.5 | 363 | 10,667 | 1,604,938 | TRUE |  |
| B_per_function | 100,000 | prop_select_multiple | 4 |  38.36 | 0.22 | 174.4x | 1145.2 | 216 | 10,428 | 1,818,182 | TRUE |  |
| B_per_function | 100,000 | mean | 6 |  17.84 | 0.11 | 162.2x | 1143.4 | 123 | 33,632 | 5,454,545 | TRUE |  |
| B_per_function | 100,000 | median | 6 |  28.00 | 0.40 | 70x | 1143.4 | 278 | 21,429 | 1,500,000 | TRUE |  |
| B_per_function | 100,000 | ratio | 2 |   7.00 | 0.10 | 70x | 921 | 101 | 28,571 | 2,000,000 | TRUE |  |
