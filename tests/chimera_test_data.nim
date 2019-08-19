## These are p-read overlaps from RAPTOR in the I.P.A./M.A.A. pipeline.
## These specific a-reads are either truly chimeric or false positives.

import streams

let tp_5600* = newStringStream(

                """5600 10005 -639 90.8333 0 1944 3264 9722 1 4349 5671 8935 u
5600 10219 -6264 99.9363 0 3437 9722 9722 1 245 6528 7272 u
5600 10419 -639 90.8333 0 1944 3264 9722 0 1877 3199 8492 u
5600 10793 -3406 99.739 0 0 3451 9722 0 5225 8673 9876 u
5600 11037 -639 90.8333 0 1944 3264 9722 0 10666 11988 18326 u
5600 11289 -639 90.8333 0 1944 3264 9722 1 308 1630 9682 u
5600 11356 -1922 99.5459 0 7733 9722 9722 1 4317 6299 6299 3
5600 11981 -607 90.6818 0 1944 3264 9722 0 6309 7630 9949 u
5600 12109 -1296 99.6918 0 2151 3451 9722 1 6411 7709 7709 u
5600 12196 -634 90.6818 0 1944 3264 9722 0 8779 10099 10912 u
5600 12653 -639 90.8333 0 1944 3264 9722 0 724 2046 9689 u
5600 12703 -639 90.8333 0 1944 3264 9722 1 9643 10965 14581 u
5600 12728 -3143 99.8428 0 0 3180 9722 0 16123 19306 19306 5
5600 1299 -2220 99.8655 0 3437 5670 9722 1 0 2230 10218 u
5600 1323 -639 90.8333 0 1944 3264 9722 0 4281 5603 12949 u
5600 13408 -5586 99.9465 0 3437 9044 9722 0 3188 8794 8794 u
5600 13619 -6273 99.9682 0 3437 9722 9722 0 1088 7373 8251 u
5600 1376 -635 90.7576 0 1944 3264 9722 0 2950 4271 9739 u
5600 13947 -7184 99.9307 0 95 7310 9722 1 0 7214 7214 contains
5600 14113 -638 90.7576 0 1944 3264 9722 0 8728 10049 11881 u
5600 14230 -2937 100 0 3437 6375 9722 1 0 2938 7896 u
5600 14367 -639 90.8333 0 1944 3264 9722 1 3052 4374 12491 u
5600 14398 -5889 99.899 0 3779 9722 9722 0 0 5939 9198 3
5600 14399 -639 90.8333 0 1944 3264 9722 1 1543 2865 13603 u
5600 1476 -6273 99.9841 0 3437 9722 9722 1 1450 7734 9772 u
5600 14907 -3407 99.8551 0 0 3451 9722 0 5466 8918 22715 u
5600 14939 -635 90.7576 0 1944 3264 9722 0 2828 4149 11232 u
5600 1519 -639 90.8333 0 1944 3264 9722 1 10835 12157 15591 u
5600 1523 -1030 99.9008 0 2442 3451 9722 0 0 1008 10022 u
5600 15580 -3420 99.7972 0 0 3451 9722 1 12276 15728 19603 u
5600 15817 -3417 99.7391 0 0 3451 9722 0 3681 7131 15521 u
5600 15934 -3413 99.7972 0 0 3451 9722 0 4871 8323 13293 u
5600 16476 -3424 99.9131 0 0 3451 9722 0 1589 5043 7750 u
5600 16601 -3413 99.6249 0 0 3473 9722 0 8025 11491 11491 5
5600 16867 -1192 99.507 0 8501 9722 9722 1 6894 8111 8111 3
5600 17610 -3637 99.7848 0 0 3717 9722 0 3360 7077 7077 5
5600 17804 -3419 99.7972 0 0 3451 9722 0 3854 7306 10514 u
5600 19094 -3187 99.8463 0 0 3253 9722 1 0 3256 8486 5
5600 19298 -2608 99.8491 0 0 2650 9722 0 10430 13082 13082 5
5600 19373 -3401 99.7102 0 0 3451 9722 0 197 3648 7429 u
5600 19461 -3970 99.9498 0 3437 7422 9722 1 0 3983 8026 u
5600 19504 -3244 99.4516 0 169 3451 9722 0 0 3294 6750 u
5600 19516 -5034 99.4327 0 4413 9722 9722 1 19687 24975 24975 3
5600 19607 -2298 99.8704 0 3437 5753 9722 0 8029 10344 10344 u
5600 19683 -634 90.6818 0 1944 3264 9722 0 5601 6921 8914 u
5600 3497 -2128 99.8131 0 1311 3451 9722 1 9391 11531 11531 u
5600 3631 -3410 99.7971 0 0 3451 9722 1 7898 11348 13248 u
5600 4486 -3406 99.6813 0 0 3451 9722 0 2833 6284 9350 u
5600 4600 -6247 99.9045 0 3437 9722 9722 0 3065 9346 14335 u
5600 462 -1085 99.4652 0 8598 9722 9722 1 8194 9316 9316 3
5600 5161 -4605 99.7662 0 5012 9722 9722 1 6525 11230 11230 3
5600 5175 -3432 99.7971 0 0 3451 9722 0 6851 10301 13989 u
5600 5176 -3700 99.8137 0 5960 9722 9722 1 9507 13265 13265 3
5600 5347 -1363 100 0 2106 3451 9722 0 0 1345 7626 u
5600 5451 -5988 99.8503 0 3706 9722 9722 0 0 6013 7125 3
5600 5711 -2305 99.7016 0 3437 5788 9722 1 0 2346 13816 u
5600 5790 -3429 99.8551 0 0 3451 9722 1 2592 6044 15686 u
5600 6030 -639 90.8333 0 1944 3264 9722 1 13443 14765 18650 u
5600 6178 -635 90.7576 0 1944 3264 9722 0 715 2036 8346 u
5600 6257 -1427 96.4666 0 7939 9722 9722 1 4433 6234 6234 3
5600 6317 -1975 98.5335 0 0 2185 9722 1 0 2182 28362 5
5600 6353 -1908 99.639 0 3437 5377 9722 1 0 1939 7745 u
5600 6537 -1969 99.6486 0 1456 3451 9722 0 0 1992 12289 u
5600 6704 -634 90.3788 0 1944 3264 9722 0 6661 7981 10295 u
5600 6884 -3408 99.7972 0 0 3451 9722 0 4881 8333 8697 u
5600 7149 -3421 99.7972 0 0 3451 9722 0 3394 6846 10011 u
5600 7316 -4623 99.8706 0 5082 9722 9722 1 5376 10012 10012 3
5600 7442 -4877 99.8981 0 4816 9722 9722 1 4494 9399 9399 3
5600 7496 -6076 99.7221 0 3437 9566 9722 0 2449 8566 8566 u
5600 7636 -6272 99.9682 0 3437 9722 9722 1 2384 8669 11685 u
5600 7648 -3414 99.7682 0 0 3451 9722 0 1210 4661 15278 u
5600 7650 -3415 99.7972 0 0 3451 9722 0 4410 7862 9230 u
5600 7838 -2262 99.4015 0 3437 5790 9722 1 0 2339 12133 u
5600 7883 -6263 99.9045 0 3437 9722 9722 0 1760 8041 10570 u
5600 8266 -638 90.7576 0 1944 3264 9722 0 5317 6638 7379 u
5600 8281 -5260 99.8868 0 4416 9722 9722 1 6746 12048 12048 3
5600 8317 -1383 99.449 0 0 1454 9722 1 0 1452 13342 5
5600 8493 -638 90.7576 0 1944 3264 9722 1 8335 9656 12085 u
5600 8558 -639 90.8333 0 1944 3264 9722 0 2584 3906 7895 u
5600 8614 -638 90.7576 0 1944 3264 9722 0 7231 8552 10856 u
5600 8822 -4414 99.8205 0 5265 9722 9722 1 2937 7395 7395 3
5600 8846 -4130 99.8557 0 5564 9722 9722 1 1772 5930 5930 3
5600 9144 -1963 99.8483 0 1474 3451 9722 1 4905 6883 6883 u
5600 9637 -1132 99.9117 0 3437 4570 9722 1 0 1132 7263 u
5600 9742 -6313 99.5277 0 0 6570 9722 1 0 6564 12240 5
""")

let fp_32757* = newStringStream(

                """32757 10030 -88 70.5137 0 9470 12284 16890 0 6365 9129 11852 u
32757 10150 -88 70.5819 0 9470 12284 16890 1 5052 7819 12414 u
32757 10304 -576 89.1658 0 14841 16710 16890 0 2004 3850 4474 u
32757 10882 -822 96.9352 0 7946 9088 16890 1 5172 6314 12879 u
32757 11159 -822 96.8449 0 7946 9088 16890 0 8374 9515 12210 u
32757 11355 -586 89.2199 0 14841 16710 16890 0 4442 6288 13480 u
32757 11504 -2325 95.5389 0 10521 14130 16890 0 3051 6667 13588 u
32757 11651 -1814 77.3585 0 9034 16319 16890 0 0 6996 12426 u
32757 11871 -88 70.5712 0 9470 12284 16890 0 9781 12547 17535 u
32757 11896 -578 89.154 0 14841 16710 16890 1 8222 10066 14093 u
32757 12238 -1680 83.7716 0 10539 16319 16890 0 543 6362 12795 u
32757 12264 -2322 95.5389 0 10521 14130 16890 0 9133 12749 23097 u
32757 12350 -917 97.5986 0 7963 9133 16890 1 177 1343 18002 u
32757 12489 -88 70.5819 0 9470 12284 16890 0 5533 8300 11821 u
32757 1256 -88 70.4241 0 9470 12284 16890 0 12319 15078 16605 u
32757 12774 -822 96.9352 0 7946 9088 16890 1 1224 2366 14720 u
32757 13446 -88 70.5564 0 9470 12284 16890 1 10766 13534 15670 u
32757 13449 -2312 95.5389 0 10521 14130 16890 1 495 4111 12691 u
32757 13906 -586 89.2258 0 14841 16710 16890 1 6620 8467 13548 u
32757 14247 -822 96.9352 0 7946 9088 16890 0 5776 6918 11658 u
32757 14280 -1748 99.1699 0 0 1807 16890 0 10303 12117 12117 5
32757 15889 -586 89.2258 0 14841 16710 16890 1 2906 4753 10084 u
32757 16130 -9142 99.9021 0 7697 16890 16890 0 0 9196 12879 3
32757 16224 -1911 96.2976 0 11321 14130 16890 0 0 2818 13388 u
32757 16517 -2055 78.275 0 8659 16319 16890 1 4213 11587 11587 u
32757 16592 -822 96.8476 0 7946 9088 16890 1 7126 8269 20636 u
32757 17165 -2315 95.5112 0 10521 14130 16890 0 3039 6656 11696 u
32757 17205 -2331 95.5943 0 10521 14130 16890 0 4459 8077 13338 u
32757 17315 -6849 99.5995 0 0 6992 16890 1 0 7004 14337 5
32757 17351 -6526 99.7425 0 10287 16890 16890 1 4756 11360 11360 3
32757 17453 -945 97.774 0 7963 9133 16890 0 8534 9702 11124 u
32757 1746 -1009 94.9858 0 12375 14130 16890 0 0 1765 4378 u
32757 17791 -822 96.9352 0 7946 9088 16890 1 7123 8265 12746 u
32757 17841 -2331 95.5943 0 10521 14130 16890 1 8454 12072 12194 u
32757 18150 -1513 77.3576 0 7963 12088 16890 1 1599 5900 12029 u
32757 18257 -1233 83.0407 0 11649 16319 16890 1 10910 15630 15630 u
32757 18261 -9581 99.9169 0 7261 16890 16890 0 0 9633 10658 3
32757 18863 -822 96.9352 0 7946 9088 16890 1 5081 6223 10463 u
32757 18899 -483 88.1883 0 15431 16710 16890 0 149 1402 15504 u
32757 19248 -822 96.6725 0 7946 9088 16890 0 10049 11192 12040 u
32757 19411 -2281 78.6333 0 7963 15704 16890 1 0 7975 11259 u
32757 19865 -88 70.5819 0 9470 12284 16890 1 7307 10074 13260 u
32757 20746 -88 70.5606 0 9470 12284 16890 1 6892 9657 12231 u
32757 20836 -2557 79.5588 0 7963 16123 16890 1 138 8528 15153 u
32757 21832 -586 89.2258 0 14841 16710 16890 1 961 2808 14723 u
32757 21993 -6728 99.4662 0 0 6931 16890 0 4627 11572 11572 5
32757 22643 -558 92.4899 0 14841 16319 16890 1 326 1804 12279 u
32757 23084 -88 70.5819 0 9470 12284 16890 1 2762 5529 10233 u
32757 23126 -1115 87.6054 0 13947 16319 16890 0 1287 3700 15613 u
32757 23452 -2627 79.8348 0 7963 16319 16890 0 3511 12100 17477 u
32757 23607 -1513 77.3818 0 7963 12088 16890 0 6264 10565 10950 u
32757 24160 -3770 99.6138 0 13005 16890 16890 1 8809 12693 12693 3
32757 24511 -6174 99.697 0 0 6271 16890 1 0 6286 11824 5
32757 24604 -88 70.5351 0 9470 12284 16890 0 5358 8124 16182 u
32757 24657 -931 97.774 0 7963 9133 16890 1 904 2072 13979 u
32757 25182 -5344 99.7964 0 11486 16890 16890 0 0 5407 12202 3
32757 2556 -2330 95.5666 0 10521 14130 16890 1 1229 4846 12969 u
32757 26424 -2568 79.787 0 7963 16319 16890 1 1668 10251 16646 u
32757 26604 -88 70.5499 0 9470 12284 16890 1 5831 8595 9917 u
32757 2667 -2331 95.5943 0 10521 14130 16890 1 7160 10778 14434 u
32757 26774 -938 97.774 0 7963 9133 16890 1 921 2089 13084 u
32757 2687 -562 89.1057 0 14841 16710 16890 0 6276 8121 12733 u
32757 27131 -959 96.0916 0 14834 16319 16890 1 10918 12402 12402 u
32757 27261 -8868 99.8767 0 7967 16890 16890 0 0 8930 10983 3
32757 27360 -1068 87.552 0 13917 16319 16890 0 0 2440 14487 u
32757 27660 -586 89.2258 0 14841 16710 16890 1 1754 3601 12469 u
32757 2786 -6237 99.6843 0 0 6336 16890 0 4449 10799 10799 5
32757 27940 -2529 79.7989 0 7963 16319 16890 1 202 8787 9136 u
32757 28401 -1951 77.5079 0 7963 15250 16890 0 3255 10775 10775 u
32757 28800 -822 96.9352 0 7946 9088 16890 0 2256 3398 13840 u
32757 29149 -822 96.8449 0 7946 9088 16890 1 8190 9331 13073 u
32757 29221 -586 89.2258 0 14841 16710 16890 1 4316 6163 11604 u
32757 30080 -1112 95.2607 0 10521 12209 16890 1 0 1694 9793 u
32757 30716 -2324 95.5112 0 10521 14130 16890 0 8500 12115 32199 u
32757 30780 -88 70.49 0 9470 12284 16890 1 12104 14859 16308 u
32757 30852 -1783 77.5591 0 8911 16319 16890 1 11126 18238 18238 u
32757 30912 -2543 79.7989 0 7963 16319 16890 0 420 9006 14048 u
32757 30963 -2330 95.5666 0 10521 14130 16890 0 3421 7038 8275 u
32757 31250 -674 96.7308 0 7973 9015 16890 0 757 1797 1797 u
32757 31513 -584 89.0456 0 14841 16710 16890 1 582 2426 11054 u
32757 31726 -586 89.2258 0 14841 16710 16890 0 2855 4702 12553 u
32757 31886 -733 95.7746 0 15183 16319 16890 0 0 1136 11758 u
32757 31888 -88 70.5393 0 9470 12284 16890 1 11823 14586 16492 u
32757 32071 -940 72.3238 0 8682 12088 16890 1 1325 4389 4630 u
32757 32164 -88 70.5819 0 9470 12284 16890 0 2218 4985 14647 u
32757 32179 -1513 77.2848 0 7963 12088 16890 0 6791 11090 11935 u
32757 32758 -1513 77.3333 0 7963 12088 16890 0 8997 13297 14777 u
32757 33153 -2592 79.8109 0 7963 16319 16890 1 3138 11726 12643 u
32757 33401 -2278 78.9507 0 8333 16319 16890 0 0 8219 12518 u
32757 33428 -2303 95.7078 0 10612 14130 16890 1 7364 10890 10890 u
32757 34303 -88 70.5819 0 9470 12284 16890 0 2592 5359 12662 u
32757 34385 -452 92.9358 0 14841 15932 16890 0 11628 12718 12718 u
32757 34656 -1655 84.403 0 11068 16319 16890 0 143 5446 14035 u
32757 34681 -586 89.2258 0 14841 16710 16890 1 2319 4166 13269 u
32757 35037 -1460 77.2606 0 7963 12088 16890 1 1069 5370 11533 u
32757 35432 -583 89.1658 0 14841 16710 16890 0 9152 10998 11209 u
32757 36249 -7872 99.8108 0 8957 16890 16890 1 4932 12862 12862 3
32757 36473 -948 97.8614 0 7963 9133 16890 0 10466 11635 12672 u
32757 36720 -822 96.9352 0 7946 9088 16890 0 5060 6202 13910 u
32757 36740 -1670 76.2763 0 7963 14897 16890 1 0 7157 14967 u
32757 37118 -586 89.2258 0 14841 16710 16890 1 7261 9108 11159 u
32757 37147 -537 93.7743 0 10521 11549 16890 0 10384 11412 11412 u
32757 38074 -7823 99.8351 0 9008 16890 16890 0 0 7891 12495 3
32757 38843 -822 96.9352 0 7946 9088 16890 0 6879 8021 19915 u
32757 39007 -581 89.1481 0 14841 16710 16890 0 8498 10341 11153 u
32757 39114 -780 96.6637 0 7946 9088 16890 1 414 1553 15295 u
32757 39305 -88 70.5925 0 9470 12284 16890 1 7919 10687 12020 u
32757 39326 -88 70.42 0 9470 12284 16890 1 1879 4641 12269 u
32757 39410 -2332 79.0805 0 8293 16319 16890 1 3034 11295 11295 u
32757 40013 -458 92.5303 0 14841 15914 16890 0 8470 9541 9646 u
32757 40161 -572 89.0456 0 14841 16710 16890 0 1396 3240 12522 u
32757 40652 -1469 77.2364 0 7963 12088 16890 1 135 4430 12592 u
32757 40678 -513 88.3786 0 15248 16710 16890 1 12948 14385 14494 u
32757 40912 -2304 95.5112 0 10521 14130 16890 0 7464 11079 15033 u
32757 41734 -1906 95.501 0 10521 13455 16890 1 0 2939 14414 u
32757 4218 -88 70.5819 0 9470 12284 16890 1 2644 5411 11133 u
32757 42265 -578 89.1658 0 14841 16710 16890 1 8998 10844 11594 u
32757 42381 -948 97.8614 0 7963 9133 16890 0 11732 12901 13865 u
32757 43268 -2314 95.4835 0 10521 14130 16890 1 1685 5299 11401 u
32757 43390 -529 88.5337 0 15248 16710 16890 0 130 1569 11552 u
32757 43680 -88 70.5244 0 9470 12284 16890 0 8413 11178 13769 u
32757 4370 -2174 95.4205 0 10521 13993 16890 0 10355 13833 13833 u
32757 43778 -8716 99.5624 0 7976 16890 16890 1 9898 18811 18811 3
32757 44144 -578 89.1658 0 14841 16710 16890 1 10559 12405 14380 u
32757 4479 -2331 95.5943 0 10521 14130 16890 1 5814 9432 12200 u
32757 44992 -574 89.0396 0 14841 16710 16890 1 724 2567 11494 u
32757 44994 -88 70.5712 0 9470 12284 16890 0 715 3481 13235 u
32757 45165 -2628 79.8468 0 7963 16319 16890 0 928 9519 11728 u
32757 45547 -578 89.1057 0 14841 16710 16890 1 2622 4467 11329 u
32757 46261 -2607 79.8468 0 7963 16319 16890 1 119 8706 11062 u
32757 4646 -1580 84.2887 0 11068 16319 16890 1 11764 17064 17342 u
32757 46951 -1513 77.3333 0 7963 12088 16890 1 1727 6029 11507 u
32757 47135 -586 89.2258 0 14841 16710 16890 0 6435 8282 12518 u
32757 47458 -9728 99.8372 0 7059 16890 16890 1 462 10291 10291 3
32757 47742 -1680 96.176 0 10521 12955 16890 1 126 2558 11145 u
32757 49803 -1400 95.9453 0 10521 12642 16890 0 2132 4253 4253 u
32757 49865 -88 70.5819 0 9470 12284 16890 1 6276 9043 14507 u
32757 50087 -88 70.5606 0 9470 12284 16890 1 2797 5562 8800 u
32757 50151 -814 96.8476 0 7946 9088 16890 0 1368 2511 15088 u
32757 50356 -927 97.6864 0 7963 9133 16890 1 910 2077 15142 u
32757 50806 -586 89.1057 0 14841 16710 16890 0 1032 2877 13602 u
32757 51881 -578 89.2141 0 14841 16710 16890 0 527 2372 11781 u
32757 52024 -586 89.2741 0 14841 16710 16890 0 8837 10683 20869 u
32757 5413 -6372 99.8133 0 10464 16890 16890 1 1302 7732 7732 3
32757 5416 -948 97.8614 0 7963 9133 16890 0 11209 12378 12751 u
32757 6004 -88 70.4628 0 9470 12284 16890 1 10717 13483 13702 u
32757 6093 -12852 99.915 0 2312 15253 16890 0 0 12942 12942 contains
32757 6470 -2312 95.5389 0 10521 14130 16890 1 1226 4842 13346 u
32757 7581 -530 92.5743 0 14841 16255 16890 0 11637 13051 13051 u
32757 7813 -1135 82.3978 0 11990 16319 16890 0 0 4384 5292 u
32757 8167 -10315 99.7416 0 0 10448 16890 0 2162 12622 12622 5
32757 8493 -1586 83.6851 0 10539 16319 16890 0 521 6331 15080 u
32757 869 -1547 76.7434 0 7963 13957 16890 0 4532 10722 10961 u
32757 8815 -88 70.5712 0 9470 12284 16890 0 1535 4301 13428 u
32757 8858 -1936 77.5647 0 7963 15264 16890 0 5216 12749 12749 u
32757 9072 -88 70.5819 0 9470 12284 16890 0 6488 9255 11554 u
32757 9197 -1115 87.5632 0 13947 16319 16890 1 7446 9859 11114 u
32757 927 -7340 99.8783 0 9494 16890 16890 0 0 7401 11463 3
32757 9644 -11063 99.76 0 1276 12529 16890 0 0 11248 11248 contains
32757 9650 -1055 91.6882 0 14382 16319 16890 0 327 2287 17865 u
32757 9834 -88 70.5137 0 9470 12284 16890 1 3364 6128 9734 u
32757 9869 -1629 84.2492 0 11005 16319 16890 0 0 5375 11913 u
""")

let tp_15210* = newStringStream(

                """15210 10043 -4436 99.9102 0 8199 12654 16697 0 8427 12886 12886 u
15210 10669 -1867 99.9467 0 8199 10075 16697 1 0 1875 10491 u
15210 10963 -1494 99.7324 0 14613 16108 16697 0 3143 4638 13253 u
15210 11191 -1484 99.6653 0 14613 16108 16697 0 2972 4466 7939 u
15210 11683 -1483 99.597 0 14613 16104 16697 0 455 1944 13576 u
15210 11706 -1487 99.6653 0 14613 16108 16697 1 3090 4584 10942 u
15210 11838 -1133 99.3266 0 14613 15805 16697 0 11462 12650 12650 u
15210 12156 -8451 99.8235 0 8199 16697 16697 1 4979 13480 22212 u
15210 12328 -1032 100 0 15074 16082 16697 0 9485 10493 10759 u
15210 1234 -6904 99.8847 0 9759 16697 16697 1 4647 11591 11591 3
15210 12350 -1417 99.7195 0 14682 16108 16697 0 1878 3304 13309 u
15210 12361 -1037 100 0 15074 16087 16697 0 3342 4355 7171 u
15210 1243 -1020 99.9012 0 15074 16087 16697 0 9185 10197 10362 u
15210 12563 -1485 99.6669 0 14613 16115 16697 1 7843 9344 10264 u
15210 12632 -1472 99.6656 0 14613 16108 16697 1 10907 12402 13166 u
15210 12956 -1036 99.9012 0 15074 16087 16697 1 3166 4178 4336 u
15210 13130 -1424 99.7896 0 14682 16108 16697 1 4054 5481 12052 u
15210 13148 -1205 98.5562 0 569 1885 16697 1 0 1317 14155 u
15210 13202 -4940 99.8792 0 11730 16697 16697 1 6511 11478 11478 3
15210 13364 -1498 99.8003 0 14613 16115 16697 1 4598 6101 12223 u
15210 13448 -1453 98.6979 0 14613 16149 16697 1 0 1543 14817 u
15210 1360 -1410 99.4382 0 14682 16108 16697 0 440 1864 16615 u
15210 13742 -1037 100 0 15074 16087 16697 0 2537 3550 11341 u
15210 13807 -1487 99.7324 0 14613 16108 16697 1 8919 10415 13748 u
15210 13979 -1037 100 0 15074 16087 16697 1 11461 12474 18379 u
15210 14397 -1437 99.2642 0 14613 16115 16697 0 582 2077 13260 u
15210 14542 -1492 99.7317 0 14613 16104 16697 0 1689 3180 9839 u
15210 14763 -1424 99.7896 0 14682 16108 16697 1 5040 6467 10820 u
15210 151 -1484 99.5296 0 14613 16104 16697 0 12365 13853 17058 u
15210 15120 -1413 99.6491 0 14682 16108 16697 0 12081 13506 14069 u
15210 15411 -1467 99.6656 0 14613 16108 16697 1 3995 5490 11894 u
15210 15604 -2411 99.5903 0 0 2441 16697 1 0 2445 11655 5
15210 15607 -1412 99.6491 0 14682 16108 16697 0 255 1680 13027 u
15210 15609 -6202 99.7136 0 8199 14483 16697 1 0 6286 20520 u
15210 1621 -1494 99.7324 0 14613 16108 16697 0 2730 4225 11714 u
15210 16294 -1037 100 0 15074 16087 16697 1 6733 7746 11954 u
15210 16451 -1487 99.7324 0 14613 16108 16697 0 7388 8884 10493 u
15210 16499 -1367 99.568 0 14717 16108 16697 1 13634 15023 15023 u
15210 16854 -1495 99.7988 0 14613 16104 16697 1 451 1943 2966 u
15210 16970 -1417 99.7195 0 14682 16108 16697 1 7613 9039 11074 u
15210 1747 -3141 99.8737 0 8199 11367 16697 1 0 3168 13606 u
15210 17514 -3444 99.8847 0 8199 11667 16697 0 9907 13379 13379 u
15210 17721 -1485 99.6669 0 14613 16115 16697 0 2104 3605 13453 u
15210 17884 -1424 99.7896 0 14682 16108 16697 0 8037 9464 9988 u
15210 17943 -1495 99.7988 0 14613 16104 16697 1 7847 9339 10158 u
15210 17952 -1475 99.5984 0 14613 16108 16697 0 6549 8043 9346 u
15210 17965 -1423 99.7195 0 14682 16108 16697 0 3921 5347 11195 u
15210 18041 -1021 99.8012 0 15074 16082 16697 1 2079 3085 15370 u
15210 18068 -1492 99.4674 0 14613 16115 16697 1 4832 6339 14458 u
15210 18176 -1484 98.5488 0 14592 16108 16697 0 0 1528 11038 u
15210 18416 -8469 99.9059 0 8199 16697 16697 1 4272 12776 14186 u
15210 18545 -1486 99.5981 0 14613 16108 16697 0 1900 3393 10861 u
15210 18612 -1505 100 0 8199 9708 16697 1 0 1509 12394 u
15210 18701 -1487 99.7324 0 14613 16108 16697 1 9108 10604 16370 u
15210 18764 -1490 99.6669 0 14613 16115 16697 0 2104 3605 12141 u
15210 18851 -8467 99.8941 0 8199 16697 16697 0 2356 10859 11604 u
15210 19023 -1487 99.7324 0 14613 16108 16697 0 10274 11770 13998 u
15210 19080 -7381 99.8654 0 9266 16697 16697 0 0 7433 12626 3
15210 19213 -1424 99.7896 0 14682 16108 16697 0 7036 8463 12823 u
15210 19584 -8465 99.8941 0 8199 16697 16697 0 1157 9658 10905 u
15210 196 -5166 99.9036 0 11508 16697 16697 0 0 5192 12462 3
15210 19711 -1472 99.6669 0 14613 16115 16697 1 307 1808 12665 u
15210 1980 -1438 98.9218 0 14613 16104 16697 1 3994 5478 12502 u
15210 20008 -8441 99.8235 0 8199 16697 16697 0 112 8609 13951 u
15210 2066 -1415 99.5792 0 14682 16108 16697 0 2717 4143 15176 u
15210 20752 -1424 99.7896 0 14682 16108 16697 0 4726 6153 12533 u
15210 20775 -1424 99.7896 0 14682 16108 16697 1 9524 10951 14828 u
15210 20792 -5496 99.7838 0 8199 13753 16697 1 0 5550 13110 u
15210 20872 -1424 99.7896 0 14682 16108 16697 1 3119 4546 11152 u
15210 21052 -1464 97.8764 0 14613 16170 16697 1 0 1554 8858 u
15210 21140 -1037 100 0 15074 16087 16697 0 5957 6970 10728 u
15210 21207 -5974 99.7181 0 10667 16697 16697 0 0 6031 10068 3
15210 21248 -1472 99.6656 0 14613 16108 16697 0 703 2198 10215 u
15210 21257 -1037 100 0 15074 16087 16697 1 9866 10879 12606 u
15210 21258 -1487 99.7324 0 14613 16108 16697 0 10957 12453 13374 u
15210 21269 -1493 99.6647 0 14613 16104 16697 1 8568 10062 13395 u
15210 21377 -1417 99.7195 0 14682 16108 16697 1 10872 12298 15024 u
15210 21738 -3274 99.727 0 0 3297 16697 0 8655 11957 11957 5
15210 21880 -1424 99.7896 0 14682 16108 16697 1 15137 16564 20863 u
15210 21957 -1881 99.7898 0 14794 16697 16697 1 8247 10150 10150 3
15210 21991 -5695 99.5726 0 8199 14051 16697 1 0 5849 13429 u
15210 22802 -1037 100 0 15074 16087 16697 0 480 1493 10712 u
15210 22851 -1037 100 0 15074 16087 16697 0 5416 6429 16709 u
15210 22992 -1037 100 0 15074 16087 16697 0 1980 2993 7446 u
15210 23003 -8464 99.8941 0 8199 16697 16697 1 4986 13489 17225 u
15210 23122 -1424 99.7896 0 14682 16108 16697 0 10275 11702 12801 u
15210 23158 -1028 99.8043 0 15080 16104 16697 1 10331 11353 11353 u
15210 23283 -7506 99.8676 0 9142 16697 16697 1 4865 12424 12424 3
15210 23734 -1015 99.6094 0 15076 16104 16697 1 10432 11456 11456 u
15210 23790 -1472 99.5296 0 14613 16104 16697 0 200 1688 10045 u
15210 2381 -1476 99.6656 0 14613 16108 16697 0 7600 9095 9853 u
15210 23981 -7734 99.8841 0 8199 15964 16697 1 0 7768 9992 u
15210 23993 -1483 99.7305 0 14624 16108 16697 1 11251 12735 17907 u
15210 24076 -1472 99.6656 0 14613 16108 16697 1 9582 11077 11319 u
15210 2429 -1487 99.7324 0 14613 16108 16697 0 2015 3511 11705 u
15210 24475 -4798 99.8755 0 11879 16697 16697 0 0 4820 13281 3
15210 24646 -1027 100 0 15084 16087 16697 0 2902 3905 9681 u
15210 24692 -892 98.1007 0 14629 15682 16697 1 0 1057 2839 u
15210 24884 -3543 99.6651 0 0 3583 16697 1 0 3585 10203 5
15210 24955 -1498 99.8003 0 14613 16115 16697 1 5686 7189 13519 u
15210 25034 -8466 99.8941 0 8199 16697 16697 0 492 8995 11651 u
15210 25047 -1480 99.7317 0 14613 16104 16697 1 3198 4689 11012 u
15210 2549 -1498 99.8003 0 14613 16115 16697 0 4692 6195 10726 u
15210 25882 -1498 99.8003 0 14613 16115 16697 0 2406 3909 10467 u
15210 25963 -6393 99.6606 0 10212 16697 16697 1 7921 14403 14403 3
15210 25994 -1076 99.1877 0 14989 16104 16697 1 16364 17472 17472 u
15210 26059 -1411 99.6491 0 14682 16108 16697 0 1680 3105 10670 u
15210 264 -1495 99.7337 0 14613 16115 16697 0 2755 4257 9999 u
15210 26511 -7866 99.8102 0 8199 16105 16697 1 0 7903 9477 u
15210 26545 -7518 99.9072 0 8199 15742 16697 1 0 7548 10597 u
15210 27348 -8447 99.8706 0 8199 16697 16697 1 4208 12709 16756 u
15210 27389 -1423 99.7195 0 14682 16108 16697 0 348 1774 10421 u
15210 27414 -1441 99.3988 0 14613 16115 16697 0 9449 10946 13126 u
15210 27747 -1022 99.8012 0 15074 16082 16697 1 2075 3081 13969 u
15210 27762 -5289 99.7938 0 8199 13536 16697 1 0 5334 14084 u
15210 27901 -1037 100 0 15074 16087 16697 0 5409 6422 13926 u
15210 2793 -1424 99.7896 0 14682 16108 16697 1 895 2322 13259 u
15210 281 -1467 99.6656 0 14613 16108 16697 1 4898 6393 12598 u
15210 28154 -1219 98.6186 0 14804 16115 16697 1 24698 26001 26001 u
15210 28839 -3112 99.5025 0 12664 15888 16697 1 0 3216 3216 contains
15210 2901 -1036 99.9012 0 15074 16087 16697 0 1318 2330 10491 u
15210 2910 -1047 99.6327 0 14613 15702 16697 1 0 1089 11146 u
15210 29439 -1037 100 0 15074 16087 16697 1 7561 8574 12652 u
15210 29975 -1100 99.3092 0 14613 15773 16697 0 12939 14097 14097 u
15210 30010 -1486 99.6656 0 14613 16108 16697 0 9985 11480 11867 u
15210 30087 -1489 99.7317 0 14613 16104 16697 0 5812 7303 11669 u
15210 30241 -3654 99.8374 0 8199 11889 16697 1 0 3692 10639 u
15210 30298 -1487 99.7324 0 14613 16108 16697 0 6118 7614 12265 u
15210 30316 -8427 99.7647 0 8199 16697 16697 1 3623 12121 15196 u
15210 30381 -2938 99.966 0 8199 11143 16697 0 10315 13260 13260 u
15210 30630 -1492 99.7317 0 14613 16104 16697 0 8190 9681 10800 u
15210 30763 -1424 99.7896 0 14682 16108 16697 1 8933 10360 13125 u
15210 31119 -5082 99.9021 0 11590 16697 16697 0 0 5108 14417 3
15210 31126 -1488 99.6653 0 14613 16108 16697 0 4282 5776 8697 u
15210 31417 -1494 99.7324 0 14613 16108 16697 0 5480 6975 12849 u
15210 31463 -1270 98.6419 0 14698 16104 16697 1 4586 5985 5985 u
15210 31672 -8469 99.9059 0 8199 16697 16697 0 1738 10242 14725 u
15210 31778 -1487 99.6653 0 14613 16108 16697 0 5986 7480 15834 u
15210 31842 -1490 99.6653 0 14613 16108 16697 1 10501 11995 13171 u
15210 31934 -7940 99.8871 0 8199 16171 16697 1 0 7979 10993 u
15210 31942 -5629 99.8236 0 11029 16697 16697 1 8275 13945 13945 3
15210 31960 -3106 99.8406 0 8199 11337 16697 1 0 3137 19285 u
15210 32039 -1493 99.6647 0 14613 16104 16697 1 8687 10179 11955 u
15210 32068 -3081 99.8396 0 13580 16697 16697 0 0 3120 8325 3
15210 32121 -1491 99.6644 0 14613 16104 16697 1 9526 11016 19884 u
15210 32211 -1482 99.5302 0 14613 16104 16697 0 3489 4979 18446 u
15210 32368 -1487 99.6653 0 14613 16108 16697 1 4783 6277 12280 u
15210 32374 -1424 99.7896 0 14682 16108 16697 1 5134 6561 12499 u
15210 32631 -1037 100 0 15074 16087 16697 1 10339 11352 15156 u
15210 32669 -1488 99.6653 0 14613 16108 16697 1 10680 12174 16474 u
15210 32891 -1344 99.6361 0 14730 16104 16697 1 10245 11621 11621 u
15210 32901 -1482 99.6 0 14613 16115 16697 1 8742 10242 12690 u
15210 33024 -1400 98.8598 0 14613 16104 16697 1 14180 15674 17940 u
15210 33185 -4588 99.8915 0 12088 16697 16697 1 6427 11037 11037 3
15210 3320 -1078 98.9492 0 14613 15761 16697 1 0 1142 13905 u
15210 33328 -1026 99.8022 0 15074 16087 16697 1 5115 6126 9935 u
15210 33424 -1373 99.639 0 0 1385 16697 1 0 1388 10882 5
15210 33931 -1475 99.4008 0 14613 16115 16697 0 2846 4351 11526 u
15210 34125 -1494 99.7317 0 14613 16104 16697 1 3832 5323 11678 u
15210 34249 -1064 99.726 0 14613 15708 16697 1 0 1096 11519 u
15210 3432 -1494 99.7324 0 14613 16108 16697 0 3014 4509 13300 u
15210 34412 -1416 99.6491 0 14682 16108 16697 1 6227 7652 9581 u
15210 34561 -1967 100 0 8199 10167 16697 0 8158 10126 10126 u
15210 34569 -1037 100 0 15074 16087 16697 0 5867 6880 11159 u
15210 34613 -1498 99.8003 0 14613 16115 16697 1 1989 3492 11495 u
15210 34638 -1480 99.6653 0 14613 16108 16697 0 4529 6023 13960 u
15210 34722 -1480 99.5318 0 14613 16108 16697 1 415 1912 11750 u
15210 34784 -1037 100 0 15074 16087 16697 0 3290 4303 11649 u
15210 35043 -7169 99.875 0 8199 15401 16697 0 3833 11038 11038 u
15210 35160 -1398 99.7193 0 14682 16108 16697 1 928 2353 17110 u
15210 35285 -1037 100 0 15074 16087 16697 1 8540 9553 10740 u
15210 35374 -1498 99.8003 0 14613 16115 16697 0 4756 6259 8246 u
15210 35418 -1040 99.8058 0 15077 16108 16697 0 0 1030 12188 u
15210 35660 -1478 99.5962 0 14620 16108 16697 1 14213 15699 17900 u
15210 35866 -1424 99.7896 0 14682 16108 16697 0 12890 14317 19842 u
15210 3594 -1412 99.6494 0 14682 16108 16697 0 11608 13035 14569 u
15210 36047 -1033 99.9012 0 15074 16087 16697 0 6440 7452 8688 u
15210 3605 -6642 99.1536 0 8199 15200 16697 0 6603 13574 13574 u
15210 36139 -1451 99.4634 0 14613 16104 16697 1 9309 10801 13069 u
15210 36263 -1400 99.5787 0 14682 16108 16697 0 13578 15002 16538 u
15210 36380 -1485 99.5976 0 14613 16104 16697 0 10018 11509 12221 u
15210 36407 -1487 99.7324 0 14613 16108 16697 0 1901 3397 10747 u
15210 36574 -1486 99.5976 0 14613 16104 16697 1 1375 2866 10428 u
15210 36688 -1495 99.7988 0 14613 16104 16697 0 2219 3711 8562 u
15210 36885 -6943 99.8571 0 8199 15195 16697 0 3780 10778 10778 u
15210 36989 -1379 99.4378 0 14682 16108 16697 1 544 1967 11123 u
15210 36992 -1420 99.7195 0 14682 16108 16697 0 10620 12046 16485 u
15210 37042 -1424 99.7896 0 14682 16108 16697 1 1054 2481 14093 u
15210 3706 -1496 99.7337 0 14613 16115 16697 1 2578 4080 10875 u
15210 37068 -1490 99.6669 0 14613 16115 16697 0 6836 8337 11951 u
15210 37 -1037 100 0 15074 16087 16697 0 3536 4549 9535 u
15210 37241 -1487 99.6653 0 14613 16108 16697 0 5520 7014 13952 u
15210 37542 -3412 99.767 0 0 3433 16697 1 0 3439 7608 5
15210 3758 -1027 97.6537 0 14637 15873 16697 0 754 1996 1996 u
15210 37759 -8446 99.8353 0 8199 16697 16697 1 1961 10463 11530 u
15210 37891 -1579 100 0 8199 9778 16697 1 0 1579 12086 u
15210 37942 -1037 100 0 15074 16087 16697 0 4170 5183 6873 u
15210 38098 -1495 99.7988 0 14613 16104 16697 0 8942 10434 13272 u
15210 38128 -1467 99.533 0 14613 16115 16697 0 9022 10521 10836 u
15210 38175 -1494 99.7324 0 14613 16108 16697 0 4591 6086 13945 u
15210 3871 -3039 99.8694 0 13634 16697 16697 1 8839 11904 11904 3
15210 3873 -7301 99.8773 0 8199 15536 16697 1 0 7340 12935 u
15210 4019 -1493 99.7317 0 14613 16104 16697 1 8486 9979 13050 u
15210 4061 -1037 100 0 15074 16087 16697 1 4529 5542 10686 u
15210 4144 -2052 99.6146 0 0 2076 16697 1 0 2082 5422 5
15210 4254 -1495 99.7988 0 14613 16104 16697 0 8447 9939 13088 u
15210 4283 -5374 99.8704 0 11294 16697 16697 0 0 5408 14144 3
15210 4334 -1414 99.6491 0 14682 16108 16697 0 440 1865 11317 u
15210 4375 -3661 99.7559 0 0 3687 16697 0 7618 11310 11310 5
15210 4578 -1037 100 0 15074 16087 16697 1 8749 9762 11604 u
15210 5089 -1494 99.7324 0 14613 16108 16697 0 4478 5973 11907 u
15210 5109 -6513 99.8779 0 8199 14749 16697 0 5462 12014 12014 u
15210 5445 -1494 99.7324 0 14613 16108 16697 0 3765 5260 15518 u
15210 5544 -7936 99.8132 0 8665 16697 16697 0 0 8032 15229 3
15210 5945 -1492 99.5981 0 14613 16108 16697 1 438 1931 10393 u
15210 6145 -1424 99.7896 0 14682 16108 16697 0 5757 7184 12322 u
15210 619 -1413 99.2632 0 14613 16108 16697 1 1616 3109 8008 u
15210 6231 -4523 99.8461 0 12148 16697 16697 0 0 4548 14346 3
15210 6294 -1487 99.7324 0 14613 16108 16697 0 7318 8814 11632 u
15210 6316 -1487 99.7324 0 14613 16108 16697 0 7324 8820 12702 u
15210 6448 -5513 99.9096 0 8199 13731 16697 0 5982 11517 11517 u
15210 6850 -5993 99.7852 0 0 6051 16697 0 5652 11708 11708 5
15210 7156 -1460 99.5981 0 14613 16108 16697 1 5601 7094 10135 u
15210 7401 -8191 99.7352 0 8199 16509 16697 0 8282 16589 16589 u
15210 7678 -1036 99.9012 0 15074 16087 16697 0 11315 12327 13562 u
15210 7764 -4994 99.8805 0 11678 16697 16697 0 0 5019 7819 3
15210 8025 -8432 99.8588 0 8199 16697 16697 1 1082 9580 11878 u
15210 8278 -1030 99.9012 0 15074 16087 16697 1 9369 10381 11690 u
15210 8679 -1424 99.7896 0 14682 16108 16697 1 3516 4943 13801 u
15210 8843 -7966 99.8501 0 8199 16205 16697 1 0 8006 12656 u
15210 8845 -1037 100 0 15074 16087 16697 0 585 1598 12971 u
15210 8864 -1487 99.7324 0 14613 16108 16697 0 5164 6660 10234 u
15210 8975 -1424 99.7896 0 14682 16108 16697 0 148 1575 11345 u
15210 9108 -1494 99.7317 0 14613 16104 16697 0 6552 8043 11260 u
15210 9175 -1493 99.6653 0 14613 16108 16697 1 8652 10146 11972 u
15210 9735 -1025 99.7038 0 15074 16087 16697 0 8695 9711 11843 u
""")

let fp_17238* = newStringStream(

                """17238 10178 -1466 98.0911 0 11260 12884 14755 1 1391 3017 5766 u
17238 1019 -1436 98.0087 0 11260 12867 14755 0 5969 7577 11604 u
17238 10483 -1465 98.0296 0 11260 12884 14755 0 8576 10201 11499 u
17238 11119 -1459 98.0296 0 11260 12884 14755 1 6309 7934 10657 u
17238 11264 -1453 98.0296 0 11260 12884 14755 1 569 2194 10535 u
17238 11422 -5226 98.8526 0 9177 14755 14755 0 0 5591 12825 3
17238 11434 -1464 97.968 0 11260 12884 14755 0 3183 4807 20025 u
17238 11752 -1441 97.7791 0 11260 12884 14755 0 648 2269 19009 u
17238 12441 -1370 97.0321 0 11281 12935 14755 0 12892 14543 14543 u
17238 12783 -1466 98.0911 0 11260 12884 14755 1 8164 9790 13049 u
17238 12901 -1010 97.8336 0 11260 12414 14755 0 10583 11737 11737 u
17238 13073 -5097 98.8106 0 9290 14755 14755 1 6680 12159 12159 3
17238 1368 -1713 98.8839 0 0 1792 14755 0 9836 11628 11628 5
17238 14 -1466 98.0911 0 11260 12884 14755 0 2212 3838 11266 u
17238 14508 -1464 98.0911 0 11260 12884 14755 1 841 2465 12726 u
17238 14615 -1466 98.0911 0 11260 12884 14755 1 4941 6567 10385 u
17238 1469 -2510 99.0476 0 0 2628 14755 1 0 2625 10547 5
17238 15582 -1943 98.506 0 0 2084 14755 1 0 2075 17768 5
17238 16387 -1175 97.0447 0 11260 12721 14755 1 0 1455 5714 u
17238 16837 -1466 98.0911 0 11260 12884 14755 1 5468 7094 11065 u
17238 17525 -1449 98.0709 0 11260 12867 14755 0 3113 4722 12470 u
17238 18180 -1422 98.0087 0 11260 12867 14755 0 1036 2644 10326 u
17238 18660 -1449 97.968 0 11260 12884 14755 1 1952 3576 12120 u
17238 19038 -1466 98.0911 0 11260 12884 14755 1 8320 9946 13626 u
17238 19418 -7457 98.9011 0 6838 14755 14755 1 2712 10639 10639 3
17238 20314 -1456 97.9038 0 11260 12884 14755 0 652 2274 15805 u
17238 20712 -1449 98.0709 0 11260 12867 14755 0 8945 10554 12998 u
17238 21217 -1411 97.9439 0 11260 12867 14755 0 1587 3192 15449 u
17238 21442 -4313 98.0928 0 10036 14755 14755 0 0 4754 12982 3
17238 21848 -1453 98.0296 0 11260 12884 14755 1 1465 3090 11008 u
17238 23047 -6607 98.9735 0 7741 14755 14755 1 3332 10360 10360 3
17238 2480 -1465 98.0296 0 11260 12884 14755 0 14412 16037 18700 u
17238 25498 -1442 98.0087 0 11260 12867 14755 0 1707 3314 13341 u
17238 26570 -1459 98.0296 0 11260 12884 14755 0 3129 4754 10662 u
17238 26680 -1459 98.0296 0 11260 12884 14755 1 8522 10147 10295 u
17238 27330 -1422 98.0087 0 11260 12867 14755 0 736 2344 10447 u
17238 27761 -1463 97.968 0 11260 12884 14755 0 3761 5387 8733 u
17238 29049 -1456 98.0296 0 11260 12884 14755 1 2912 4539 13276 u
17238 29747 -8845 99.1416 0 5435 14755 14755 0 0 9325 10379 3
17238 30138 -1459 98.0296 0 11260 12884 14755 0 2897 4522 10987 u
17238 30149 -1466 98.0911 0 11260 12884 14755 0 6646 8272 13635 u
17238 30193 -1466 98.0911 0 11260 12884 14755 1 8111 9737 12105 u
17238 30316 -1464 98.0296 0 11260 12884 14755 1 16564 18189 20741 u
17238 30604 -1417 97.818 0 11260 12867 14755 0 1017 2621 12082 u
17238 32274 -1435 97.9439 0 11260 12867 14755 0 7049 8654 10299 u
17238 33443 -1461 97.968 0 11260 12884 14755 0 3264 4888 17138 u
17238 3347 -1459 98.0296 0 11260 12884 14755 1 5798 7423 10998 u
17238 33641 -5985 99.3252 0 0 6224 14755 1 0 6233 10099 5
17238 34495 -1466 98.0911 0 11260 12884 14755 1 7439 9065 12561 u
17238 34586 -2525 98.6296 0 0 2711 14755 0 12520 15220 15220 5
17238 3494 -1385 97.622 0 11281 12884 14755 0 6175 7773 7903 u
17238 35082 -4496 98.5519 0 9852 14755 14755 0 0 4916 5893 3
17238 35379 -1176 97.7695 0 11260 12606 14755 0 10116 11461 11461 u
17238 35899 -1338 97.1375 0 11260 12867 14755 0 12571 14180 15464 u
17238 35935 -1464 98.0296 0 11260 12884 14755 0 4591 6216 11684 u
17238 36203 -12912 99.2718 0 0 13457 14755 0 1427 14895 14895 5
17238 36903 -1252 97.9695 0 11502 12884 14755 0 0 1379 4365 u
17238 36960 -1422 98.0087 0 11260 12867 14755 0 623 2231 11078 u
17238 37535 -1459 98.0296 0 11260 12884 14755 0 3267 4892 10495 u
17238 37925 -1459 98.0296 0 11260 12884 14755 0 2367 3992 11138 u
17238 38451 -1106 97.7474 0 11260 12503 14755 1 0 1243 19959 u
17238 38793 -1455 98.0296 0 11260 12884 14755 1 1333 2958 12399 u
17238 39103 -1176 97.9104 0 11260 12603 14755 1 0 1340 10756 u
17238 39115 -1446 97.9051 0 11260 12884 14755 0 8111 9734 12778 u
17238 39362 -1466 98.0911 0 11260 12884 14755 0 3625 5251 22036 u
17238 39732 -1466 98.0911 0 11260 12884 14755 1 4789 6415 10736 u
17238 39875 -9102 99.201 0 5243 14755 14755 0 0 9519 13278 3
17238 4683 -10281 99.4422 0 0 10578 14755 0 392 10975 10975 5
17238 5080 -1484 98.1527 0 11260 12884 14755 1 10345 11970 12718 u
17238 5250 -1459 98.0296 0 11260 12884 14755 1 9557 11182 11391 u
17238 6327 -1439 98.0296 0 11260 12884 14755 0 3746 5371 12949 u
17238 6772 -1448 97.9465 0 11260 12867 14755 0 10737 12348 12555 u
17238 7313 -1074 99.0196 0 0 1122 14755 1 0 1123 15031 5
17238 7771 -5524 98.8796 0 8864 14755 14755 0 0 5902 13405 3
17238 9749 -1464 98.0296 0 11260 12884 14755 0 6824 8449 10559 u
17238 9975 -1466 98.0911 0 11260 12884 14755 1 5548 7174 7699 u
""")
