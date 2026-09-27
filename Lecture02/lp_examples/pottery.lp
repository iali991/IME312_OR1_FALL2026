\ Beaver Creek Pottery Company -- maximisation
\ x1 = bowls per day, x2 = mugs per day
\ Expected: x1 = 24, x2 = 8, Z = 1360
Maximize
 Z: 40 x1 + 50 x2
Subject To
 labour: 1 x1 + 2 x2 <= 40
 clay:   4 x1 + 3 x2 <= 120
End
