\ Multiple optimal solutions: the objective is parallel to the clay constraint
\ Expected: every point of the segment (30,0)-(24,8) is optimal, Z = 1200
Maximize
 Z: 40 x1 + 30 x2
Subject To
 labour: x1 + 2 x2 <= 40
 clay:   4 x1 + 3 x2 <= 120
End
