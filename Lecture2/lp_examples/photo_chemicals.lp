\ Photo-developing chemicals -- minimisation with >= constraints
\ x1 = tons of black-and-white chemical, x2 = tons of colour chemical
\ Expected: x1 = 40, x2 = 20, Z = 160000
Minimize
 cost: 2500 x1 + 3000 x2
Subject To
 bw_requirement:     x1 >= 30
 colour_requirement: x2 >= 20
 total_tonnage:      x1 + x2 >= 60
End
