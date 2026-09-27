\ Hot-dog mixture -- blending, with an EQUALITY and a RATIO constraint
\ x1 = kg of chicken, x2 = kg of beef, 1 ton batch
\ the ratio  x1/x2 >= 2/1  is written linearly as  x1 - 2 x2 >= 0
\ Expected: x1 = 800, x2 = 200, Z = 3400
Minimize
 Z: 3 x1 + 5 x2
Subject To
 batch_size:      x1 + x2 = 1000
 chicken_minimum: x1 >= 500
 beef_minimum:    x2 >= 200
 ratio_2_to_1:    x1 - 2 x2 >= 0
End
