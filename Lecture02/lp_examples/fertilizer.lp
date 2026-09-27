\ Fertilizer blending -- minimisation
\ x1 = bags of Super-Gro, x2 = bags of Crop-Quick
\ Expected: x1 = 0, x2 = 8, Z = 24
Minimize
 Z: 6 x1 + 3 x2
Subject To
 nitrogen:  2 x1 + 4 x2 >= 16
 phosphate: 4 x1 + 3 x2 >= 24
End
