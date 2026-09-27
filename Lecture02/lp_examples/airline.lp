\ Airline flight planning -- maximisation
\ x1 = daily flights to Rome, x2 = daily flights to Dublin
\ Expected: x1 = 6, x2 = 6, Z = 27000
Maximize
 Z: 2500 x1 + 2000 x2
Subject To
 gate:   x1 + x2 <= 12
 labour: 15 x1 + 10 x2 <= 150
 market: x1 <= 9
End
