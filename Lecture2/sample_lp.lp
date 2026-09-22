\ ===========================================================================
\  sample_lp.lp  --  the product-mix problem of the lecture
\
\  An electronics manufacturer makes DVD players (x1) and mobile phones (x2).
\  Each DVD player needs 4 hours of electronic time and 2 hours of assembly
\  and earns $7 profit; each mobile phone needs 3 and 1 hours and earns $5.
\  240 hours of electronic time and 100 hours of assembly are available.
\
\  Expected answer:  x1 = 30, x2 = 40, Z = 410
\
\  FORMAT NOTES
\    * lines beginning with a backslash are comments
\    * sections: Maximize | Minimize , Subject To , Bounds , Free , End
\    * operators: <= >= =   (also =< => < >)
\    * every variable is >= 0 unless you say otherwise in Bounds / Free
\    * name a constraint by prefixing it with  name:
\ ===========================================================================

Maximize
 profit: 7 x1 + 5 x2

Subject To
 electronics: 4 x1 + 3 x2 <= 240
 assembly:    2 x1 + 1 x2 <= 100

Bounds
 x1 >= 0
 x2 >= 0

End
