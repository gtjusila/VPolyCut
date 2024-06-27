\ Model Seperating LP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C2 + C3
Subject To
 point1: - 0.50001 C1 = 1
 point2: - 125003.0000091811 C0 - 0.50001 C1 = 1
 R2: - C0 - C2 <= 0
 R3: C0 - C2 <= 0
 R4: - C1 - C3 <= 0
 R5: C1 - C3 <= 0
Bounds
 -10000 <= C0 <= 10000
 -10000 <= C1 <= 10000
 -10000 <= C2 <= 10000
 -10000 <= C3 <= 10000
End
