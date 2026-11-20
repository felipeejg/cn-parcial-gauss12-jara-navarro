una flop por suma, resta, multiplicación o división del núcleo Las copias, asignaciones y mediciones se excluyen

## sin_pivoteo 

| paso | mult/div | sumas/restas | subtotal |
|---|---:|---:|---:|
| 1 | 143 | 132 | 275 |
| 2 | 120 | 110 | 230 |
| 3 | 99 | 90 | 189 |
| 4 | 80 | 72 | 152 |
| 5 | 63 | 56 | 119 |
| 6 | 48 | 42 | 90 |
| 7 | 35 | 30 | 65 |
| 8 | 24 | 20 | 44 |
| 9 | 15 | 12 | 27 |
| 10 | 8 | 6 | 14 |
| 11 | 3 | 2 | 5 |
| regresiva | 78 | 66 | 144 |
| total | 716 | 638 | 1354 |

total aritmético principal = 1354 flops (-8.88% respecto a la cota 1486; ahorro de 132 flops por no operar los 66 ceros subdiagonales)

sobrecosto separado: CostoPivote(comparaciones_escalas=0, comparaciones_candidatos=0, cocientes=0, intercambios=0).

## parcial_escalado 

| paso | mult/div | sumas/restas | subtotal |
|---|---:|---:|---:|
| 1 | 143 | 132 | 275 |
| 2 | 120 | 110 | 230 |
| 3 | 99 | 90 | 189 |
| 4 | 80 | 72 | 152 |
| 5 | 63 | 56 | 119 |
| 6 | 48 | 42 | 90 |
| 7 | 35 | 30 | 65 |
| 8 | 24 | 20 | 44 |
| 9 | 15 | 12 | 27 |
| 10 | 8 | 6 | 14 |
| 11 | 3 | 2 | 5 |
| regresiva | 78 | 66 | 144 |
| total | 716 | 638 | 1354 |

total aritmético principal = 1354 flops (-8.88% respecto a la cota 1486; ahorro de 132 flops por no operar los 66 ceros subdiagonales)

sobrecosto separado: CostoPivote(comparaciones_escalas=132, comparaciones_candidatos=66, cocientes=77, intercambios=8).

