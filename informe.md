# Prueba Parcial: El Pivote — Eliminación Gaussiana con Mínimo Costo Operacional y Máxima Exactitud de Mantisa

## 1. Identificación

| Campo | Detalle |
|:---|:---|
| **Estudiantes** | Felipe Jara Guzmán · Kevin Navarro |
| **Asignatura** | Computación Numérica (INF-321) — Sección S2 |
| **Carrera** | Ingeniería Civil Informática · Universidad Católica del Maule |
| **Docentes** | Juan Manuel Cabezas Contardo · Pedro Pesante |
| **Fecha** | Talca, noviembre de 2026 |
| **Repositorio Git** | `https://github.com/felipeejg/cn-parcial-gauss12-jara-navarro` |

---

## 2. Planteamiento

Se resuelve un sistema lineal denso $Ax = b$ de dimensión $n = 12$ ($144$ de $144$ entradas no nulas en $A$, entregadas con 9 cifras significativas exactas, y vector $b$ entregado con 20 dígitos significativos). Al no existir ceros estructurales en $A$, el ahorro operacional no puede venir de saltarlos, sino exclusivamente de evitar trabajo redundante y reutilizar una única factorización.

El sistema presenta tres dificultades numéricas que ponen a prueba la aritmética IEEE 754 de doble precisión ($\varepsilon_{\text{mach}} \approx 2.2204 \times 10^{-16}$, 53 bits de mantisa):

1. **Escalas de fila incompatibles:** las magnitudes por fila conviven entre $\sim 10^{-5}$ (filas 1, 5 y 9) y $\sim 10^{5}$ (filas 2, 6, 8 y 12). Al restar un múltiplo de una fila grande a una fila pequeña, el alineamiento de exponentes desplaza los bits bajos de la mantisa fuera de los 53 bits disponibles y los destruye.
2. **Pivote inicial minúsculo:** $a_{11} = 1.3407 \times 10^{-8}$ está doce órdenes de magnitud por debajo de $a_{21} = 1.1923 \times 10^{4}$ y trece por debajo de $a_{61} = 3.6904 \times 10^{5}$. Ningún pivote se anula sin pivoteo, por lo que el sistema se resuelve y permite comparar cuantitativamente ambas corridas.
3. **Casi dependencia lineal oculta:** la fila 7 se encuentra a una perturbación de orden $10^{-4}$ de ser combinación lineal de las filas 3 y 5, elevando el número de condición a $\kappa_2(A) \approx 3.0847 \times 10^{13}$ aunque $\det(A) \approx 1.3646 \times 10^{9} \neq 0$.

---

## 3. Estrategia de eliminación

Se implementó eliminación gaussiana densa en Python puro (`src/gauss.py`) con **pivoteo parcial escalado** y almacenamiento de una única descomposición $PA = LU$:

* **Escalamiento previo por fila:** antes de la eliminación se calcula una sola vez $s_i = \max_{1 \le j \le n} |a_{ij}|$ (`escalas_filas`). En cada paso $k$, se elige la fila $p \ge k$ que maximiza $|a_{ik}^{(k)}|/s_i$ (`elegir_pivote`), evitando que las filas de escala $10^5$ dominen artificialmente la selección de pivote.
* **Un multiplicador por fila:** en cada paso $k = 1, \dots, 11$ y fila $i = k+1, \dots, 12$, el multiplicador $m = a_{ik}^{(k)} / a_{kk}^{(k)}$ se calcula una sola vez, se guarda en $L_{ik}$ y se reutiliza en toda la fila $i$ y en el término independiente $c_i$.
* **Barrido estricto:** las entradas subdiagonales anuladas se fijan directamente en `0.0` y el bucle interno recorre únicamente las columnas $j = k+1, \dots, 12$.
* **Única factorización sin inversa:** no se calcula $A^{-1}$. La función `correccion()` reutiliza los factores $L$, $U$ y el vector de permutación `orden` ($P$) ya calculados para ejecutar el refinamiento iterativo en $\mathcal{O}(n^2)$ operaciones.

---

## 4. Conteo de operaciones

**Convención declarada:** 1 flop equivale a 1 operación aritmética elemental (`+`, `-`, `*`, `/`) ejecutada dentro del núcleo de eliminación y sustitución regresiva (`Cuenta.dividir` y `Cuenta.actualizar` en `src/conteo.py`). Las comparaciones, cocientes de escala e intercambios de fila se miden aparte en `CostoPivote`.

### Desglose paso a paso (idéntico en flops aritméticos para ambas corridas)

| Paso | Mult/div | Sumas/restas | Subtotal Flops |
|:---|---:|---:|---:|
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
| Regresiva ($Ux = c$) | 78 | 66 | 144 |
| **Total aritmético** | **716** | **638** | **1354** |

* **Comparación con la cota:** el término dominante teórico $2n^3/3$ para $n=12$ es $1152$ flops (solo eliminación de $A$). Frente a la cota completa del enunciado ($\approx 1486$ flops incluyendo $b$ y sustitución regresiva), el algoritmo ejecuta **1354 flops** ($-8.88\%$, dentro de $\pm 10\%$). La diferencia exacta de $132$ flops ($1354 + 132 = 1486$) proviene de no ejecutar la multiplicación y resta redundantes sobre los $66$ elementos subdiagonales que se fijan en `0.0`.
* **Sobrecosto de estabilidad (contabilizado aparte en `CostoPivote`):** `132` comparaciones para calcular las escalas iniciales $s_i$, `77` divisiones $|a_{ik}|/s_i$, `66` comparaciones de candidatos y `8` intercambios de fila.

---

## 5. Análisis de mantisa y pivoteo

Comparando las dos corridas del mismo código frente a la referencia independiente de `mpmath` (120 dígitos):

| Indicador | Sin pivoteo | Pivoteo parcial escalado |
|:---|---:|---:|
| Multiplicador máximo $\max \|m_{ik}\|$ | $2.7526 \times 10^{13}$ | $2.4656 \times 10^{9}$ |
| Factor de crecimiento $\rho$ | $2641.2243$ | $1.7832$ |
| Error relativo en norma infinito (inicial) | $2.0803 \times 10^{-10}$ | $2.3573 \times 10^{-13}$ |
| Bits de mantisa perdidos $\approx \log_2(e_{\text{rel}} / \varepsilon_{\text{mach}})$ | $\approx 19.8$ bits (de 53) | $\approx 10.1$ bits (de 53) |
| Cifras verificadas por incógnita (inicial) | 7 a 10 | 9 a 13 |

### Por qué el pivote $a_{11} = 1.3407 \times 10^{-8}$ es inaceptable

En el paso 1 sin pivoteo, eliminar $a_{61} = 3.6904 \times 10^{5}$ con $a_{11} = 1.3407 \times 10^{-8}$ produce un multiplicador $m_{61} \approx 2.7526 \times 10^{13}$. Al restar $m_{61} \cdot a_{1j}$ (del orden de $10^{9}$) a $a_{6j}$ (del orden de $10^{5}$), el alineamiento de exponentes en IEEE 754 desplaza y descarta los bits menos significativos de $a_{6j}$. La matriz intermedia crece por un factor $\rho = 2641.22$ y se pierden casi 20 de los 53 bits de mantisa. Con **pivoteo parcial escalado**, el paso 1 intercambia inmediatamente las filas para elegir un pivote de magnitud acorde a la escala de su fila, manteniendo el factor de crecimiento en $\rho = 1.7832$ y reduciendo el error relativo en casi tres órdenes de magnitud.

---

## 6. Condicionamiento

Mediante descomposición en valores singulares (`mp.svd` a 80 dígitos) y verificación matricial se obtuvo:

* $\det(A) = 1.364551 \times 10^{9}$ (sistema no singular).
* $\kappa_2(A) = 3.084740 \times 10^{13}$.
* **Casi dependencia lineal de la fila 7:** proyectando $F_7$ sobre $F_3$ y $F_5$ (`analizar_fila7` en `src/verificacion.py`) se confirma que:
  $$F_7 \approx 4.177250 \cdot F_3 + 583.246740 \cdot F_5$$
  con perturbación residual $\|F_7 - \text{comb}\|_2 = 1.942375 \times 10^{-4}$ (error relativo $2.283352 \times 10^{-3}$).

### Cota clásica vs. resultado real

La cota clásica en doble precisión ($\varepsilon_{\text{mach}} \approx 2.2204 \times 10^{-16}$) predice:

$$\frac{\|\delta x\|}{\|x\|} \lesssim \kappa_2(A) \cdot \varepsilon_{\text{mach}} \approx 6.85 \times 10^{-3} \quad (\approx 2.16 \text{ dígitos garantizados})$$

En la práctica se obtienen **9 a 13 dígitos correctos** (antes de refinar) porque: (1) gran parte de $\kappa_2(A)$ refleja la disparidad de $10^{10}$ entre escalas de fila y el pivoteo escalado es invariante frente a dicho escalamiento, manteniendo pequeño el error hacia atrás por fila; y (2) la cota $\kappa \cdot \varepsilon$ asume el peor alineamiento direccional posible con el vector singular asociado a $\sigma_{\min}$, escenario que no se satura aquí ($\rho = 1.7832$).

---

## 7. Resultados

Valores obtenidos en la corrida directa con **pivoteo parcial escalado** (verificados contra `mpmath` a 80 y 120 dígitos en `resultados/solucion.csv`):

| Incógnita | Valor reportado (inicial) | Cifras verificadas | Valor tras 1 pasada de refinamiento |
|:---|---:|:---:|---:|
| $x_1$ | `0.738491562700` | 12 | `0.7384915627000000` |
| $x_2$ | `-1.20593741680` | 12 | `-1.205937416800000` |
| $x_3$ | `0.004387195200` | 10 | `0.004387195200000000` |
| $x_4$ | `-0.091726453800` | 11 | `-0.09172645380000000` |
| $x_5$ | `2.846109372500` | 13 | `2.8461093725000000` |
| $x_6$ | `-0.00628154090` | 9 | `-0.006281540900000000` |
| $x_7$ | `1.473820591600` | 13 | `1.4738205916000000` |
| $x_8$ | `-0.038492716400` | 11 | `-0.03849271640000000` |
| $x_9$ | `0.915264803700` | 12 | `0.9152648037000000` |
| $x_{10}$ | `-2.607139485200` | 13 | `-2.607139485200000` |
| $x_{11}$ | `0.051873629400` | 11 | `0.05187362940000000` |
| $x_{12}$ | `-1.08493716` | 9 | `-1.084937162500000` |

Todas las incógnitas tienen entre 9 y 13 cifras significativas exactas verificadas en la corrida directa (y 16 a 17 cifras representables en `float64` tras el refinamiento iterativo).

---

## 8. Verificación

* **Residuo en precisión extendida (`mpmath`, 80 dígitos):**
  * Pivoteo parcial escalado (inicial): $\|r\|_\infty = 2.290876 \times 10^{-10}$ (error relativo de la solución $2.3573 \times 10^{-13}$).
  * Sin pivoteo (inicial): $\|r\|_\infty = 1.771752 \times 10^{-7}$ (error relativo $2.0803 \times 10^{-10}$, casi mil veces peor).
* **Refinamiento iterativo (reutilizando $PA = LU$ sin refactorizar):** tras una pasada con el residuo exacto en `mpmath`, el error relativo de la solución cae de $2.3573 \times 10^{-13}$ a $6.0160 \times 10^{-17}$ (límite de `float64`), dejando $\|r\|_\infty = 2.603187 \times 10^{-10}$ (acotado por el redondeo de representar $x$ en `float64` frente a filas de orden $10^5$).
* **Control con biblioteca (`numpy.linalg.solve`):** el error relativo de `numpy.linalg.solve` frente a la referencia exacta es $5.2322 \times 10^{-13}$, ligeramente mayor que nuestra implementación con pivoteo escalado ($2.3573 \times 10^{-13}$) porque LAPACK `dgesv` utiliza pivoteo parcial sin escalamiento previo de filas.

---

## 9. Conclusiones

El compromiso elegido —eliminación gaussiana densa con pivoteo parcial escalado, barrido estricto desde $j = k+1$ y una única factorización $PA = LU$— cumple el mínimo operacional (`1354` flops, $-8.88\%$ bajo la cota de `1486`) y alcanza 9 a 13 dígitos significativos exactos (16+ tras una pasada de refinamiento).

**Qué haríamos distinto con $n = 1000$:**
* Reemplazar los bucles puros de Python por factorización por bloques vectorizada (BLAS-3 / LAPACK) y equilibrado previo de filas/columnas antes de pivotar, evitando el costo cuadrático de reevaluar cocientes de escala en Python interpretado.
* Si la matriz $n = 1000$ fuera dispersa (común en EDPs o redes), migrar a métodos iterativos precondicionados (Gauss-Seidel/SOR o GMRES/CG con ILU) para no destruir la estructura dispersa ni almacenar $\mathcal{O}(n^2)$ entradas densas.

---

## 10. Reproducibilidad

* **Entorno verificado:** Python `3.11.2` (y `3.12.14`), `mpmath==1.3.0`, `numpy==2.1.3`.
* **Comando ejecutado:**
  ```bash
  pip install -r requirements.txt
  python src/main.py
  ```

### Declaración de uso de asistente de IA

De acuerdo con la sección 9 del enunciado, se declara el uso de asistentes de IA (Gemini en la revisión de rúbrica, contraste con NumPy, y redacción/validación del informe). Todo el código fue ejecutado, verificado contra la referencia de `mpmath` a 80/120 dígitos y revisado línea por línea por ambos integrantes

### Referencias

* Burden, R. L., y Faires, J. D. *Análisis Numérico* (9.ª ed.). Cengage Learning, Cap. 6 (Secciones 6.1–6.2).
* Golub, G. H., y Van Loan, C. F. *Matrix Computations* (4.ª ed.). Johns Hopkins University Press, Cap. 3.
* Higham, N. J. *Accuracy and Stability of Numerical Algorithms* (2.ª ed.). SIAM, 2002.
* Trefethen, L. N., y Bau, D. *Numerical Linear Algebra*. SIAM, 1997 (Lecciones 20–22).
* Pesante, P., y Cabezas Contardo, J. M. Apuntes de Cátedra *Computación Numérica (INF-321) — Unidad III*. Universidad Católica del Maule, 2026.
