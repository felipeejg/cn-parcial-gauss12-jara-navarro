# El pivote · INF-321

Eliminación gaussiana densa de 12 incógnitas con mínimo costo operacional y máxima exactitud de mantisa, comparando ejecución sin pivoteo y con pivoteo parcial escalado.

## Autoría y asignatura

* **Estudiantes:** Felipe Jara Guzmán · Kevin Navarro
* **Asignatura:** Computación Numérica (INF-321) — Sección S2
* **Institución:** Ingeniería Civil Informática · Universidad Católica del Maule (Talca, noviembre de 2026)

## Requisitos

* Python 3.10 o superior (probado en Python 3.11.2 y 3.12.14).
* Instalación de dependencias:

```bash
pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repositorio:

```bash
python src/main.py
```

El comando regenera automáticamente todos los archivos en `resultados/`:
* `resultados/solucion.csv` — incógnitas (sin pivoteo y pivoteo parcial escalado, inicial y refinada) con cifras verificadas y error relativo.
* `resultados/traza_eliminacion.txt` — matriz aumentada tras cada uno de los 11 pasos en ambas estrategias.
* `resultados/conteo_flops.md` — desglose paso a paso de flops y sobrecosto de pivoteo.
* `resultados/residuo.txt` — residuo en precisión extendida (`mpmath`), norma infinito, condicionamiento, fila 7 y control con NumPy.

## Estructura de carpetas

```text
cn-parcial-gauss12-jara-navarro/
├─ datos/         # Matriz A (12×12) y vector b (20 dígitos) sin truncar en formato CSV
├─ src/           # Código fuente: eliminación, contador de flops, verificación y punto de entrada
├─ resultados/    # Salidas generadas automáticamente al ejecutar python src/main.py
├─ informe.md     # Informe técnico completo (10 secciones)
├─ readme.md      # Guía de instalación, ejecución y resultado resumido
└─ requirements.txt
```

## Resultado resumido

Corrida con **pivoteo parcial escalado** y una pasada de refinamiento iterativo:

| Incógnita | Valor a 10 dígitos | Incógnita | Valor a 10 dígitos |
|:---|---:|:---|---:|
| $x_1$ | `0.7384915627` | $x_7$ | `1.473820592` |
| $x_2$ | `-1.205937417` | $x_8$ | `-0.03849271640` |
| $x_3$ | `0.004387195200` | $x_9$ | `0.9152648037` |
| $x_4$ | `-0.09172645380` | $x_{10}$ | `-2.607139485` |
| $x_5$ | `2.846109372` | $x_{11}$ | `0.05187362940` |
| $x_6$ | `-0.006281540900` | $x_{12}$ | `-1.084937162` |

* **Total de flops aritméticos:** `1354` (`716` mult/div + `638` sum/res; sobrecosto de pivoteo escalado contabilizado aparte).
* **Norma infinito del residuo:** $\|r\|_\infty = 2.290876 \times 10^{-10}$ (inicial con pivoteo escalado) y $2.603187 \times 10^{-10}$ (refinada).

## Informe técnico

El análisis completo de estrategia, conteo paso a paso, pérdida de mantisa, condicionamiento y verificación se encuentra en [informe.md](informe.md).

