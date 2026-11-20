import csv
from pathlib import Path
from mpmath import mp
from gauss import resolver
from verificacion import analizar_fila7, cifras, control, control_biblioteca, evaluar, refinar

RAIZ = Path(__file__).resolve().parents[1]


def leer_csv(nombre):
    with (RAIZ / "datos" / nombre).open(encoding="utf-8", newline="") as archivo:
        return list(csv.reader(archivo))


def escribir_paso(archivo, k, p, U, c):
    print(f"paso {k+1} | fila pivote antes del intercambio: {p+1}", file=archivo)
    for fila, termino in zip(U, c):
        print(" ".join(map(repr, fila)), "|", repr(termino), file=archivo)


def tabla_costo(archivo, factores):
    print("| paso | mult/div | sumas/restas | subtotal |", file=archivo)
    print("|---|---:|---:|---:|", file=archivo)
    filas = list(enumerate(factores.cuentas, 1)) + [('regresiva', factores.regresiva)]
    for paso, cuenta in filas:
        print(f"| {paso} | {cuenta.md} | {cuenta.sr} | {cuenta.md + cuenta.sr} |", file=archivo)
    md = sum(c.md for _, c in filas)
    sr = sum(c.sr for _, c in filas)
    total = md + sr
    desvio = ((total - 1486) / 1486) * 100.0
    print(f"| total | {md} | {sr} | {total} |\n", file=archivo)
    print(f"total aritmético principal = {total} flops "
          f"({desvio:+.2f}% respecto a la cota 1486; ahorro de 132 flops por no operar los 66 ceros subdiagonales)",
          file=archivo)
    print(f"\nsobrecosto separado: {factores.pivote}.\n", file=archivo)


def main():
    mp.dps = 80
    matriz = leer_csv("matriz_a.csv")
    vector = leer_csv("vector_b.csv")
    if len(matriz) != 12 or any(len(f) != 12 for f in matriz):
        raise ValueError("Se esperaban 12 filas y 12 columnas")
    if len(vector) != 12 or any(len(f) != 1 for f in vector):
        raise ValueError("se esperaban 12 términos independientes")
    terminos = [f[0] for f in vector]
    A, b, referencia, condicion = control(matriz, terminos)
    A64 = [list(map(float, fila)) for fila in matriz]
    b64 = list(map(float, terminos))
    x_np = control_biblioteca(A64, b64)
    c3, c5, pert_f7, rel_f7, det_a = analizar_fila7(A64)

    carpeta = RAIZ / "resultados"
    carpeta.mkdir(exist_ok=True)
    with ((carpeta / "traza_eliminacion.txt").open("w", encoding="utf-8") as traza,
          (carpeta / "conteo_flops.md").open("w", encoding="utf-8") as costos,
          (carpeta / "residuo.txt").open("w", encoding="utf-8") as residuos,
          (carpeta / "solucion.csv").open("w", encoding="utf-8", newline="") as salida):
        print("una flop por suma, resta, multiplicación o división del núcleo "
              "Las copias, asignaciones y mediciones se excluyen\n", file=costos)
        print(f"det_A = {det_a:.6e}", file=residuos)
        print("kappa_2 =", mp.nstr(condicion, 18), file=residuos)
        print(f"fila7_aprox = {c3:.6f} * F3 + {c5:.6f} * F5 | "
              f"||pert||_2 = {pert_f7:.6e} (rel = {rel_f7:.6e})", file=residuos)
        csv_salida = csv.writer(salida)
        csv_salida.writerow(["estrategia", "estado", "i", "x_reportada", "cifras_verificadas",
                             "x_binary64_repr", "error_absoluto", "error_relativo"])
        for nombre, escalado in [("sin_pivoteo", False), ("parcial_escalado", True)]:
            print(nombre, file=traza)
            factores = resolver(A64, b64, escalado,
                                lambda k, p, U, c: escribir_paso(traza, k, p, U, c))
            print("##", nombre, "\n", file=costos)
            tabla_costo(costos, factores)
            print(nombre, "max_m =", factores.max_m, "crecimiento =", factores.crecimiento,
                  file=residuos)
            x = factores.solucion
            for estado in ("inicial", "refinada"):
                r, norma, error = evaluar(A, b, x, referencia)
                print(estado, "r =", [mp.nstr(v, 20) for v in r], file=residuos)
                print("norma_inf =", mp.nstr(norma, 20), "error_relativo_inf =",
                      mp.nstr(error, 20), file=residuos)
                for i, (valor, ref) in enumerate(zip(x, referencia), 1):
                    d = cifras(valor, ref)
                    err = abs(mp.mpf(valor) - ref)
                    csv_salida.writerow([nombre, estado, i,
                                         mp.nstr(mp.mpf(valor), d, strip_zeros=False) if d else "",
                                         d, repr(valor), mp.nstr(err, 18),
                                         mp.nstr(err / abs(ref), 18) if ref else ""])
                if estado == "inicial":
                    x = refinar(factores, r)

        _, _, err_np = evaluar(A, b, list(map(float, x_np)), referencia)
        print("control_numpy error_relativo_inf =", mp.nstr(err_np, 20), file=residuos)


if __name__ == "__main__":
    main()
