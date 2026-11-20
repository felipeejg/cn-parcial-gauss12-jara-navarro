import numpy as np
from mpmath import mp
from gauss import correccion


def control(matriz_texto, vector_texto):
    # referencia independiente del sistema decimal publicado a 80 y 120 cifras
    referencias = []
    for precision in (80, 120):
        with mp.workdps(precision):
            A, b = mp.matrix(matriz_texto), mp.matrix(vector_texto)
            referencias.append(mp.lu_solve(A, b))
    x = referencias[-1]
    if max(abs(t - s) for t, s in zip(*referencias)) > mp.mpf("1e-60"):
        raise ArithmeticError("Las referencias no coinciden a la precisión prevista")
    singulares = mp.svd(A, compute_uv=False)
    return A, b, x, singulares[0] / singulares[len(singulares) - 1]


def control_biblioteca(A64, b64):
    # control independiente mediante numpy.linalg.solve
    return np.linalg.solve(np.array(A64, dtype=np.float64), np.array(b64, dtype=np.float64))


def analizar_fila7(A64):
    # cuantifica la casi dependencia lineal de la fila 7 respecto a las filas 3 y 5
    A_np = np.array(A64, dtype=np.float64)
    r3, r5, r7 = A_np[2], A_np[4], A_np[6]
    coefs, _, _, _ = np.linalg.lstsq(np.column_stack([r3, r5]), r7, rcond=None)
    pert = r7 - (coefs[0] * r3 + coefs[1] * r5)
    norma_pert = float(np.linalg.norm(pert, 2))
    norma_rel = float(norma_pert / np.linalg.norm(r7, 2))
    return float(coefs[0]), float(coefs[1]), norma_pert, norma_rel, float(np.linalg.det(A_np))


def evaluar(A, b, x, referencia):
    aproximacion = mp.matrix(x)
    r = b - A * aproximacion
    error = aproximacion - referencia
    return r, max(map(abs, r)), max(map(abs, error)) / max(map(abs, referencia))


def refinar(factores, residuo):
    delta = correccion(factores, residuo)
    return [x + d for x, d in zip(factores.solucion, delta)]


def cifras(x, referencia):
    # mayor prefijo de precisiones 1..17 con igual redondeo significativo
    for d in range(1, 18):
        if mp.nstr(mp.mpf(x), d) != mp.nstr(referencia, d):
            return d - 1
    return 17
