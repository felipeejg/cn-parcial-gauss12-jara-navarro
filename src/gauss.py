from dataclasses import dataclass
from conteo import Cuenta, CostoPivote


@dataclass
class Factorizacion:
    L: list
    U: list
    orden: list
    solucion: list
    cuentas: list
    regresiva: Cuenta
    pivote: CostoPivote
    max_m: float
    crecimiento: float


def escalas_filas(A, costo):
    escalas = []
    for fila in A:
        mayor = abs(fila[0])
        for coef in fila[1:]:
            costo.comparaciones_escalas += 1
            mayor = max(mayor, abs(coef))
        if mayor == 0:
            raise ValueError('Fila nula en la matriz')
        escalas.append(mayor)
    return escalas


def elegir_pivote(U, escalas, k, costo):
    elegido = k
    mejor = abs(U[k][k]) / escalas[k]
    costo.cocientes += 1
    for i in range(k + 1, len(U)):
        razon = abs(U[i][k]) / escalas[i]
        costo.cocientes += 1
        costo.comparaciones_candidatos += 1
        if razon > mejor:
            elegido, mejor = i, razon
    return elegido


def regresiva(U, c, cuenta):
    x = c.copy()
    for i in reversed(range(len(c))):
        for j in range(i + 1, len(c)):
            x[i] = cuenta.actualizar(x[i], U[i][j], x[j])
        if U[i][i] == 0:
            raise ValueError(f'Pivote cero en la fila {i + 1}')
        x[i] = cuenta.dividir(x[i], U[i][i])
    return x


def resolver(A, b, escalado, registrar=None):
    #Resuelve Ax=b sin modificar los datos. registrar recibe cada paso terminado.
    n = len(b)
    if not n or len(A) != n or any(len(fila) != n for fila in A):
        raise ValueError('Se requiere una matriz cuadrada compatible con b')
    U = [fila.copy() for fila in A]
    L = [[float(i == j) for j in range(n)] for i in range(n)]
    c, orden = b.copy(), list(range(n))
    costo = CostoPivote()
    escalas = escalas_filas(A, costo) if escalado else []
    inicial = max(abs(v) for fila in A for v in fila)
    pico, max_m = inicial, 0.0
    cuentas = []
    for k in range(n - 1):
        cuenta = Cuenta()
        cuentas.append(cuenta)
        p = elegir_pivote(U, escalas, k, costo) if escalado else k
        if p != k:
            U[k], U[p] = U[p], U[k]
            c[k], c[p] = c[p], c[k]
            orden[k], orden[p] = orden[p], orden[k]
            escalas[k], escalas[p] = escalas[p], escalas[k]
            L[k][:k], L[p][:k] = L[p][:k], L[k][:k]
            costo.intercambios += 1
        if U[k][k] == 0:
            raise ValueError(f'Pivote cero en el paso {k + 1}')
        for i in range(k + 1, n):
            m = cuenta.dividir(U[i][k], U[k][k])
            L[i][k] = m
            max_m = max(max_m, abs(m))
            U[i][k] = 0.0
            for j in range(k + 1, n):
                U[i][j] = cuenta.actualizar(U[i][j], m, U[k][j])
                pico = max(pico, abs(U[i][j]))
            c[i] = cuenta.actualizar(c[i], m, c[k])
        if registrar:
            registrar(k, p, U, c)
    final = Cuenta()
    x = regresiva(U, c, final)
    return Factorizacion(L, U, orden, x, cuentas, final, costo, max_m, pico / inicial)


def correccion(factores, residuo):
    #Aplica P al residuo y resuelve Ly=Pr, U delta=y con los mismos factores.
    y = [float(residuo[i]) for i in factores.orden]
    for i in range(len(y)):
        for j in range(i):
            y[i] -= factores.L[i][j] * y[j]
    return regresiva(factores.U, y, Cuenta())
