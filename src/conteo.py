from dataclasses import dataclass


@dataclass
class Cuenta:
    #Registra la operación al ejecutarla, incluida la sustitución regresiva.
    md: int = 0
    sr: int = 0

    def dividir(self, a, b):
        self.md += 1
        return a / b

    def actualizar(self, valor, multiplicador, coeficiente):
        self.md += 1
        self.sr += 1
        return valor - multiplicador * coeficiente


@dataclass
class CostoPivote:
    comparaciones_escalas: int = 0
    comparaciones_candidatos: int = 0
    cocientes: int = 0
    intercambios: int = 0
