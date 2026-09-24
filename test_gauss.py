import unittest
from fractions import Fraction

from gauss import (DETERMINADO, INCOMPATIBLE, INDETERMINADO, formato_solucion,
                   resolver_gauss)


def evaluar(expr, parametros):
    """Evalúa una expresión afín {j: coef, None: cte} con valores de parámetros."""
    return expr[None] + sum(c * parametros[j] for j, c in expr.items() if j is not None)


def comprobar_solucion(test, matriz, resultado, valores_parametros):
    n = len(matriz[0]) - 1
    params = dict(zip(resultado["libres"], valores_parametros))
    x = [evaluar(e, params) for e in resultado["solucion"]]
    for fila in matriz:
        izquierda = sum(Fraction(fila[j]) * x[j] for j in range(n))
        test.assertEqual(izquierda, Fraction(fila[n]))


class TestGauss(unittest.TestCase):
    def test_solucion_unica(self):
        matriz = [[2, 1, -1, 8],
                  [-3, -1, 2, -11],
                  [-2, 1, 2, -3]]
        r = resolver_gauss(matriz)
        self.assertEqual(r["tipo"], DETERMINADO)
        self.assertEqual([e[None] for e in r["solucion"]], [2, 3, -1])

    def test_solucion_fraccionaria(self):
        r = resolver_gauss([[3, 2, 1], [1, 4, 2]])
        self.assertEqual(r["tipo"], DETERMINADO)
        self.assertEqual([e[None] for e in r["solucion"]],
                         [Fraction(0), Fraction(1, 2)])

    def test_pivote_cero_requiere_intercambio(self):
        r = resolver_gauss([[0, 1, 2], [1, 0, 3]])
        self.assertEqual([e[None] for e in r["solucion"]], [3, 2])

    def test_menos_ecuaciones_que_variables(self):
        matriz = [[1, 2, -1, 3],
                  [2, 4, 1, 3]]
        r = resolver_gauss(matriz)
        self.assertEqual(r["tipo"], INDETERMINADO)
        self.assertEqual(r["rango"], 2)
        self.assertEqual(len(r["libres"]), 1)
        for t in (0, 1, -5, Fraction(7, 3)):
            comprobar_solucion(self, matriz, r, [t])

    def test_ecuaciones_dependientes(self):
        matriz = [[1, 1, 1, 6],
                  [2, 2, 2, 12],
                  [1, -1, 0, 0]]
        r = resolver_gauss(matriz)
        self.assertEqual(r["tipo"], INDETERMINADO)
        self.assertEqual(r["rango"], 2)
        for t in (0, 3, -2):
            comprobar_solucion(self, matriz, r, [t])

    def test_dos_parametros(self):
        matriz = [[1, 2, 3, 4, 10]]
        r = resolver_gauss(matriz)
        self.assertEqual(r["tipo"], INDETERMINADO)
        self.assertEqual(len(r["libres"]), 3)
        comprobar_solucion(self, matriz, r, [1, 2, 3])
        self.assertIn("x1 = 10 - 2·t1 - 3·t2 - 4·t3", formato_solucion(r))

    def test_incompatible(self):
        r = resolver_gauss([[1, 1, 2], [1, 1, 3]])
        self.assertEqual(r["tipo"], INCOMPATIBLE)
        self.assertIsNone(r["solucion"])

    def test_decimales_y_fracciones_texto(self):
        r = resolver_gauss([["0.5", "1/2", "1"], [1, -1, 0]])
        self.assertEqual([e[None] for e in r["solucion"]], [1, 1])

    def test_rechaza_mas_ecuaciones_que_variables(self):
        with self.assertRaises(ValueError):
            resolver_gauss([[1, 2], [3, 4]])

    def test_rechaza_filas_irregulares(self):
        with self.assertRaises(ValueError):
            resolver_gauss([[1, 2, 3], [1, 2]])


if __name__ == "__main__":
    unittest.main()
