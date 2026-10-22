"""Genera los iconos PNG de la app (sin dependencias externas).

Dibuja |A|: dos barras de determinante y una matriz 2x2 cuyo pivote va en ámbar.
Uso: python3 herramientas/generar_iconos.py
"""
import os
import struct
import zlib

FONDO = (45, 70, 196)
BLANCO = (255, 255, 255)
AMBAR = (246, 196, 69)
SUAVE = (196, 206, 245)


def png(ruta, tam):
    pix = [[FONDO] * tam for _ in range(tam)]

    def rect(x0, y0, x1, y1, color, radio=0.0):
        X0, Y0, X1, Y1 = (round(v * tam) for v in (x0, y0, x1, y1))
        r = radio * tam
        for y in range(Y0, Y1):
            for x in range(X0, X1):
                cx = min(max(x + 0.5, X0 + r), X1 - r)
                cy = min(max(y + 0.5, Y0 + r), Y1 - r)
                if r == 0 or (x + 0.5 - cx) ** 2 + (y + 0.5 - cy) ** 2 <= r * r:
                    pix[y][x] = color

    # Barras del determinante
    rect(0.22, 0.22, 0.27, 0.78, BLANCO, 0.02)
    rect(0.73, 0.22, 0.78, 0.78, BLANCO, 0.02)
    # Matriz 2x2 (el pivote en ámbar)
    rect(0.335, 0.30, 0.475, 0.44, AMBAR, 0.03)
    rect(0.525, 0.30, 0.665, 0.44, BLANCO, 0.03)
    rect(0.335, 0.56, 0.475, 0.70, SUAVE, 0.03)
    rect(0.525, 0.56, 0.665, 0.70, BLANCO, 0.03)

    crudo = b"".join(b"\x00" + bytes(c for p in fila for c in p) for fila in pix)

    def bloque(tipo, datos):
        return struct.pack(">I", len(datos)) + tipo + datos + struct.pack(">I", zlib.crc32(tipo + datos) & 0xFFFFFFFF)

    with open(ruta, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(bloque(b"IHDR", struct.pack(">IIBBBBB", tam, tam, 8, 2, 0, 0, 0)))
        f.write(bloque(b"IDAT", zlib.compress(crudo, 9)))
        f.write(bloque(b"IEND", b""))


if __name__ == "__main__":
    carpeta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "app", "icons")
    for tam in (180, 192, 512):
        png(os.path.join(carpeta, f"icon-{tam}.png"), tam)
