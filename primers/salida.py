"""
salida.py
Funciones para escribir los archivos de salida del pipeline.
"""


def escribir_tabla_primers(resultados, ruta_salida):
    """
    Escribe el archivo TSV con la información de cada par de primers.

    Columnas: Nombre; seq_fw; Tm_fw; Ta_fw; seq_rv; Tm_rv; Ta_rv; Tamaño_producto
    """
    encabezado = "\t".join([
        "Nombre_gen", "Seq_primer_fw", "Tm_fw", "Ta_fw",
        "Seq_primer_rv", "Tm_rv", "Ta_rv", "Tamaño_producto_pb"
    ])
    with open(ruta_salida, "w") as archivo:
        archivo.write(encabezado + "\n")
        for r in resultados:
            fila = "\t".join([
                r["nombre"],
                r["seq_fw"],
                str(r["tm_fw"]),
                str(r["ta_fw"]),
                r["seq_rv"],
                str(r["tm_rv"]),
                str(r["ta_rv"]),
                str(r["tamaño"]),
            ])
            archivo.write(fila + "\n")


def escribir_multifasta_amplicones(resultados, ruta_salida):
    """
    Escribe un multifasta con la secuencia de cada amplicón.

    Formato del identificador: >Amplicón_[nombre_gen]_[tamaño]
    """
    with open(ruta_salida, "w") as archivo:
        for r in resultados:
            header = f">Amplicón_{r['nombre']}_{r['tamaño']}"
            archivo.write(header + "\n")
            archivo.write(r["seq_amplicón"] + "\n")


def escribir_tabla_restriccion(filas, ruta_salida):
    """
    Escribe el archivo TSV con los resultados del análisis de restricción.

    Columnas: Amplicón; Enzima; Sitio; Cantidad
    """
    encabezado = "\t".join(["Amplicón", "Enzima", "Sitio", "Cantidad"])
    with open(ruta_salida, "w") as archivo:
        archivo.write(encabezado + "\n")
        for fila in filas:
            linea = "\t".join([
                fila["amplicón"],
                fila["enzima"],
                fila["sitio"],
                str(fila["cantidad"]),
            ])
            archivo.write(linea + "\n")
