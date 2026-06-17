"""
Parseo de archivos GFF3 para extraer información de genes predichos.
"""

import re


def parsear_gff(ruta_gff):
    """
    Lee un archivo GFF3 y devuelve una lista de diccionarios con la info de cada gen.

    Cada diccionario tiene: nombre, seqid, inicio (0-based), fin (0-based), hebra.
    Solo se procesan líneas con feature 'gene'.
    """
    genes = []

    with open(ruta_gff) as archivo:
        for linea in archivo:
            if linea.startswith("#") or linea.strip() == "":
                continue

            campos = linea.strip().split("\t")
            if len(campos) < 9:
                continue

            seqid, _, feature, start, end, _, hebra, _, atributos = campos

            if feature not in ("gene", "CDS"):
                continue

            nombre = _extraer_nombre(atributos)
            genes.append({
                "nombre":  nombre,
                "seqid":   seqid,
                "inicio":  int(start) - 1,   # GFF es 1-based → convertir a 0-based
                "fin":     int(end),          # el fin en GFF es inclusivo; en slice de Python queda correcto
                "hebra":   hebra
            })

    return genes


def _extraer_nombre(atributos):
    """Extrae el valor de ID o Name del campo de atributos GFF."""
    for clave in ("ID", "Name"):
        match = re.search(rf"{clave}=([^;]+)", atributos)
        if match:
            return match.group(1)
    return atributos.split(";")[0]
