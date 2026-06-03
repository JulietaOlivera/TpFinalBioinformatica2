"""
restriccion.py
Análisis de sitios de restricción en los amplicones generados.

Las enzimas con sitios degenerados (como AvaII: GGWCC, donde W = A o T)
se manejan expandiendo el patrón con expresiones regulares.
"""

import re


ENZIMAS = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "AvaII": "GG[AT]CC",   # W = A o T
}


def contar_sitios(secuencia, patron):
    """
    Cuenta ocurrencias NO solapadas del patrón en la secuencia.
    Usa re.findall para manejar bases degeneradas.
    """
    return len(re.findall(patron, secuencia, re.IGNORECASE))


def analizar_amplicones(resultados_primers):
    """
    Recorre los amplicones y cuenta los sitios de cada enzima.

    Parámetros
    ----------
    resultados_primers : list[dict]
        Lista de diccionarios retornados por disenar_primers_gen.

    Retorna
    -------
    list[dict] con claves: amplicón, enzima, sitio, cantidad.
    """
    filas = []
    for resultado in resultados_primers:
        nombre    = resultado["nombre"]
        secuencia = resultado["seq_amplicón"]
        for enzima, patron in ENZIMAS.items():
            cantidad = contar_sitios(secuencia, patron)
            filas.append({
                "amplicón": f"Amplicón_{nombre}",
                "enzima":   enzima,
                "sitio":    patron,
                "cantidad": cantidad,
            })
    return filas
