

import re #modulo para trabajr con expresiones regulares


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

            seqid, _, feature, start, end, _, hebra, _, atributos = campos #dada la forma de un archivo .GFF, hay campos que no resultan relevantes para la consigna y se pasan por alto

            if feature not in ("gene", "CDS"): #Prokka da resultados de CDSs no de Genes, son los unicos resultados que vamos  a poder analizar 
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
        match = re.search(rf"{clave}=([^;]+)", atributos) #search busca la primera coincidencia, r (rawstring para que no se interpreten carácteres especiales de python)
        """Los paréntesis indican que se capture el patron dentro de ellos. cualquier caracter (^) excepto
        ; + significa 1 o mas veces para que vaya leyendo todo lo correspondiente al nombre del gen deteniendose en el ;"""
        if match:
            return match.group(1) #match.group(0) devolveria todo lo capturado incluyendo ID, (1) solo el id o nombre limpio
    return atributos.split(";")[0] #si no funciona va a devolver al menos el primer fragmento de información. 
