import sys #guarda argumentos de línea de comandos (ingresar los archivos de entrada)
from Bio import SeqIO #para leer archivos FASTA. generará objetos SeqRecord

from parseo_gff      import parsear_gff
from diseno_primers  import disenar_primers_gen
from restriccion     import analizar_amplicones
from salida          import (
    escribir_tabla_primers,
    escribir_multifasta_amplicones,
    escribir_tabla_restriccion,
)


#Archivos de salida
SALIDA_PRIMERS     = "primers.tab"
SALIDA_AMPLICONES  = "amplicones.fasta"
SALIDA_RESTRICCION = "res_enzimas.tab"


def cargar_genoma(ruta_fasta):
    """
    Carga el FASTA como secuencia genómica.
    """
    registros = list(SeqIO.parse(ruta_fasta, "fasta")) #parse devuelve un iterador (Eficiencia de memoria) que luego se hace una lista para verificar cuantas secuencias hay
    if not registros:
        raise ValueError(f"No se encontraron secuencias en {ruta_fasta}")
    if len(registros) > 1:
        print(f"[AVISO] El FASTA contiene {len(registros)} secuencias. "
              "Se usará la primera.")
    return registros[0].seq


def main():
    if len(sys.argv) != 3:
        print("Uso: python diseño_kit_diagnostico.py <ensamblado.fasta> <genes.gff>")
        sys.exit(1)

    ruta_fasta = sys.argv[1]
    ruta_gff   = sys.argv[2]

    #Se cargan los datos de entrada
    print("Cargando genoma...")
    secuencia_genomica = cargar_genoma(ruta_fasta)

    print("Parseando GFF...")
    genes = parsear_gff(ruta_gff)
    print(f"  → {len(genes)} genes encontrados.")

    # Se diseñan los primers
    print("Diseñando primers...")
    resultados = []
    omitidos   = 0
    for gen in genes:
        resultado = disenar_primers_gen(gen, secuencia_genomica)
        if resultado is None:
            print(f"  [OMITIDO] {gen['nombre']}: coordenadas fuera de rango.")
            omitidos += 1
        else:
            resultados.append(resultado)
    print(f"  → {len(resultados)} pares de primers generados ({omitidos} omitidos).")

    #Se escribe tabla de primers y multifasta de amplicones
    escribir_tabla_primers(resultados, SALIDA_PRIMERS)
    print(f"Tabla de primers guardada en: {SALIDA_PRIMERS}")

    escribir_multifasta_amplicones(resultados, SALIDA_AMPLICONES)
    print(f"Multifasta de amplicones guardado en: {SALIDA_AMPLICONES}")

    #Análisis de sitios de restricción
    print("Analizando sitios de restricción...")
    filas_restriccion = analizar_amplicones(resultados)
    escribir_tabla_restriccion(filas_restriccion, SALIDA_RESTRICCION)
    print(f"Tabla de restricción guardada en: {SALIDA_RESTRICCION}")

    print("\n¡Pipeline completado exitosamente!")


if __name__ == "__main__":
    main()
