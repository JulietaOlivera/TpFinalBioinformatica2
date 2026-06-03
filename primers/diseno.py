"""
diseno_primers.py
Diseño de primers forward y reverse de 20 nt para cada gen predicho.

Convención de posiciones:
    - Los primers "pegan" a 5 nt del inicio/fin del gen.
    - Primer FW: comienza en (inicio_gen - 5), longitud 20 nt, hebra directa.
    - Primer RV: termina en  (fin_gen   + 5), longitud 20 nt, complementario inverso.
    - Para genes en hebra negativa los roles se intercambian respecto a la secuencia
      genómica, pero la lógica de extracción es la misma; el amplicón siempre se
      reporta en 5'→3' de la hebra codificante.
"""

from Bio.SeqUtils import MeltingTemp as mt


LARGO_PRIMER    = 20
OFFSET_EXTREMO  =  5


def calcular_tm(secuencia):
    """Calcula la Tm con el método de Wallace (regla 2/4)."""
    return round(mt.Tm_Wallace(secuencia), 2)


def calcular_ta(tm):
    """Temperatura de annealing: Tm - 5 °C."""
    return round(tm - 5, 2)


def disenar_primers_gen(gen, secuencia_genomica):
    """
    Diseña el par de primers para un gen dado.

    Parámetros
    ----------
    gen : dict
        Diccionario con claves 'nombre', 'inicio', 'fin', 'hebra'.
    secuencia_genomica : Bio.Seq
        Secuencia completa del genoma ensamblado.

    Retorna
    -------
    dict con primers, Tm, Ta y datos del amplicón, o None si las coordenadas
    quedan fuera de rango.
    """
    inicio_gen = gen["inicio"]
    fin_gen    = gen["fin"]
    hebra      = gen["hebra"]

    # Coordenadas de los primers sobre la hebra directa del genoma
    inicio_fw = inicio_gen - OFFSET_EXTREMO
    fin_rv     = fin_gen   + OFFSET_EXTREMO

    if inicio_fw < 0 or fin_rv > len(secuencia_genomica):
        return None

    seq_fw_raw = secuencia_genomica[inicio_fw : inicio_fw + LARGO_PRIMER]
    seq_rv_raw = secuencia_genomica[fin_rv - LARGO_PRIMER : fin_rv]

    if hebra == "+":
        primer_fw = seq_fw_raw
        primer_rv = seq_rv_raw.reverse_complement()
    else:
        # En hebra negativa el primer que "pega" al inicio del gen
        # es el complementario inverso del extremo 3' del amplicón genómico
        primer_fw = seq_rv_raw.reverse_complement()
        primer_rv = seq_fw_raw

    tm_fw = calcular_tm(primer_fw)
    tm_rv = calcular_tm(primer_rv)

    amplicón = secuencia_genomica[inicio_fw:fin_rv]
    if hebra == "-":
        amplicón = amplicón.reverse_complement()

    return {
        "nombre":       gen["nombre"],
        "seq_fw":       str(primer_fw),
        "tm_fw":        tm_fw,
        "ta_fw":        calcular_ta(tm_fw),
        "seq_rv":       str(primer_rv),
        "tm_rv":        tm_rv,
        "ta_rv":        calcular_ta(tm_rv),
        "tamaño":       len(amplicón),
        "seq_amplicón": str(amplicón),
    }

