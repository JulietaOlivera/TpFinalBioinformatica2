from Bio.SeqUtils import MeltingTemp as mt


LARGO_PRIMER = 20   #consigna
OFFSET_EXTREMO = 5


def calcular_tm(secuencia):
    """Calcula la Tm con el método de Wallace."""
    return round(mt.Tm_Wallace(secuencia), 2) #calcula con regla de Wallace 2(A+T)+4(G+C)


def calcular_ta(tm):
    """Temperatura de annealing: Tm - 5 °C."""
    return round(tm - 5, 2)


def disenar_primers_gen(gen, secuencia_genomica):
    """
    Diseña el par de primers para un gen dado.
    """
    inicio_gen = gen["inicio"] #la lista genes se definió en parseo_gff.py
    fin_gen    = gen["fin"]
    hebra      = gen["hebra"]

    # Coordenadas de los primers sobre la hebra directa del genoma
    inicio_fw = inicio_gen - OFFSET_EXTREMO
    fin_rv = fin_gen   + OFFSET_EXTREMO

    if inicio_fw < 0 or fin_rv > len(secuencia_genomica):
        return None

    seq_fw_raw = secuencia_genomica[inicio_fw : inicio_fw + LARGO_PRIMER] #forward
    seq_rv_raw = secuencia_genomica[fin_rv - LARGO_PRIMER : fin_rv] #reverse

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

    amplicón = secuencia_genomica[inicio_fw:fin_rv] #incluye los 5 nt antes
    if hebra == "-":
        amplicón = amplicón.reverse_complement() #para reportar el producto en la orientación correcta

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

