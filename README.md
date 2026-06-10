# Trabajo Práctico Final Integrador Análisis y Alineamiento de Secuencias
# Kit Diagnóstico PCR — Pipeline Bioinformático

Pipeline en Python para el diseño automatizado de primers de PCR a partir de un ensamblado genómico y un archivo de predicción de genes (GFF3).


## Requisitos

- Python ≥ 3.8
- [Biopython](https://biopython.org/)

```bash
pip install biopython
```


## Uso

```bash
python diseño_kit_diagnostico.py <ensamblado.fasta> <genes.gff>
```

Argumentos:
"ensamblado.fasta" : Secuencia genómica resultante del ensamblado (FASTA) 
"genes.gff": Archivo GFF3 con los genes predichos

Ejemplo:

```bash
python diseño_kit_diagnostico.py genoma_hibrido.fasta genes_predichos.gff
```


## Archivos de salida

**"primers.tab"**
Tabla separada por tabulaciones con un par de primers por gen.

"Columna" : Explicación
"Nombre_gen" : Identificador del gen
"Seq_primer_fw" : Secuencia del primer forward
"Tm_fw" : Temperatura de melting del primer forward
"Ta_fw" : Temperatura de annealing del fw
"Seq_primer_rv" : Secuencia del primer reverse
"Tm_rv" : Temperatura de melting del primer reverse
"Ta_rv" : Temperatura de annealing del rv
"Tamaño_producto_pb" : Longitud del amplicón en pares de bases

**"amplicones.fasta"**
Multifasta con la secuencia de cada producto de PCR.

```
>Amplicón_[nombre_gen]_[tamaño_pb]
ATCGATCG...
```

**"res_enzimas.tab"**
Tabla separada por tabulaciones con el recuento de sitios de restricción por amplicón y enzima.

"Columna" : Explicación
"Amplicón" : Nombre del amplicón (`Amplicón_[nombre_gen]`)
"Enzima" : Nombre de la enzima (EcoRI, BamHI, AvaII)
"Sitio" : Patrón de reconocimiento
"Cantidad" : Número de veces que aparece el sitio en el amplicón

---

## Lógica de diseño de primers

- Los primers tienen 20 nt de longitud.
- Se posicionan a 5 nt del inicio y fin del gen, de modo que el amplicón incluye la secuencia codificante completa.
- Para genes en hebra negativa (−), el primer forward corresponde al complementario inverso del extremo 3' en la secuencia genómica, y viceversa.

### Temperatura de melting

Se usa la **regla de Wallace** (Tm = 2·(A+T) + 4·(G+C)), implementada por `Bio.SeqUtils.MeltingTemp.Tm_Wallace`.



## Enzimas de restricción analizadas

| Enzima | Sitio de reconocimiento| Patrón regex |
|--------|------------------------|--------------|
| EcoRI  | GAATTC                 | `GAATTC`     |
| BamHI  | GGATCC                 | `GGATCC`     |
| AvaII  | GGWCC (W = A o T)      | `GG[AT]CC`   |

