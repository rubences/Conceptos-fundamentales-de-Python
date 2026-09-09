# -*- coding: utf-8 -*-
"""
MUBIO07 - Programación en Python
Actividad 1: Conceptos fundamentales de Python

Autor: Rubén Juárez Cádiz

El script resuelve los tres bloques de la actividad:
1) Preparación del entorno Conda/Spyder (comandos documentados abajo).
2) Manipulación del fichero PDB 1TUP y análisis de la secuencia de p53.
3) Manipulación de los conjuntos de datos con Pandas y visualización.

IMPORTANTE
----------
Los ficheros 1tup.pdb, actividad.csv y ciudades.tsv deben estar en la misma
carpeta que este script. Las figuras se guardan en la subcarpeta "resultados".

===============================================================================
1. PREPARACIÓN DEL ENTORNO DE TRABAJO (2 puntos)
===============================================================================

A) Crear el entorno "actividad1" e instalar Pandas, Matplotlib, Seaborn y Spyder.
   Se fija Python 3.12 para trabajar con una versión estable y ampliamente
   compatible con el ecosistema científico de Conda:

    conda create -n actividad1 -c conda-forge python=3.12 pandas matplotlib seaborn spyder
    conda activate actividad1

B) Instalar/actualizar Spyder a la versión más reciente compatible con el Python
   disponible en el entorno. Conda resuelve automáticamente las dependencias:

    python --version
    conda update -c conda-forge spyder
    spyder --version

C) Comprobar los entornos creados y las librerías del entorno actividad1:

    conda env list
    conda list -n actividad1

   Si el entorno ya está activado, también es suficiente:

    conda list

===============================================================================
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# -----------------------------------------------------------------------------
# CONFIGURACIÓN GENERAL
# -----------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
PDB_FILE = BASE_DIR / "1tup.pdb"
ACTIVIDAD_FILE = BASE_DIR / "actividad.csv"
CIUDADES_FILE = BASE_DIR / "ciudades.tsv"
RESULTADOS_DIR = BASE_DIR / "resultados"
RESULTADOS_DIR.mkdir(exist_ok=True)

# Los 20 aminoácidos proteicos estándar, en código de tres letras.
# Este conjunto se utiliza exclusivamente como filtro para impedir que entren
# nucleótidos (DA, DT, DG, DC) u otros tokens del registro SEQRES.
AMINOACIDOS_ESTANDAR = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

sns.set_theme(style="whitegrid", context="notebook")


# -----------------------------------------------------------------------------
# 2. MANIPULACIÓN DE ARCHIVOS DE TEXTO PLANO (3,5 puntos)
# -----------------------------------------------------------------------------
def leer_metadatos_pdb(ruta_pdb: Path) -> tuple[str, str]:
    """Extrae TITLE y AUTHOR del fichero PDB utilizando un contexto ``with``.

    Los registros PDB pueden continuar en varias líneas, por lo que se concatenan
    todos los fragmentos TITLE y AUTHOR encontrados.
    """
    titulo_partes: list[str] = []
    autor_partes: list[str] = []

    with ruta_pdb.open("r", encoding="utf-8") as fichero:
        for linea in fichero:
            if linea.startswith("TITLE"):
                # En PDB, el contenido del registro comienza después del campo
                # de identificación/continuación. strip() elimina espacios fijos.
                titulo_partes.append(linea[10:].strip())
            elif linea.startswith("AUTHOR"):
                autor_partes.append(linea[10:].strip())

    titulo = " ".join(titulo_partes)
    autores = " ".join(autor_partes)
    return titulo, autores


def extraer_secuencia_p53(ruta_pdb: Path, cadena: str = "A") -> list[str]:
    """Extrae la secuencia SEQRES de p53 para una cadena proteica.

    1TUP contiene tres copias de la proteína p53 (cadenas A, B y C) con la misma
    secuencia. Se utiliza una única cadena (A por defecto) para representar la
    secuencia de p53 una sola vez y evitar triplicar artificialmente el conteo.

    La extracción usa las posiciones fijas del formato PDB:
    - cadena: columna 12 (índice 11 en Python),
    - número total declarado: columnas 14-17,
    - residuos: desde la columna 20.
    """
    secuencia: list[str] = []
    total_declarado: int | None = None

    with ruta_pdb.open("r", encoding="utf-8") as fichero:
        for linea in fichero:
            if not linea.startswith("SEQRES"):
                continue

            cadena_linea = linea[11].strip()
            if cadena_linea != cadena:
                continue

            # En la primera línea SEQRES de la cadena se recoge el total esperado.
            if total_declarado is None:
                total_declarado = int(linea[13:17].strip())

            # Los residuos ocupan el tramo fijo de SEQRES. El filtro excluye ADN
            # y cualquier token que no corresponda a un aminoácido estándar.
            residuos = linea[19:70].split()
            secuencia.extend(
                residuo for residuo in residuos if residuo in AMINOACIDOS_ESTANDAR
            )

    if not secuencia:
        raise ValueError(f"No se encontró una secuencia proteica SEQRES en la cadena {cadena!r}.")

    if total_declarado is not None and len(secuencia) != total_declarado:
        raise ValueError(
            "La secuencia extraída no coincide con el número de residuos declarado "
            f"por el PDB: extraídos={len(secuencia)}, declarados={total_declarado}."
        )

    return secuencia


def comprobar_cadenas_p53_identicas(ruta_pdb: Path) -> bool:
    """Comprueba que las cadenas A, B y C contienen la misma secuencia p53."""
    secuencias = [extraer_secuencia_p53(ruta_pdb, cadena) for cadena in ("A", "B", "C")]
    return secuencias[0] == secuencias[1] == secuencias[2]


def contar_aminoacidos(secuencia: list[str]) -> dict[str, int]:
    """Cuenta aminoácidos sin definir manualmente las claves del diccionario."""
    conteo: dict[str, int] = {}

    for aminoacido in secuencia:
        # La clave se crea dinámicamente al aparecer por primera vez.
        conteo[aminoacido] = conteo.get(aminoacido, 0) + 1

    return conteo


def graficar_frecuencia_aminoacidos(conteo: dict[str, int]) -> Path:
    """Genera y guarda el gráfico de barras solicitado en el apartado 2.D."""
    frecuencias = (
        pd.DataFrame(conteo.items(), columns=["aminoacido", "frecuencia"])
        .sort_values("frecuencia", ascending=False)
        .reset_index(drop=True)
    )

    fig, ax = plt.subplots(figsize=(12, 6))

    # Usar el propio aminoácido como hue permite asignar automáticamente un color
    # distinto a cada barra sin codificar los colores uno a uno.
    sns.barplot(
        data=frecuencias,
        x="aminoacido",
        y="frecuencia",
        hue="aminoacido",
        palette="husl",
        legend=False,
        dodge=False,
        ax=ax,
    )

    ax.set_title("Frecuencia de aminoácidos en la secuencia de p53 (PDB 1TUP, cadena A)")
    ax.set_xlabel("Aminoácido (código de tres letras)")
    ax.set_ylabel("Frecuencia absoluta")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    salida = RESULTADOS_DIR / "frecuencia_aminoacidos_p53.png"
    fig.savefig(salida, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return salida


# -----------------------------------------------------------------------------
# 3. MANIPULACIÓN DE CONJUNTOS DE DATOS (4,5 puntos)
# -----------------------------------------------------------------------------
def analizar_datos_actividad() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resuelve los apartados 3.A-3.F y devuelve datos completos y resumen."""

    # 3.A. Leer actividad.csv. El separador real del fichero proporcionado es ';'.
    datos = pd.read_csv(ACTIVIDAD_FILE, sep=";")

    if datos.shape[1] != 5:
        raise ValueError(
            f"Se esperaban 5 columnas en actividad.csv y se han encontrado {datos.shape[1]}."
        )

    # Renombrado exactamente como solicita el enunciado.
    datos.columns = ["id", "dieta", "pulsaciones", "tiempo", "actividad"]

    print("\n" + "=" * 79)
    print("3.A - DATOS DE ACTIVIDAD Y RENOMBRADO DE COLUMNAS")
    print("=" * 79)
    print(f"Dimensiones iniciales: {datos.shape[0]} filas x {datos.shape[1]} columnas")
    print("Columnas:", list(datos.columns))

    # 3.B. Detectar celdas vacías y eliminar las filas que las contengan.
    nulos_por_columna = datos.isna().sum()
    print("\n3.B - Celdas vacías por columna:")
    print(nulos_por_columna.to_string())

    filas_antes = len(datos)
    datos = datos.dropna().copy()
    filas_despues = len(datos)
    print(f"Filas eliminadas por contener valores vacíos: {filas_antes - filas_despues}")

    # 3.C. value_counts() devuelve en una sola llamada los niveles y su frecuencia.
    frecuencia_dieta = datos["dieta"].value_counts()
    print("\n3.C - Niveles de dieta y frecuencia (value_counts):")
    print(f"Número de niveles: {len(frecuencia_dieta)}")
    print(frecuencia_dieta.to_string())

    # 3.D. Agrupar por actividad y convertir el GroupBy a lista.
    agrupado_actividad = datos.groupby("actividad", sort=False)
    grupos_lista = list(agrupado_actividad)

    print("\n3.D - Agrupación por nivel de actividad:")
    print(f"Número de elementos de la lista: {len(grupos_lista)}")
    print(
        "Hay tres elementos porque la variable 'actividad' tiene tres niveles: "
        + ", ".join(nombre for nombre, _ in grupos_lista)
        + "."
    )
    print(
        "Cada elemento es una tupla (nombre_del_grupo, DataFrame_del_grupo). "
        "Tamaños de los grupos:"
    )
    for nombre, grupo in grupos_lista:
        print(f"  - {nombre!r}: {len(grupo)} filas")

    # 3.E. Estadísticos con agg sobre el objeto GroupBy original.
    resumen_pulsaciones = agrupado_actividad.agg(
        pulsaciones_media=("pulsaciones", "mean"),
        pulsaciones_desviacion_estandar=("pulsaciones", "std"),
    )

    print("\n3.E - Frecuencia cardíaca media y desviación estándar por actividad:")
    print(resumen_pulsaciones.round(3).to_string())

    # 3.F. El segundo fichero es TSV (tabuladores), aunque el texto del enunciado
    # habla de dos CSV. Se lee con sep='\\t' conforme a su formato real.
    ciudades = pd.read_csv(CIUDADES_FILE, sep="\t").rename(columns={"city": "ciudad"})

    datos_completos = datos.merge(
        ciudades,
        on="id",
        how="left",
        validate="many_to_one",  # cada ID debe corresponder a una única ciudad
    )

    if datos_completos["ciudad"].isna().any():
        ids_sin_ciudad = datos_completos.loc[
            datos_completos["ciudad"].isna(), "id"
        ].unique()
        raise ValueError(f"Existen IDs sin ciudad asociada: {ids_sin_ciudad.tolist()}")

    print("\n3.F - Merge con ciudades:")
    print(f"Dimensiones tras el merge: {datos_completos.shape}")
    print("Ciudades incorporadas:", ", ".join(datos_completos["ciudad"].unique()))

    # Se guardan resultados tabulares para facilitar la trazabilidad del análisis.
    datos_completos.to_csv(
        RESULTADOS_DIR / "actividad_con_ciudades.csv", index=False, encoding="utf-8-sig"
    )
    resumen_pulsaciones.to_csv(
        RESULTADOS_DIR / "resumen_pulsaciones_por_actividad.csv",
        encoding="utf-8-sig",
    )

    return datos_completos, resumen_pulsaciones


def graficar_relacion_pulsaciones_tiempo(datos: pd.DataFrame) -> Path:
    """Genera la figura multifacetada solicitada en 3.G.

    Se conserva la columna textual ``tiempo`` y se crea ``tiempo_min`` únicamente
    para que el eje X sea numérico y mantenga correctamente el orden temporal
    1 -> 15 -> 30 minutos.

    Cada faceta representa un estrato dieta x actividad. La línea muestra la
    media de pulsaciones y la banda de error representa ±1 desviación estándar.
    """
    datos_grafico = datos.copy()
    datos_grafico["tiempo_min"] = (
        datos_grafico["tiempo"].str.extract(r"(\d+)", expand=False).astype(int)
    )

    # Orden lógico explícito de los estratos, obtenido de los valores del fichero.
    orden_actividad = ["rest", "walking", "running"]
    orden_dieta = ["low fat", "no fat"]

    g = sns.relplot(
        data=datos_grafico,
        x="tiempo_min",
        y="pulsaciones",
        col="actividad",
        row="dieta",
        col_order=orden_actividad,
        row_order=orden_dieta,
        kind="line",
        estimator="mean",
        errorbar="sd",
        marker="o",
        height=3.2,
        aspect=1.15,
        facet_kws={"margin_titles": True},
    )

    g.set_axis_labels("Tiempo (min)", "Pulsaciones (latidos/min)")
    g.set_titles(row_template="Dieta: {row_name}", col_template="Actividad: {col_name}")
    g.set(xticks=[1, 15, 30])
    g.fig.suptitle(
        "Relación entre pulsaciones y tiempo según actividad y dieta",
        y=1.02,
        fontsize=14,
        fontweight="bold",
    )

    salida = RESULTADOS_DIR / "pulsaciones_tiempo_por_actividad_y_dieta.png"
    g.fig.savefig(salida, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    return salida


def main() -> None:
    # Comprobación temprana de que todos los ficheros necesarios están presentes.
    for ruta in (PDB_FILE, ACTIVIDAD_FILE, CIUDADES_FILE):
        if not ruta.exists():
            raise FileNotFoundError(
                f"No se encuentra {ruta.name}. Sitúalo en la misma carpeta que el script."
            )

    print("=" * 79)
    print("2. MANIPULACIÓN DEL FICHERO PDB 1TUP")
    print("=" * 79)

    # 2.A. TITLE y AUTHOR.
    titulo, autores = leer_metadatos_pdb(PDB_FILE)
    print("2.A - TITLE:", titulo)
    print("2.A - AUTHOR:", autores)

    # 2.B. Secuencia de p53.
    secuencia_p53 = extraer_secuencia_p53(PDB_FILE, cadena="A")
    print("\n2.B - Longitud de la secuencia de p53 (cadena A):", len(secuencia_p53))
    print("Primeros 20 aminoácidos:", secuencia_p53[:20])
    print("Cadenas A, B y C idénticas:", comprobar_cadenas_p53_identicas(PDB_FILE))

    # 2.C. Conteo dinámico en diccionario.
    conteo_aa = contar_aminoacidos(secuencia_p53)
    print("\n2.C - Conteo de aminoácidos:")
    for aminoacido, frecuencia in sorted(
        conteo_aa.items(), key=lambda elemento: (-elemento[1], elemento[0])
    ):
        print(f"  {aminoacido}: {frecuencia}")
    print("Total comprobado:", sum(conteo_aa.values()))

    # 2.D. Gráfico de barras.
    grafico_aa = graficar_frecuencia_aminoacidos(conteo_aa)
    print("\n2.D - Gráfico guardado en:", grafico_aa)

    # 3.A-3.F. Pandas.
    datos_completos, _ = analizar_datos_actividad()

    # 3.G. Figura multifacetada.
    grafico_relacion = graficar_relacion_pulsaciones_tiempo(datos_completos)
    print("\n3.G - Figura multifacetada guardada en:", grafico_relacion)

    print("\n" + "=" * 79)
    print("ACTIVIDAD COMPLETADA CORRECTAMENTE")
    print("=" * 79)


if __name__ == "__main__":
    main()
