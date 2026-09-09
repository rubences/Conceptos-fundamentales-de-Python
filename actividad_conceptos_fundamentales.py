"""
Script de apoyo para la actividad "Conceptos fundamentales de Python".

Incluye:
1) Comandos de Conda para preparar el entorno.
2) Procesado de archivo PDB (TITLE, AUTHOR, SEQRES, conteo y gráfico).
3) Procesado de dos CSV con pandas (limpieza, groupby, agg, merge y gráfico).
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def comandos_conda() -> list[str]:
    """Devuelve comandos para los apartados 1A, 1B y 1C."""
    return [
        "conda create -n actividad1 python=3.12 pandas matplotlib seaborn spyder -y",
        "conda activate actividad1",
        "conda install -n actividad1 -c conda-forge spyder -y",
        "conda env list",
        "conda list -n actividad1",
    ]


def extraer_title_y_author_pdb(ruta_pdb: str | Path) -> tuple[str, str]:
    """Extrae TITLE y AUTHOR del archivo PDB."""
    titulo_partes: list[str] = []
    autores_partes: list[str] = []

    with open(ruta_pdb, "r", encoding="utf-8") as archivo:
        for linea in archivo:
            if linea.startswith("TITLE"):
                titulo_partes.append(linea[10:].strip())
            elif linea.startswith("AUTHOR"):
                autores_partes.append(linea[10:].strip())

    return " ".join(titulo_partes), " ".join(autores_partes)


def extraer_aa_seqres(ruta_pdb: str | Path) -> list[str]:
    """
    Lee el PDB con `with` y extrae únicamente aminoácidos de líneas SEQRES.
    Devuelve una lista de códigos de 3 letras.
    """
    aa_validos = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
    secuencia: list[str] = []

    with open(ruta_pdb, "r", encoding="utf-8") as archivo:
        for linea in archivo:
            if not linea.startswith("SEQRES"):
                continue
            tokens = linea.strip().split()
            for token in tokens[4:]:
                if token in aa_validos:
                    secuencia.append(token)

    return secuencia


def contar_aminoacidos(secuencia: list[str]) -> dict[str, int]:
    """Cuenta frecuencia de aminoácidos sin definir claves manualmente."""
    return dict(Counter(secuencia))


def graficar_frecuencia_aminoacidos(frecuencias: dict[str, int], ruta_salida: str | Path) -> None:
    """Genera gráfico de barras de frecuencia de aminoácidos con Seaborn."""
    if not frecuencias:
        raise ValueError("No hay datos de aminoácidos para graficar.")

    datos = pd.DataFrame(
        sorted(frecuencias.items(), key=lambda x: x[0]),
        columns=["aminoacido", "frecuencia"],
    )

    plt.figure(figsize=(12, 5))
    palette = sns.color_palette("husl", n_colors=len(datos))
    sns.barplot(data=datos, x="aminoacido", y="frecuencia", hue="aminoacido", palette=palette, legend=False)
    plt.title("Frecuencia de aminoácidos en la proteína P53 (1TUP)")
    plt.xlabel("Aminoácido")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plt.savefig(ruta_salida, dpi=150)
    plt.close()


def cargar_y_preparar_datos(ruta_csv: str | Path) -> pd.DataFrame:
    """Lee el CSV principal, renombra columnas y elimina filas con vacíos."""
    df = pd.read_csv(ruta_csv)

    # Si pandas usó la primera fila de datos como cabecera, se relee sin cabecera.
    columnas_originales = [str(c) for c in df.columns]
    tokens_numericos = sum(token.replace(".", "", 1).isdigit() for token in columnas_originales)
    if tokens_numericos >= 2:
        df = pd.read_csv(ruta_csv, header=None)

    df.columns = ["id", "dieta", "pulsaciones", "tiempo", "actividad"]
    df = df.dropna()
    return df


def niveles_y_frecuencia_dieta(df: pd.DataFrame) -> pd.Series:
    """Devuelve niveles y frecuencias de dieta con un único método."""
    return df["dieta"].value_counts()


def agrupar_por_actividad(df: pd.DataFrame) -> tuple[pd.core.groupby.DataFrameGroupBy, list[tuple[str, pd.DataFrame]]]:
    """
    Genera groupby por actividad y también lista de agrupaciones.
    La lista tendrá tres elementos si existen tres niveles en 'actividad'.
    """
    agrupado = df.groupby("actividad")
    lista_grupos = list(agrupado)
    return agrupado, lista_grupos


def resumen_pulsaciones(agrupado: pd.core.groupby.DataFrameGroupBy) -> pd.DataFrame:
    """Calcula media y desviación estándar de pulsaciones por actividad."""
    return agrupado.agg(
        pulsaciones_media=("pulsaciones", "mean"),
        pulsaciones_std=("pulsaciones", "std"),
    )


def unir_con_ciudades(df_principal: pd.DataFrame, ruta_csv_ciudades: str | Path) -> pd.DataFrame:
    """Completa dataset haciendo merge con archivo de ciudades por id."""
    df_ciudades = pd.read_csv(ruta_csv_ciudades)

    # Normalización mínima para soportar encabezados frecuentes.
    columnas = [c.strip().lower() for c in df_ciudades.columns]
    df_ciudades.columns = columnas
    if "id" not in df_ciudades.columns:
        df_ciudades = df_ciudades.rename(columns={df_ciudades.columns[0]: "id"})
    if "ciudad" not in df_ciudades.columns and len(df_ciudades.columns) > 1:
        df_ciudades = df_ciudades.rename(columns={df_ciudades.columns[1]: "ciudad"})

    return df_principal.merge(df_ciudades, on="id", how="left")


def graficar_pulsaciones_tiempo_por_actividad_dieta(df: pd.DataFrame, ruta_salida: str | Path) -> None:
    """Crea figura multi-facetada de pulsaciones vs tiempo por actividad y dieta."""
    g = sns.relplot(
        data=df,
        x="tiempo",
        y="pulsaciones",
        hue="dieta",
        col="actividad",
        kind="scatter",
        height=4,
        aspect=1,
    )
    g.figure.suptitle("Relación entre pulsaciones y tiempo por actividad y dieta", y=1.05)
    g.set_axis_labels("Tiempo", "Pulsaciones")
    g.savefig(ruta_salida, dpi=150)
    plt.close(g.figure)


def main() -> None:
    """
    Ejemplo de uso.
    Ajusta los nombres de archivo si son distintos.
    """
    print("=== Comandos de Conda (actividad 1) ===")
    for comando in comandos_conda():
        print(comando)

    ruta_pdb = Path("1TUP.pdb")
    ruta_csv_principal = Path("nombre_archivo.csv")
    ruta_csv_ciudades = Path("ciudades.csv")

    if ruta_pdb.exists():
        titulo, autores = extraer_title_y_author_pdb(ruta_pdb)
        print("\n=== Apartado 2A ===")
        print(f"TITLE: {titulo}")
        print(f"AUTHOR: {autores}")

        secuencia = extraer_aa_seqres(ruta_pdb)
        frecuencias = contar_aminoacidos(secuencia)
        graficar_frecuencia_aminoacidos(frecuencias, "frecuencia_aminoacidos.png")
        print(f"Total de aminoácidos SEQRES extraídos: {len(secuencia)}")
    else:
        print(f"\nNo se encontró el archivo PDB: {ruta_pdb}")

    if ruta_csv_principal.exists() and ruta_csv_ciudades.exists():
        df = cargar_y_preparar_datos(ruta_csv_principal)
        print("\n=== Apartado 3C: niveles y frecuencia de dieta ===")
        print(niveles_y_frecuencia_dieta(df))

        agrupado, lista_grupos = agrupar_por_actividad(df)
        print("\n=== Apartado 3D: lista de grupos ===")
        print(f"Número de elementos en la lista: {len(lista_grupos)}")
        for actividad, bloque in lista_grupos:
            print(f"- Grupo '{actividad}' con {len(bloque)} filas")

        print("\n=== Apartado 3E: media y desviación estándar de pulsaciones ===")
        print(resumen_pulsaciones(agrupado))

        df_completo = unir_con_ciudades(df, ruta_csv_ciudades)
        print("\n=== Apartado 3F: resultado merge (primeras filas) ===")
        print(df_completo.head())

        graficar_pulsaciones_tiempo_por_actividad_dieta(df_completo, "pulsaciones_tiempo_actividad_dieta.png")
    else:
        print(
            f"\nNo se encontraron ambos CSV requeridos: {ruta_csv_principal} y {ruta_csv_ciudades}. "
            "Ajusta los nombres en main() si fuera necesario."
        )


if __name__ == "__main__":
    main()
