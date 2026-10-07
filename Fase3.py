from datetime import datetime
import json
import os
import shutil
import matplotlib.pyplot as plt


class GruposangDescriptor:
    """Descriptor para validar que el grupo sanguíneo introducido sea válido."""

    GRUPOS_VALIDOS = {"A", "B", "AB", "O"}

    def __init__(self, nombre_atributo):
        self.nombre_atributo = f"_{nombre_atributo}"

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return getattr(instance, self.nombre_atributo, None)

    def __set__(self, instance, value):
        if not isinstance(value, str):
            raise TypeError("El grupo sanguíneo debe ser una cadena de texto.")

        valor_normalizado = value.strip().upper()
        if valor_normalizado not in self.GRUPOS_VALIDOS:
            raise ValueError(
                f"Grupo sanguíneo '{value}' no válido. Opciones permitidas: {', '.join(sorted(self.GRUPOS_VALIDOS))}"
            )

        setattr(instance, self.nombre_atributo, valor_normalizado)


class CalculadoraHerenciaSangre:
    """Clase principal para calcular las probabilidades de grupo sanguíneo en la descendencia,

    gestionar persistencia JSON y procesar flujos de trabajo en carpetas.
    """

    # Descriptores para la validación de los datos
    padre_grupo = GruposangDescriptor("padre_grupo")
    madre_grupo = GruposangDescriptor("madre_grupo")

    # Atributo de clase
    GENOTIPOS_POSIBLES = {
        "A": ["AA", "AO"],
        "B": ["BB", "BO"],
        "AB": ["AB"],
        "O": ["OO"],
    }

    def __init__(self, padre_grupo: str, madre_grupo: str):
        self.padre_grupo = padre_grupo
        self.madre_grupo = madre_grupo
        self._porcentajes = {}

    def _obtener_alelos_posibles(self, grupo: str) -> list[tuple[str, str]]:
        genotipos = self.GENOTIPOS_POSIBLES[grupo]
        return [tuple(g) for g in genotipos]

    def _determinar_fenotipo(self, alelo1: str, alelo2: str) -> str:
        alelos = sorted([alelo1, alelo2])
        genotipo = "".join(alelos)

        if genotipo in ("AA", "AO"):
            return "A"
        elif genotipo in ("BB", "BO"):
            return "B"
        elif genotipo == "AB":
            return "AB"
        elif genotipo == "OO":
            return "O"

    def calcular_probabilidades(self) -> dict[str, float]:
        genotipos_padre = self._obtener_alelos_posibles(self.padre_grupo)
        genotipos_madre = self._obtener_alelos_posibles(self.madre_grupo)

        total_combinaciones = 0
        conteo_fenotipos = {"A": 0, "B": 0, "AB": 0, "O": 0}

        for g_padre in genotipos_padre:
            for g_madre in genotipos_madre:
                for alelo_p in g_padre:
                    for alelo_m in g_madre:
                        fenotipo = self._determinar_fenotipo(alelo_p, alelo_m)
                        conteo_fenotipos[fenotipo] += 1
                        total_combinaciones += 1

        self._porcentajes = {
            grupo: round((conteo / total_combinaciones) * 100, 2)
            for grupo, conteo in conteo_fenotipos.items()
            if conteo > 0
        }

        return self._porcentajes

    def mostrar_resultados(self):
        if not self._porcentajes:
            self.calcular_probabilidades()

        print(
            f"\nResultados para el cruce: Padre ({self.padre_grupo}) x Madre ({self.madre_grupo})"
        )
        print("-" * 55)
        for grupo, porcentaje in self._porcentajes.items():
            print(f"  • Grupo Sanguíneo {grupo}: {porcentaje}%")
        print("-" * 55)

    def graficar_resultados(self):
        if not self._porcentajes:
            self.calcular_probabilidades()

        etiquetas = [
            f"Grupo {grupo}" for grupo in self._porcentajes.keys()
        ]
        valores = list(self._porcentajes.values())
        colores = ["#ff9999", "#66b3ff", "#99ff99", "#ffcc99"]

        plt.figure(figsize=(7, 7))
        plt.pie(
            valores,
            labels=etiquetas,
            autopct="%1.1f%%",
            startangle=140,
            colors=colores[: len(valores)],
            explode=[0.05] * len(valores),
        )
        plt.title(
            f"Probabilidad de Grupo Sanguíneo en la Descendencia\nPadres: {self.padre_grupo} x {self.madre_grupo}"
        )
        plt.show()

    # ==========================================
    # MÉTODOS DE LA FASE 3: FLUJO DE ARCHIVOS
    # ==========================================

    @classmethod
    def procesar_flujo_archivos(cls, base_dir: str = "data"):
        """Procesa todos los archivos JSON en data/pending:

        1. Genera los resultados del análisis en data/results.
        2. Copia el archivo procesado a data/done.
        3. Elimina el archivo original de data/pending.
        """
        # Definición de rutas
        pending_dir = os.path.join(base_dir, "pending")
        done_dir = os.path.join(base_dir, "done")
        results_dir = os.path.join(base_dir, "results")

        # Asegurar que las tres carpetas existan
        for carpeta in [pending_dir, done_dir, results_dir]:
            os.makedirs(carpeta, exist_ok=True)

        # Listar archivos JSON pendientes
        archivos_pendientes = [
            f for f in os.listdir(pending_dir) if f.endswith(".json")
        ]

        if not archivos_pendientes:
            print(f"📂 No hay archivos pendientes de procesar en '{pending_dir}'.")
            return

        print(f"\n🔄 Iniciando el procesamiento de {len(archivos_pendientes)} archivo(s)...")

        for nombre_archivo in archivos_pendientes:
            ruta_pending = os.path.join(pending_dir, nombre_archivo)
            ruta_done = os.path.join(done_dir, nombre_archivo)
            ruta_results = os.path.join(
                results_dir, f"resultado_{nombre_archivo}"
            )

            try:
                # 1. Leer registros del archivo en pending
                with open(ruta_pending, "r", encoding="utf-8") as f:
                    datos = json.load(f)

                if isinstance(datos, dict):
                    datos = [datos]

                resultados_totales = []

                # 2. Procesar cada pareja dentro del archivo
                for elemento in datos:
                    padre = elemento.get("padre")
                    madre = elemento.get("madre")

                    calc = cls(padre_grupo=padre, madre_grupo=madre)
                    porcentajes = calc.calcular_probabilidades()

                    resultados_totales.append(
                        {
                            "fecha_analisis": datetime.now().strftime(
                                "%Y-%m-%d %H:%M:%S"
                            ),
                            "padres": {"padre": padre, "madre": madre},
                            "resultados_descendencia": porcentajes,
                        }
                    )

                # 3. Guardar resultados en data/results
                with open(ruta_results, "w", encoding="utf-8") as f:
                    json.dump(resultados_totales, f, indent=4, ensure_ascii=False)

                # 4. Copiar el archivo original a data/done
                shutil.copy2(ruta_pending, ruta_done)

                # 5. Eliminar el archivo procesado de data/pending
                os.remove(ruta_pending)

                print(
                    f"✅ Procesado con éxito: '{nombre_archivo}' -> Resultados en 'results/', movido a 'done/'"
                )

            except Exception as e:
                print(
                    f"❌ Error al procesar el archivo '{nombre_archivo}': {e}"
                )


# ==========================================
# PRUEBA Y DEMOSTRACIÓN DE LA FASE 3
# ==========================================
if __name__ == "__main__":
    # 1. Crear un lote de prueba en data/pending
    os.makedirs("data/pending", exist_ok=True)

    muestra_1 = [
        {"padre": "A", "madre": "B"},
        {"padre": "O", "madre": "AB"},
    ]
    muestra_2 = [
        {"padre": "A", "madre": "A"},
        {"padre": "B", "madre": "O"},
    ]

    with open("data/pending/lote_01.json", "w", encoding="utf-8") as f:
        json.dump(muestra_1, f, indent=4)

    with open("data/pending/lote_02.json", "w", encoding="utf-8") as f:
        json.dump(muestra_2, f, indent=4)

    # 2. Ejecutar la función de procesamiento de flujos
    CalculadoraHerenciaSangre.procesar_flujo_archivos(base_dir="data")