from datetime import datetime
import json
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
    """Clase principal para calcular las probabilidades de grupo sanguíneo en la descendencia y gestionar persistencia JSON."""

    # Descriptores para la validación de los datos
    padre_grupo = GruposangDescriptor("padre_grupo")
    madre_grupo = GruposangDescriptor("madre_grupo")

    # Atributo de clase: mapeo de fenotipo a posibles genotipos
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
        """Método privado para obtener las combinaciones de alelos de un fenotipo."""
        genotipos = self.GENOTIPOS_POSIBLES[grupo]
        return [tuple(g) for g in genotipos]

    def _determinar_fenotipo(self, alelo1: str, alelo2: str) -> str:
        """Método privado que determina el grupo sanguíneo según dominancia mendeliana."""
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
        """Calcula el porcentaje de probabilidad para cada fenotipo en los descendientes."""
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
        """Muestra por consola los porcentajes calculados."""
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
        """Muestra un gráfico circular (pie chart) con los resultados."""
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
    # MÉTODOS AÑADIDOS EN LA FASE 2
    # ==========================================

    def guardar_analisis_json(self, ruta_archivo: str = "analisis_resultados.json"):
        """Guarda el análisis actual en un archivo JSON externo agregando fecha de análisis,

        grupos de los padres y los resultados de la descendencia.
        """
        if not self._porcentajes:
            self.calcular_probabilidades()

        # Estructura del registro a guardar
        nuevo_registro = {
            "fecha_analisis": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "padres": {"padre": self.padre_grupo, "madre": self.madre_grupo},
            "resultados_descendencia": self._porcentajes,
        }

        # Cargar registros existentes si el archivo ya existe
        registros = []
        try:
            with open(ruta_archivo, "r", encoding="utf-8") as archivo:
                registros = json.load(archivo)
                if not isinstance(registros, list):
                    registros = [registros]
        except (FileNotFoundError, json.JSONDecodeError):
            registros = []

        registros.append(nuevo_registro)

        # Guardar en archivo JSON
        with open(ruta_archivo, "w", encoding="utf-8") as archivo:
            json.dump(registros, archivo, indent=4, ensure_ascii=False)

        print(f"✅ Análisis guardado correctamente en '{ruta_archivo}'.")

    @classmethod
    def procesar_desde_json(cls, ruta_archivo_entrada: str, ruta_archivo_salida: str = "analisis_resultados.json"):
        """Lee múltiples parejas desde un archivo JSON externo, realiza el análisis

        para cada una y guarda los resultados agregados.
        """
        try:
            with open(ruta_archivo_entrada, "r", encoding="utf-8") as archivo:
                datos_padres = json.load(archivo)
        except FileNotFoundError:
            print(f"❌ Error: El archivo '{ruta_archivo_entrada}' no existe.")
            return
        except json.JSONDecodeError:
            print(
                f"❌ Error: El archivo '{ruta_archivo_entrada}' no tiene un formato JSON válido."
            )
            return

        print(
            f"\n--- Procesando registros desde '{ruta_archivo_entrada}' ---"
        )
        for i, registro in enumerate(datos_padres, 1):
            padre = registro.get("padre")
            madre = registro.get("madre")

            try:
                instancia = cls(padre_grupo=padre, madre_grupo=madre)
                instancia.calcular_probabilidades()
                instancia.mostrar_resultados()
                instancia.guardar_analisis_json(
                    ruta_archivo=ruta_archivo_salida
                )
            except (ValueError, TypeError) as e:
                print(
                    f"⚠️ Error procesando el registro {i} (Padre: '{padre}', Madre: '{madre}'): {e}"
                )


# ==========================================
# EJEMPLO DE USO - FASE 2
# ==========================================
if __name__ == "__main__":
    # 1. Ejemplo de guardar un único análisis en JSON:
    analisis = CalculadoraHerenciaSangre(padre_grupo="A", madre_grupo="B")
    analisis.mostrar_resultados()
    analisis.guardar_analisis_json("resultados_fase2.json")

    # 2. Ejemplo de lectura masiva de parejas desde un JSON de entrada:
    # Creamos un archivo de entrada de prueba 'entrada_padres.json'
    datos_ejemplo = [
        {"padre": "A", "madre": "O"},
        {"padre": "AB", "madre": "B"},
        {"padre": "O", "madre": "O"},
    ]
    with open("entrada_padres.json", "w", encoding="utf-8") as f:
        json.dump(datos_ejemplo, f, indent=4)

    # Leemos del archivo de entrada y procesamos automáticamente guardando la fecha y resultados
    CalculadoraHerenciaSangre.procesar_desde_json(
        ruta_archivo_entrada="entrada_padres.json",
        ruta_archivo_salida="resultados_fase2.json",
    )