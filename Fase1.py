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
    """Clase principal para calcular las probabilidades de grupo sanguíneo en la descendencia."""

    # Descriptores para la validación de las entradas de los padres
    padre_grupo = GruposangDescriptor("padre_grupo")
    madre_grupo = GruposangDescriptor("madre_grupo")

    # Mapeo de fenotipo a posibles genotipos
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
        """Método privado que determina el grupo sanguíneo final según dominancia mendeliana."""
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

        # Evaluamos cada combinación genotípica posible entre ambos padres
        for g_padre in genotipos_padre:
            for g_madre in genotipos_madre:
                # Cuadro de Punnett para cada combinación de genotipos
                for alelo_p in g_padre:
                    for alelo_m in g_madre:
                        fenotipo = self._determinar_fenotipo(alelo_p, alelo_m)
                        conteo_fenotipos[fenotipo] += 1
                        total_combinaciones += 1

        # Generar porcentajes finales
        self._porcentajes = {
            grupo: (conteo / total_combinaciones) * 100
            for grupo, conteo in conteo_fenotipos.items()
            if conteo > 0
        }

        return self._porcentajes

        # Generar porcentajes finales
        self._porcentajes = {
            grupo: (conteo / total_combinaciones) * 100
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
            print(f"  • Grupo Sanguíneo {grupo}: {porcentaje:.2f}%")
        print("-" * 55)

    def graficar_resultados(self):
        """Muestra un gráfico circular (pie chart) con los resultados."""
        if not self._porcentajes:
            self.calcular_probabilidades()

        etiquetas = [
            f"Grupo {grupo}" for grupo in self._porcentajes.keys()
        ]
        valores = list(self._porcentajes.values())
        colores = ["#e59177", "#e091d7", "#99ff99", "#E7CF64"]

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
# EJEMPLO DE USO
# ==========================================
if __name__ == "__main__":
    # Ejemplo 1: Padre grupo A y Madre grupo B
    cruce = CalculadoraHerenciaSangre(padre_grupo="A", madre_grupo="B")
    cruce.mostrar_resultados()
    cruce.graficar_resultados()