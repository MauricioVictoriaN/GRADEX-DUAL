# GRADEX-DUAL v1.0: Marco Dual para Estimación del Parámetro GRADEX

**Autor:** Mauricio Javier Victoria Niño  
**ORCID:** [0009-0003-4328-5691](https://orcid.org/0009-0003-4328-5691)  
**Afiliación:** Investigador Independiente  
**Licencia:** MIT  

[![DOI](https://img.shields.io/badge/DOI-10.31224%2F6945-blue)](https://doi.org/10.31224/6945)
[![Licencia: MIT](https://img.shields.io/badge/Licencia-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lenguaje: R](https://img.shields.io/badge/Lenguaje-R-276DC3.svg)](https://www.r-project.org/)

---

## 📝 Descripción

Este repositorio contiene el marco computacional **GRADEX-DUAL v1.0**, desarrollado en lenguaje **R**. El software está diseñado para la **estimación del parámetro GRADEX** —fundamental en diseño hidrológico con metodologías como GRADEX-IDF y HEC-HMS— especialmente en cuencas con **redes pluviométricas escasas**.

El modelo implementa un **flujo de trabajo dual**:

- **Capa Puntual:** Ajuste de distribuciones por estación mediante **L‑momentos** y **máxima verosimilitud**, con regionalización espacial vía **IDW** u **Ordinary Kriging**.
- **Capa Regional:** Análisis de frecuencias regional basado en el enfoque de **Hosking & Wallis (1997)** utilizando L‑momentos ponderados.

La **innovación principal** radica en un **sistema de compuertas de calidad** que integra ambas estimaciones basándose en la incertidumbre (amplitud del intervalo de confianza bootstrap) y la capacidad predictiva espacial (R² de validación cruzada).

---

## 📁 Estructura del Repositorio

| Archivo | Descripción |
|---------|-------------|
| `gradex_dual_v1.0.0.R` | Script principal estructurado en 19 secciones. Contiene toda la lógica del modelo. |
| `datos_precipitacion.xlsx` | Series de precipitación máxima anual y coordenadas (caso de estudio: Alto Río Cauca, Colombia). |
| `LICENCIA` | Términos de uso bajo licencia MIT. |
| `.gitignore` | Configuración para excluir archivos temporales de R. |

---

## 🚀 Requisitos e Instalación

```r
install.packages(c("readxl", "writexl", "Lmoments", "lmomRFA", "ggplot2",
                   "nortest", "sp", "gstat", "automap", "sf", "tseries"))
```

---

## 📖 Cita Recomendada

Victoria-Niño, M. J. (2026). GRADEX-DUAL: Un marco dual que combina el análisis frecuencial puntual y regional para la estimación del parámetro GRADEX en cuencas con redes pluviométricas escasas. *EngrXiv*. DOI: [10.31224/6945](https://doi.org/10.31224/6945)

```bibtex
@unpublished{Victoria2026GRADEX,
  author = {Victoria-Niño, Mauricio Javier},
  title  = {GRADEX-DUAL: Un marco dual para la estimación del parámetro GRADEX
            en cuencas con redes pluviométricas escasas},
  year   = {2026},
  doi    = {10.31224/6945},
  note   = {Preprint en EngrXiv}
}
```

---

## ⚖️ Licencia

Este proyecto está bajo la **Licencia MIT**. Se permite el uso, copia y modificación libre del código, siempre que se mantenga el crédito al autor original.

---

## ✉️ Contacto

Para dudas, sugerencias o reportar problemas, abra un *issue* en este repositorio o contacte al autor vía ORCID: [0009-0003-4328-5691](https://orcid.org/0009-0003-4328-5691).
