# GRADEX-DUAL v1.0: Marco Dual para Estimación del Parámetro GRADEX

**Autor:** Mauricio Javier Victoria Niño  
**ORCID:** [0009-0003-4328-5691](https://orcid.org)  
**Afiliación:** Investigador Independiente  
**Licencia:** MIT

---

## 📝 Descripción
Este repositorio contiene el marco computacional **GRADEX-DUAL v1.0**, desarrollado en lenguaje **R**. El software está diseñado para la estimación del parámetro GRADEX —fundamental en el diseño hidrológico con metodologías como GRADEX-IDF y HEC-HMS— especialmente en cuencas con redes pluviométricas escasas.

El modelo implementa un flujo de trabajo dual:
1.  **Capa Puntual:** Ajuste de distribuciones por estación mediante L-momentos y máxima verosimilitud, con regionalización espacial vía IDW u Ordinary Kriging.
2.  **Capa Regional:** Análisis de frecuencias regional basado en el enfoque de Hosking & Wallis (1997) utilizando L-momentos ponderados.

La innovación principal radica en un **sistema de compuertas de calidad** que integra ambas estimaciones basándose en la incertidumbre (amplitud del intervalo de confianza bootstrap) y la capacidad predictiva espacial ($R^2$ de validación cruzada).

## 📁 Estructura del Repositorio
Para garantizar la ejecución directa del script, los archivos se encuentran en el directorio raíz:
*   `gradex_dual_v1.0.R`: Script principal estructurado en 19 secciones que contiene toda la lógica del modelo.
*   `datos_precipitacion.xlsx`: Archivo Excel con las series de precipitación máxima anual y coordenadas del caso de estudio (Alto Río Cauca, Colombia).
*   `LICENCIA`: Términos de uso bajo la licencia MIT.
*   `.gitignore`: Archivo de configuración para excluir archivos temporales de R.

## 🚀 Requisitos e Instalación
Para ejecutar el script, asegúrese de tener instaladas las siguientes librerías en R:
```R
install.packages(c("readxl", "writexl", "Lmoments", "lmomRFA", "ggplot2", 
                   "nortest", "sp", "gstat", "automap", "sf", "tseries"))
```

## 📖 Cita Recomendada
Si utiliza este código o el software para su investigación, por favor cite el preprint asociado en EarthArXiv:

> **Victoria-Niño, M. J. (2026).** *GRADEX-DUAL: Un marco dual que combina el análisis frecuencial puntual y regional para la estimación del parámetro GRADEX en cuencas con redes pluviométricas escasas.* EarthArXiv (Preprint). DOI: [El DOI será asignado tras la carga en EarthArXiv]

## ⚖️ Licencia
Este proyecto está bajo la **Licencia MIT**. Esto permite el uso, copia y modificación libre del código, siempre que se mantenga el crédito al autor original.

