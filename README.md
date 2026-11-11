# Estimación de la degradación del compresor de una turbina a gas naval

Proyecto Integrador: Regresión · Caso 1 · **Equipo 5** · Estadística Aplicada (USACH, 2026).

Equipo 5: Claudio Andrés Pérez Lagos, Daniela Sofía Sepúlveda Serrano, Ariel Felipe Garrido Armijo,
Benjamín Fernando Sanchez Vargas y Elías Andrés Rita Rodríguez.

## Contenido del repositorio

| Carpeta / archivo | Qué contiene |
|---|---|
| `codigo/00_main.R` | Punto de entrada único: ejecuta todo el análisis |
| `codigo/R/` | Configuración (`00_config.R`), utilidades y módulos |
| `codigo/pasos/` | Pasos 10 a 26 del análisis |
| `codigo/referencia/` | Referencia para verificar la reproducibilidad |
| `datos/data.txt` | Dataset (UCI, CC BY 4.0) |
| `resultados_*` | Tablas, figuras y registros de la ejecución oficial |
| `entrega/` | Informe final (PDF) y presentación (PPTX) |
| `equipo_regresion_naval.Rproj` | Proyecto de RStudio (fija el directorio de trabajo) |
| `CAMBIOS.md` | Registro de cambios del código |

Para obtenerlo: botón **Code → Download ZIP** (y descomprimir), o

```bash
git clone https://github.com/claudioperezl-usach/Grupo_5_Caso_1_de_Estudios_Estadistica_Aplicada.git
```

## 1. Datos

* Fuente: *Condition Based Maintenance of Naval Propulsion Plants*, Coraddu et al.,
  UCI Machine Learning Repository, DOI [10.24432/C5K31K](https://doi.org/10.24432/C5K31K),
  licencia CC BY 4.0.
* Archivo: `datos/data.txt` (incluido). Si falta, descargue el zip desde la página de UCI
  y copie `UCI CBM Dataset/data.txt` en `datos/`.
* Control: `codigo/R/02_datos.R` verifica que haya 11.934 filas y que la suma de todas las
  celdas coincida con `cfg$firma_datos`. Si usa otra versión del dataset, la ejecución se
  detiene con un mensaje.

## 2. Instalación (una vez)

R ≥ 4.3. En R:

```r
install.packages(c("data.table", "dplyr", "tibble", "tidyr", "ggplot2", "rlang",
                   "rsample", "recipes", "parsnip", "workflows", "tune", "yardstick",
                   "dials", "generics", "glmnet", "ranger", "xgboost", "rpart",
                   "sandwich", "mgcv", "gridExtra"))
```

Las versiones usadas quedan registradas en
`resultados_elastic_net_final/versiones_paquetes.csv` y `sessionInfo.txt`.

## 3. Ejecución (un solo comando)

Abra `equipo_regresion_naval.Rproj` (en la raíz del repositorio) en RStudio, reinicie R
(Ctrl+Shift+F10) y ejecute:

```r
source("codigo/00_main.R")
```

Desde la terminal, en la carpeta raíz:

```bash
Rscript codigo/00_main.R
```

* Ejecuta los pasos 10 a 26 en orden, borra los resultados anteriores y regenera todas
  las tablas y figuras. Tiempo: entre 60 y 100 minutos según el computador con `anidar_ensembles = TRUE`
  (60,8 min en la ejecución oficial; la validación
  anidada de Random Forest y XGBoost es la parte lenta); unos 5 minutos con `FALSE`.
* Semilla global `20261016` y un solo hilo, en `codigo/R/00_config.R`. Es el único lugar
  donde se cambian parámetros.
* `anidar_ensembles = TRUE`: árbol, Random Forest y XGBoost eligen sus hiperparámetros con
  LOSO anidado y la regla 1-SE, como los demás modelos. Con `FALSE` se usan configuraciones
  fijas (`arbol_fijo`, `xgb_fijo`, `rf_arboles`) y la ejecución es más rápida.
* Reproducibilidad: el paso 11 compara el ranking anidado con
  `codigo/referencia/ranking_anidado_referencia.csv` (tolerancia 1e-3). Si cambia la
  configuración, se crea una referencia nueva y hay que ejecutar dos veces.
* La referencia y los resultados incluidos son de la ejecución oficial del equipo (Windows 11,
  R 4.6.1; versiones en `resultados_elastic_net_final/versiones_paquetes.csv`). Con versiones
  antiguas de xgboost (1.7) solo cambia XGBoost en el tercer decimal (RMSE 1,776 frente a 1,778),
  y la comparación del paso 11 lo señalará.

### Advertencias esperadas

glmnet avisa que no converge en el extremo de su trayectoria interna. Solo ocurre con Lasso
y Elastic-net y en λ ≤ 4·10⁻⁴: la trayectoria se trunca allí. Esos λ no se eligen: el menor
λ elegido por la regla 1-SE en cualquier partición es 0,0022, y el final es 0,0060. Por eso,
en los dos λ más pequeños de la grilla (1e-4 y 2,8e-4) la predicción usa el último λ
convergido; esto solo afecta al paso 23 (CV aleatoria), donde se reporta como "penalty = 1e-04".
El paso 17 reajusta el modelo final con tolerancia 1e-10 y registra la diferencia en
`17_convergencia_estricta.csv`. No hay otras advertencias.

## 4. Dónde está cada cifra del informe

| Informe | Archivo (carpeta `tablas/`) | Paso |
|---|---|---|
| §1 Auditoría, grilla y partición | `01_*`, `02_particion.csv` | 10 |
| §2 Relación T2–pérdida con IC (HC3 y clásico) | `04_relacion_ic.csv` | 10 |
| §3 CV aleatoria vs LOSO (solo entrenamiento) | `23_cv_aleatoria_vs_loso.csv` | 23 |
| §4 Hiperparámetros finales (λ, α) | `12_parametros_elegidos.csv` | 12 |
| §4 OLS: VIF, residuos y Cook | `05_ols_*`, `18_*` | 18 |
| §5 Tabla 1 y Tabla 2 (validación LOSO) | `21_tabla1_presentacion.csv`, `21_tabla1_validacion.csv`, `21_comparacion_en_ridge.csv` | 21 |
| §5 Figura 1 (RMSE por velocidad) | `figuras/13_rmse_por_velocidad.png` | 13 |
| §5 Test, IC bootstrap y subgrupo | `12_resultados_test.csv`, `12_test_incertidumbre_subgrupo.csv` | 12 |
| §5 Intervalos de predicción | `16_resumen_intervalos.csv` | 16 |
| §5 Sensibilidad a ruido (LOSO, solo entrenamiento) | `15_resumen_sensibilidad_ruido.csv`, `15_sensibilidad_ruido_por_velocidad.csv` | 15 |
| §6 Coeficientes y estabilidad | `14_coeficientes_elastic_net.csv`, `17_estabilidad_coeficientes.csv` | 14, 17 |
| §6 Importancia por permutación, PDP e ICE | `22_*.csv`, `figuras/22_*.png` | 22 |
| §6 Umbral de inspección | `16_politica_alerta.csv` | 16 |
| §6 Multicolinealidad (correlaciones y VIF) | `24_correlaciones_predictores.csv`, `figuras/24_correlaciones_predictores.png`, `figuras/24_vif.png` | 24 |
| §6 Valores SHAP del Elastic-net | `24_shap_importancia.csv`, `figuras/24_shap_enjambre.png` | 24 |
| Presentación: figuras de las láminas 4, 6, 7, 8 y 10 | `figuras/presentacion/26_*.png` | 26 |
| §5 Error de entrenamiento vs validación (curva según λ) | `25_curva_entrenamiento_validacion.csv`, `figuras/25_curva_entrenamiento_validacion.png` | 25 |

Las tablas de los pasos 11 y 23 están en `resultados_validacion_anidada/`; las demás,
en `resultados_elastic_net_final/`.
`REGISTRO_EJECUCION.txt` resume tiempos, semilla y el modelo ganador.

## 5. Uso de los test

Los test se usan en un único punto del pipeline: el paso 12 evalúa una vez el modelo
congelado. Los pasos 13, 14 y 16 solo leen esas predicciones para graficarlas y para
verificar la cobertura y la política de alerta. Ninguna cifra calculada con los test se usa
para elegir partición, modelo, hiperparámetros, intervalos ni umbral. Todos los demás
análisis (validación, sensibilidad a ruido, importancia, PDP/ICE) usan solo el entrenamiento. Los test se examinaron antes en un análisis
exploratorio previo (no incluido en este repositorio, porque ya no se ejecuta), por lo que la
evaluación en test es confirmatoria y no totalmente independiente.
