# Registro de cambios

**Estado actual:** el proyecto se ejecutó completo desde una sesión limpia y la segunda
ejecución reproduce el ranking anidado (ver "Revisión 3" al final). Pasos 10 a 26.
Estructura: `codigo/00_main.R`, `codigo/R/` (configuración, utilidades y módulos),
`codigo/pasos/` (10 a 26) e `historico/modulos_antiguos/` (flujo antiguo, ya no se ejecuta).

---

# Revisión 1 (histórica: escrita sin ejecutar R; sus pendientes se cerraron en las revisiones 2 y 3)

## Nuevo
- `pasos/22_interpretacion_pdp_ice.R`: importancia por permutación dentro de la LOSO (sin test), PDP e ICE del Elastic-net final sobre filas de entrenamiento. Agregado a la lista de pasos del main. Salidas: `22_*.csv` y `22_*.png`.

## Config (`R/00_config.R`)
- Eliminadas las claves que solo usaba el flujo antiguo (`n_folds`, `grilla_lambda`, `nivel_intervalo`, `umbral_mantencion`, `subgrupo_ref`, `ruido_niveles`, `dir_resultados`).
- `n_hilos = 1` (ranger y xgboost dependen del número de hilos); `rf_arboles = 500` (igual que los pasos 17 y 20).
- Nuevas: `semilla_ruido` (mismo valor 20261031 que estaba fijo en el paso 15, así los números no cambian), `n_pdp_variables`, `pdp_puntos`; `n_perm` sube a 10.
- Los comentarios sobre los test ya no dicen "sellados" ni "una sola vez": dicen que se examinaron antes con el pipeline histórico.

## Reproducibilidad
- `R/01_utilidades.R`: `verificar_referencia_ranking()` reemplaza a la función antigua. Compara con la referencia **solo si la huella de la configuración coincide**; si cambió la configuración, crea una referencia nueva y lo avisa. La referencia existente (sin huella) se reemplaza en la próxima ejecución.
- `pasos/11`: usa esa función.
- Entreguen `codigo/referencia/` generada **después** del último cambio (ejecutar dos veces).

## Parametrización y coherencia
- Rutas fijas reemplazadas por `cfg$dir_final`, `cfg$dir_anidada` o `DIR_TABLAS` en los pasos 14, 15, 17, 18 y 19. El paso 14 ahora carga la config por sí mismo.
- Pasos 12 y 16: se eliminaron los conteos fijos (8486, 2122, 1326); ahora se verifica coherencia con `cfg$prop_test` y que haya una predicción LOSO por fila.
- Paso 17: la comprobación de convergencia se hace sobre el **entrenamiento** (antes usaba los test). La grilla de λ se lee de `cfg$grilla_lambda_loso`, y hay un `warning` si λ queda en el borde inferior.
- `R/03_particion.R`: comentarios corregidos y eliminada `crear_folds_aleatorios()` (no se usaba).
- `00_main.R`: el registro lee el mejor modelo del ranking anidado en vez de escribirlo a mano, y avisa si no es el Elastic-net.

## Movido a `historico/modulos_antiguos/`
`07_comparacion.R`, `08_test.R`, `09_interpretacion.R` (flujo antiguo con CV aleatoria; `08` usa `CONTADOR_TEST`, que no está definido).

## Sin cambiar (decisiones de ustedes)
- `anidar_ensembles = FALSE`: RF, XGBoost y árbol siguen con configuración fija. Cambiarlo a `TRUE` agrega 30 a 60 minutos y cambia la huella.
- Paso 15 sigue usando los test para la sensibilidad a ruido (solo diagnóstico, no selecciona nada).
- Rejilla de λ sin ampliar: primero revisen `12_parametros_elegidos.csv`.
- Los pasos 10, 13, 20 y 21, y los módulos 02, 04, 05 y 06, no se tocaron.

---

# Revisión 2 (ejecutada y verificada)

El proyecto se ejecutó completo dos veces desde una sesión limpia (R 4.3.3; versiones en
`resultados_elastic_net_final/versiones_paquetes.csv`). La segunda ejecución reproduce el
ranking anidado (diferencia máxima 4e-15).

## Reproducibilidad
- `datos/data.txt` incluido, con verificación de contenido (`cfg$filas_datos`, `cfg$firma_datos`).
- `README.md` con instalación, ejecución y origen de cada cifra; `equipo_regresion_naval.Rproj`.
- `codigo/referencia/` generada con la configuración final.
- XGBoost con `xgb.train` y `tree_method`/`base_score` explícitos: el resultado ya no depende
  de la versión de xgboost (1.7 o 3.x).
- Los pasos 11 y 12 guardan solo los hiperparámetros (antes guardaban también columnas de métricas).

## Sin fuga
- Nuevo paso 23: CV aleatoria (vfold_cv) vs LOSO solo con el entrenamiento. Reemplaza las
  cifras del §3 que venían del pipeline histórico y del test de 15 nudos.

## Materia del curso
- Paso 20: GAM (mgcv, s() con bs = "cr", REML dentro de cada partición) evaluado con LOSO; entra en la Tabla 1.
- Paso 10: IC clásico de confint() junto al IC con HC3.
- Comentarios sobre los supuestos del bootstrap (paso 12), de la t pareada (paso 21), del
  intervalo de predicción (paso 16) y de la validación anidada (paso 11).

## Legibilidad y coherencia
- Eliminado el código muerto de 01, 02 y 06. `diagnostico_ols` ahora se usa en el paso 18
  (residuos, Cook y VIF).
- Paso 18: se eliminaron los mensajes que contradecían al paso 16.
- Paso 17: se eliminó la comparación OLS/RF/XGBoost duplicada (queda solo en el paso 20).
- Pasos 14, 15, 17, 18 y 19: encabezado estándar, `guardar_tabla` y registro en el log.
- Paso 22: etiquetas con unidades y nota sobre las predicciones fuera de rango en PDP/ICE.

---

# Revisión 3 (respuesta a la segunda auditoría; ejecutada y verificada)

Dos ejecuciones completas e independientes desde sesiones limpias (82 min cada una): las 40
tablas son idénticas y la verificación del paso 11 reproduce la referencia con diferencia 0.

## Sin fuga
- Paso 15: la sensibilidad a ruido se calcula con LOSO sobre el entrenamiento (antes usaba los
  dos test 100 veces cada uno). Los test solo se usan en la evaluación única del paso 12.

## Rigor
- `anidar_ensembles = TRUE`: árbol, Random Forest y XGBoost eligen hiperparámetros por LOSO
  anidado con la regla 1-SE. Pasan a la Tabla 1 como "LOSO anidada (1-SE interno)".
  XGBoost usa `tree_method = "hist"` (igual resultado en xgboost 1.7 y 3.x).
- Paso 10: fila adicional con el efecto de +10 y su IC clásico (confint).

## Legibilidad y coherencia
- `crear_loso()`, `metricas_fila()` y `rmse_vec()` viven en `01_utilidades.R` (antes estaban
  repetidas en los pasos 11, 12, 17, 18, 19, 20 y 23).
- Valores fijos a la configuración: `B_ruido`, `rf_arboles` (ahora sí se usa), `arbol_fijo` y `xgb_fijo`.
- Pasos 14, 17, 18 y 19: `cat`/`print` reemplazados por `log_msg`; tildes en los comentarios.
- `huella_cfg()` ignora `grilla_lambda`, que los pasos copian en tiempo de ejecución: la huella
  depende solo de `00_config.R`.
- Diccionario: temperaturas en K (T1 = 288 constante solo es físico en Kelvin; +10 K = +10 °C);
  correlación lp–v = 0,99992 en vez de "1". Comentarios de rutas y de `rf_arboles` corregidos.
- Paso 16: "intervalos conformales" pasa a llamarse "cuantil empírico de residuos LOSO".

## Verificación final
Tercera ejecución completa con la referencia ya incluida: `resultados_validacion_anidada/log.txt`
registra "[OK] el ranking anidado reproduce la ejecución de referencia (dif. máx. = 4e-15)" y las
45 tablas son idénticas a las de las ejecuciones anteriores. Los resultados incluidos son de esa corrida.

---

# Revisión 4: multicolinealidad y valores SHAP

- Nuevo paso 24 (solo entrenamiento, no cambia ningún resultado anterior ni la configuración):
  matriz de correlaciones entre predictores (`cor()`, como en clases), gráfico de VIF (del paso 18)
  y valores SHAP exactos del Elastic-net final: en un modelo lineal, SHAP = coeficiente ×
  (valor estandarizado − media). Se verifica la aditividad (predicción media + suma de los SHAP =
  predicción) con tolerancia 1e-8.

---

# Revisión 5: curva de error de entrenamiento vs validación

- Nuevo paso 25 (solo entrenamiento y salidas del paso 12): para cada λ de la grilla, con el α final,
  RMSE dentro de muestra del Elastic-net ajustado con todo el entrenamiento frente al RMSE LOSO
  (media ± 1 EE) del paso 12. Los test aparecen solo en el λ elegido, leídos del paso 12: no se
  evalúa ningún λ nuevo en test. Verifica que el RMSE de entrenamiento con el λ elegido coincide con
  el del paso 17 (1,029).

---

# Revisión 6: ejecución oficial del equipo

- Ejecución completa en el computador del equipo (Windows 11, R 4.6.1, tune 2.1.0, xgboost 3.2.1.1;
  58,7 min, pasos 10 a 25). Sus resultados reemplazan a los anteriores y son la referencia.
- Frente a la ejecución en Linux (R 4.3.3, xgboost 1.7.7.1), todas las tablas coinciden salvo XGBoost
  (tercer decimal) y las simulaciones de ruido y permutación (tercer o cuarto decimal).
- `huella_cfg()` ahora se calcula con los valores de la configuración, no con `str()`, cuyo texto
  cambiaba entre versiones de R y hacía que el paso 11 creara una referencia nueva en vez de comparar.

---

# Revisión 7: figuras de la presentación dentro del proyecto

- Nuevo paso 26: genera las imágenes de la presentación (láminas 4, 6, 7, 8 y 10) a 350 dpi desde
  las tablas que ya produce el script; no calcula resultados nuevos. Verifica que la pendiente de T2
  de la lámina 4 coincide con la del paso 10. Requiere el paquete gridExtra (agregado a la lista).
