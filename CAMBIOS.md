# Estado final del proyecto · Equipo 5

- **Ejecución:** `codigo/00_main.R` corre los pasos 10 a 26 con una sola orden; todos los parámetros están en `codigo/R/00_config.R`.
- **Datos:** UCI *Condition Based Maintenance of Naval Propulsion Plants* (11.934 filas); al cargarlos se verifica que sean la versión usada en el informe.
- **Partición:** test de velocidad no vista (15 nudos, 1.326 filas), test aleatorio (20 %, 2.122 filas) y entrenamiento (8.486 filas, 8 velocidades). Semilla 20261016.
- **Validación:** dejando una velocidad completa fuera (LOSO), anidada; los hiperparámetros se eligen con la regla 1-SE en la validación interna. Estandarización, imputación y nudos de splines se estiman dentro de cada fold mediante recetas.
- **Modelos:** baseline, OLS, Ridge, Lasso, Elastic-net, splines naturales, GAM, árbol, Random Forest, XGBoost y promedio EN + Ridge. Sistema final: Elastic-net (λ = 0,0060; α = 0,5).
- **Test:** evaluación única del modelo congelado (paso 12). Ninguna decisión (modelo, hiperparámetros, intervalos ni umbral) usa los test.
- **Interpretación:** coeficientes, importancia por permutación, PDP/ICE, multicolinealidad (correlaciones y VIF), valores SHAP y curva de error de entrenamiento frente a validación.
- **Verificaciones automáticas:** versión de los datos, particiones sin filas compartidas, exclusión de variables no disponibles al estimar (kMt), aditividad de los SHAP y reproducibilidad del ranking anidado (tolerancia 1e-3).
- **Ejecución oficial:** Windows 11, R 4.6.1, 60,8 min. Reproduce la referencia con una diferencia máxima de 4e-15; dos ejecuciones del equipo dieron tablas idénticas. Versiones de paquetes en `resultados_elastic_net_final/versiones_paquetes.csv`.
