# BV_Microbiomics
Estado:
Al parecer tenemos una nueva clasificación con 4 clusters que van en sintonía con las otras cohortes.

Siguientes pasos:
- Asegurarnos de la clusterización utilizando otro método como **ConsensusClusterPlus** (k-means y jerárquico) y diferentes métodos de **normalización** (clr, alr y log2). A partir de ahí, definir los clusters que nos salen para todo el análisis posterior.
- Con las etiquetas nuevas correr los siguientes análisis (con todas las especies de cada cohorte):
  - **Co-abundancias** mediante SpiecEasi en la cohorte de Ravel
  - **Meta-análisis** de DA en todas las cohortes que tengamos disponibles. En caso de no poder hacer DA multi-class hacer una comparación de los clusters contra healthy (en las cohortes donde tengamos healthy)