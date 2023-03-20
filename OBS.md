OBS:
- A nivel de especie:
- Con filtrado: 14 
- Sin filtrado: 24
A nivel de genero:
- Con filtrado: 21
- Sin filtrado: 37

- Si se avanza con las especies filtradas se pierden bastantes especies (13) de las que se habían utilizado para ela clusterización y 3 de las 8 para los modelos de clasificación.

- *** Tiene sentido hacer el filtrado de presencia teniendo en cuenta que tenemos una cohorte donde están todas enfermas en el D0 (PRJNA3020) y seguramente especies como las Lactobacillus no van a tener presencia y se van a perder en el filtrado. Por consiguiente se van a perder en el intersect con el resto de cohortes(pero dichas especies si van a ser relevantes en esas cohortes donde hay sanos)

- Se podría hacer un merge con todos los phyloseqs (solo con las especies compartidas) y ahí hacer el filtrado?

CLUST:
Original Ravel
                ID   D   H
  high          18  78   1
  intermediate  38   6   5
  low          121   5 122

Original Sriniv
               ID  D  H
  high         28 88  1
  intermediate  8  2  3
  low          41  0 49

- Usando 14 especies:
* 5 proposed 2 as the best number of clusters 
* 12 proposed 3 as the best number of clusters **** de 15 a 12 ****
* 3 proposed 4 as the best number of clusters 
* 5 proposed 5 as the best number of clusters 
* 2 proposed 6 as the best number of clusters 
* According to the majority rule, the best number of clusters is  3 

Ravel
                 D  ID   H      
  high          73  12  12
  intermediate   7  25  17
  low            3  94 151

Sriniv
                D ID  H
  high         89  8 20
  intermediate  3  3  7
  low           1 28 61

PR3020
 D ID  H 
48 10 10 

OBS: Podemos ver como reparte mas de High e I en la clase que se correspondería con low

- Usando 24 especies:                                              
* 3 proposed 2 as the best number of clusters 
* 6 proposed 3 as the best number of clusters 
* 8 proposed 4 as the best number of clusters 
* 6 proposed 5 as the best number of clusters 
* 4 proposed 6 as the best number of clusters 
* According to the majority rule, the best number of clusters is  4 

Ravel 3C
                ID   H   D
  high          22   1  74
  intermediate  37   6   6
  low          120 124   4

Sriniv 3C
               ID  H  D
  high         29  1 87
  intermediate  9  3  1
  low          40 50  0

PR3020 3C
ID D 
18 50 

Ravel 4C
                 D  ID1  H ID2
  high          70  15   1  11
  intermediate   5  23   2  19
  low            3  29 117  99

Sriniv 4C
                D ID1 H ID2
  high         80 12  1 24
  intermediate  1  2  2  8
  low           1  1 49 39

PR3020 4C
D ID1 ID2 
42 10 16 

OBS: Aunque no nos den 3C como lo recomendado si hacemos los clusters con las 24 especies tenemos resultados similares (en terminos de numeros) de los cluster anteriores, bailan algunas samples al cluster de ID.