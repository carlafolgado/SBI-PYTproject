# SBIPYTproject
Per configurar github, obriu un terminal i feu (el $ no l'escriviu, és per indicar terminal):

```{sh} 
$ git config --global user.name "Nom Cognom"
$ git config --global user.email "correu"
```

## Clonar el repositori
Des del terminal o el Finder/Finder per windows, creeu una carpeta allà on vulgueu (en aquesta carpeta hi haurà tot el treball (es pot dir com vulgueu).Un cop estigueu dins de la carpeta en qüestió, heu de clonar el repositori:
```{sh} 
$ git clone https://github.com/carlafolgado/SBIPYTproject.git . <---(el punt és important, no us el deixeu o no funcionarà)
```
Amb això us descarregareu tot el que estigui a GitHub en aquell moment, i a més enllaçareu la vostra carpeta (repositori) local amb el repositori online a GitHub.

## Actualitzar el repositori
La idea de git és que tu edites els arxius com et dóna la gana, amb els programes que vols. I després, amb git fas "commits", és a dir, guardes els canvis conforme els vas fent (si no m'he explicat, m'ho pregunteu en persona). Quan estigueu fent canvis, sempre és útil tenir una finestra de terminal oberta al repositori. I abans de fer res, és recomanable fer
```
$ git fetch
$ git status
```
que us dirà si hi ha canvis nous online que encara no teniu a l'ordinador. Si tot està bé, l'ordinador no us cridarà i podreu fer
```{sh}
$ git merge
```
per incorporar els canvis. Si esteu segurs de que no hi ha conflictes podeu fer ```git pull```, que és com fer git fetch i git merge tot alhora.

## Editar, fer commits i push
Aleshores vosaltres anireu editant amb el programa que més us agradi. Quan hagueu fet prous canvis que tingui sentit agrupar-los (o no, en realitat podeu fer això quan vulgueu), és moment de fer un commit. Per això sempre va bé fer
```
$ git status
```
per saber com teniu el repositori. Tindreu una llista d'arxius que us sortirà com "Changes not staged for commit", o alguna cosa així. Això vol dir que git ha detectat canvis, però que de moment no fa res. Aleshores feu
```
$ git add fitxer1 fitxer2 etc
```
Amb això li esteu dient a git que us interessa guardar els canvis d'aquests fitxers. Amb
```
$ git .
```
li dieu que us interessen tots i us estalvieu afegir-los tots un a un. I aleshores, si esteu segurs (no està de més tornar a fer git status, i haurieu de veure els arxius que li heu dit com a "Changes to be commited" o alguna cosa així), heu de fer commit:
```
$ git commit -m "Missatge"
```
Amb això li dieu a git que guardi els canvis que li heu dit (quan dic guardar vull dir guardar dins de git. Òbviament els canvis els teniu a l'rdinador. Per cada commit, git guarda una versió del repositori de manera que sempre pots tornar enrer). Siusplau, poseu missatges més o menys intel·ligibles on digueu què heu fet, que va bé quan hi ha 150 commits. Oju, que aquests canvis de moment només estan al vostre repositori local. Per enviar-los a interné, heu de fer un push:
```
$ git push
```
i amb això envieu els vostres canvis per les autovies d'internet. Els altres, si ho han entès, quan es posin a treballar, abans de res actualitzaran el seu repositori (véase la secció anterior), de manera que així s'asseguren que tenen la versió més actualitzada del repositori, amb els canvis que vosaltres acabeu de publicar.

Quan feu git push per primera vegada us demanarà l'usuari i la contrasenya de GitHub. Si la poseu bé no us hauria de donar problemes.

## Aplicació
Si tot això se us fa bola, hi ha una aplicació de GitHub amb la que podeu fer tot això sense el terminal. El procés és el mateix (pull -> add -> commit -> push), però us estalvieu el terminal. No sé com funciona massa bé, però m'imagino que hi haurà tutorials per internet i tal.
