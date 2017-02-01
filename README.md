# BE_TSTR

Préambule:
Notre travail n'a pas abouti pour le domaine fréquentiel et en plus a provoqué un dysfonctionnement dans le domaine temporel. Nous avons donc décidé de fournir deux fichiers différents et de créer un exécutable pour chaque fichier.
Dans reverb_temporel.cpp, il n'y a que la convolution temporelle (qui fonctionne).
Dans reverb.cpp, il y a les deux méthodes mais aucunes ne marchent.

CONFIGURER
Pour configurer le projet, aller dans le dossier rtaudio-4.2.1 et taper:
./configure

COMPILER
Une fois la configuration finie, il faut aller rajouter dans le /tests/Makefile, toutes les lignes concernants la compilation du fichier reverb.cpp et reverb_temporel.cpp.
Pour cela il faut copier tous les endroits contenant le mot duplex et remplacer le mot duplex par reverb. De même, pour reverb_temporel.
Pour la ligne ci-dessou, il ne faut pas la recopier ni remplacer la copie par reverb ou reverb_temporel:
include ./$(DEPDIR)/duplex.Po
Pour compiler le projet, aller dans le dossier tests/ et taper make.

PROGRAMME reverb.cpp
Pour avoir la réponse impulsionnelle du filtre en entier commenter:
#define SHORT_FILTER //mettre en commentaire pour avoir la réponse impulsionnelle complète

Si cette ligne n'est pas commentée, il est possible de changer la taille du filtre à la ligne 587
Pour changer la taille du buffer du stream, modifier la ligne 601

Pour avoir la convolution dans le domaine temporelle, il faut commenter la ligne dans reverb.cpp 
#define DOMAINE_FREQUENTIEL //mettre en commentaire pour avoir la convolution add dans le domaine temporelle
