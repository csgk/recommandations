Recommandations
===============

#Données
Elles doivent être sous la forme :
personne	objet	note	autre
personne	objet	note	autre
personne	objet	note	autre

#Compilation
Sous unix, il faut compiler le fichier main.cpp à l'aide de la commande :
	
	g++ main.cpp
	
#Drapeaux
On trouvera les drapeaux obligatoires :

* `-input-file file`	avec file le nom du fichier de données
* `-persons m`		avec m le nombre de personnes
* `-objects n`		avec n le nombre d'objets
* `-entries k`		avec k tel que 1 note sur k sera consacrée à l'évaluation de l'algorithme

Et les drapeaux optionnels :

* `-convex l`		avec l le coefficient lambda défini dans le rapport (par défaut, lambda = 10^5)
* `-neighb k`		avec k le nombre de voisins (par défaut, k = 50)
* `-svd a`		avec a le nombre de valeurs singulières (par défaut, a = 10)
* `--details`		qui permet d'afficher les erreurs en fonction du nombre de votes par personne et objet

#Exécution
Par exemple, on pourra tester 1/5 des données sur le fichier u.data de MovieLens 100k à l'aide de la commande suivante :
	
	./a.out -input-file u.data -persons 943 -objects 1682 -entries 5
	
