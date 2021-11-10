# Compilateur utilisé
CC=mpicxx

CCseq=g++
#g++

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -O3 -DNDEBUG -std=c++11
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -std=c++11

# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler

SRC = main.cpp GradConj.cpp Problem.cpp BC.cpp Output.cpp Readfile.cpp

SRC1 = main.cpp GradConj.cpp Problem.cpp BC.cpp Output.cpp Readfile.cpp

# La commande complète : compile le séquentiel
$(PROG) : $(SRC)
	$(CCseq) $(SRC) $(CXX_FLAGS) -o $(PROG)
	$(warning "Attention, pas de mode de compilation renseigné, compilation du code séquentiel")
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

seq: $(SRC)
	$(CCseq) $(SRC) $(CXX_FLAGS) -o $(PROG)
# à changer le compilateur pour le seq

par: $(SRC1)
	$(CC) $(SRC1) $(CXX_FLAGS) -o $(PROG)


#éxecution adaptée
run_seq:
	./run data_file.txt

run_par:
	time mpiexec --oversubscribe -n $(N) ./run data_file.txt



# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
