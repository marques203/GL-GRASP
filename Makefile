# Makefile para o C-IGDP (GRASP/TABU para o problema IMLCM)
#
# Build padrão: apenas os modos que não exigem SDK pago
#   (grasp1, grasp2, grasp3, tabu)
#
#   make            -> compila ./C-IGDP.exe
#   make run ARGS=... -> compila (se preciso) e executa com os argumentos dados
#   make clean      -> remove o binário e objetos
#
# Para habilitar os modos "cplex" e "lsolver" (exige IBM CPLEX / Localsolver
# instalados e configurados via CPLEX_DIR / LOCALSOLVER_DIR):
#   make WITH_CPLEX=1 WITH_LOCALSOLVER=1

CXX      := g++
CXXFLAGS := -std=c++14 -O2 -w
LDFLAGS  := -static

TARGET   := C-IGDP.exe

SRCS := main.cpp HDAG.cpp GRASP.cpp GRASPv1.cpp GRASPv2.cpp GRASPv3.cpp \
        TABU.cpp EliteSet.cpp PathRelinking.cpp SolutionIMLCM.cpp SolutionMLCM.cpp

ifdef WITH_CPLEX
CXXFLAGS += -DWITH_CPLEX -I"$(CPLEX_DIR)/include"
LDFLAGS  += -L"$(CPLEX_DIR)/lib" -lilocplex -lcplex -lconcert -lm -lpthread
SRCS     += Cplex.cpp
endif

ifdef WITH_LOCALSOLVER
CXXFLAGS += -DWITH_LOCALSOLVER -I"$(LOCALSOLVER_DIR)/include"
LDFLAGS  += -L"$(LOCALSOLVER_DIR)/bin" -llocalsolver
SRCS     += Localsolver.cpp
endif

OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) $(ARGS)

clean:
	rm -f $(OBJS) $(TARGET)
