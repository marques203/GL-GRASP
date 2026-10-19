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
        GLGRASP.cpp EmbeddingDistances.cpp \
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

# Dependencia dos headers. -MMD faz o compilador escrever, ao lado de cada
# .o, um .d com a lista de headers do projeto que aquele .cpp inclui; -MP
# evita erro se um header for apagado. Sem isso, editar um .h (por exemplo
# acrescentar um membro em GRASP.h) recompilava so os .cpp editados, e o
# binario misturava duas versoes da mesma classe -- ja causou segfault.
DEPS := $(OBJS:.o=.d)

TEST_TARGET := tests/test_glgrasp_units.exe

# objetos de biblioteca (tudo menos o main) para linkar com os testes
LIB_OBJS := $(filter-out main.o,$(OBJS))

.PHONY: all clean run test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

run: $(TARGET)
	./$(TARGET) $(ARGS)

# Testes de unidade das pecas novas do GL-GRASP (itens 1.1 e 1.2).
test: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TEST_TARGET): tests/test_glgrasp_units.cpp $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ tests/test_glgrasp_units.cpp $(LIB_OBJS) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET) $(TEST_TARGET)
