
# -Wconversion

CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
ifeq ("$(wildcard $(CXX_INST_DIR))","")
  SUFFIX = 
  CXX_INST_DIR = $(HOME)/bin
  ifeq ("$(wildcard $(CXX_INST_DIR))","")
    ifneq ($(wildcard "/mingw64"),"")
      CXX_INST_DIR = /mingw64
    else
      CXX_INST_DIR = /usr
    endif
  endif
endif

GFORTRAN = $(CXX_INST_DIR)/bin/gfortran -g -Wall -Wextra -Wno-compare-reals
GCC = $(CXX_INST_DIR)/bin/gcc -g -Wall -Wextra
CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi
CXX20 = $(CXX_INST_DIR)/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi
CXXMAX = $(CXX20)
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/7.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = .

INPUT_DIR = test_input
OUTPUT_DIR = test_output
BIN_DIR = bin

BINS = $(BIN_DIR) \
  $(BIN_DIR)/test_bairstow \
  $(BIN_DIR)/test_horner \
  $(BIN_DIR)/test_jenkins_traub \
  $(BIN_DIR)/test_polynomial \
  $(BIN_DIR)/test_polynomial_root \
  $(BIN_DIR)/test_rational_polynomial \
  $(BIN_DIR)/test_solution \
  $(BIN_DIR)/test_solvers \
  $(BIN_DIR)/test_static_polynomial \
  $(BIN_DIR)/test_laguerre_step \
  $(BIN_DIR)/test_quadratic_step \
  $(BIN_DIR)/test_jacobi_roots



all: $(BINS)

ext/polynomial.h: include/ext/polynomial.tcc

ext/solver_low_degree.h: include/ext/solution.h include/ext/solver_low_degree.tcc

ext/rational_polynomial.h: include/ext/polynomial.h

$(BIN_DIR)/test_bairstow: test_bairstow.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_bairstow test_bairstow.cpp -lquadmath

$(BIN_DIR)/test_horner: include/ext/horner.h test_horner.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_horner test_horner.cpp -lquadmath

$(BIN_DIR)/test_jacobi_roots: include/ext/solution.h include/ext/solver_low_degree.h include/ext/polynomial.h test_jacobi_roots.cpp
	$(CXXMAX) -Iinclude -I../include -o $(BIN_DIR)/test_jacobi_roots test_jacobi_roots.cpp -lquadmath

$(BIN_DIR)/test_jenkins_traub: include/ext/solution.h include/ext/solver_low_degree.h include/ext/polynomial.h test_jenkins_traub.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_jenkins_traub test_jenkins_traub.cpp -lquadmath

$(BIN_DIR)/test_polynomial: include/ext/polynomial.h include/ext/polynomial.tcc test_polynomial.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_polynomial test_polynomial.cpp -lquadmath

$(BIN_DIR)/test_polynomial_root: include/ext/polynomial.h test_polynomial_root.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_polynomial_root test_polynomial_root.cpp -lquadmath

$(BIN_DIR)/test_rational_polynomial: include/ext/rational_polynomial.h test_rational_polynomial.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_rational_polynomial test_rational_polynomial.cpp -lquadmath

$(BIN_DIR)/test_solvers: include/ext/solution.h include/ext/solver_low_degree.h include/ext/solver_low_degree.tcc test_solvers.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_solvers test_solvers.cpp -lquadmath

$(BIN_DIR)/test_solution: include/ext/solution.h test_solution.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_solution test_solution.cpp -lquadmath

$(BIN_DIR)/test_static_polynomial: include/ext/static_polynomial.h test_static_polynomial.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_static_polynomial test_static_polynomial.cpp -lquadmath

$(BIN_DIR)/test_laguerre_step: test_laguerre_step.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_laguerre_step test_laguerre_step.cpp -lquadmath

$(BIN_DIR)/test_quadratic_step: test_quadratic_step.cpp
	$(CXXMAX) -Iinclude -o $(BIN_DIR)/test_quadratic_step test_quadratic_step.cpp -lquadmath

docs:
	rm -rf docs/html/*
	rm -rf docs/latex/*
	doxygen
	cd docs/latex && make

test: $(OUTPUT_DIR)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rational_polynomial > $(OUTPUT_DIR)/test_rational_polynomial.txt
	touch $(OUTPUT_DIR)/test_bairstow.prev && cp -f $(OUTPUT_DIR)/test_bairstow.prev $(OUTPUT_DIR)/test_bairstow.prev2
	touch $(OUTPUT_DIR)/test_bairstow.txt && cp -f $(OUTPUT_DIR)/test_bairstow.txt $(OUTPUT_DIR)/test_bairstow.prev
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_bairstow.in > $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver1.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver2.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver3.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver4.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver5.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver6.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver7.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver8.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_horner > $(OUTPUT_DIR)/test_horner.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver1.in > $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver2.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver3.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver4.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver5.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver6.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver7.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jenkins_traub < $(INPUT_DIR)/test_solver8.in >> $(OUTPUT_DIR)/test_jenkins_traub.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_polynomial > $(OUTPUT_DIR)/test_polynomial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_polynomial_root > $(OUTPUT_DIR)/test_polynomial_root.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_solvers > $(OUTPUT_DIR)/test_solvers.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_static_polynomial > $(OUTPUT_DIR)/test_static_polynomial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_laguerre_step > $(OUTPUT_DIR)/test_laguerre_step.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_quadratic_step > $(OUTPUT_DIR)/test_quadratic_step.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_jacobi_roots > $(OUTPUT_DIR)/test_jacobi_roots.txt


$(OUTPUT_DIR): $(OUTPUT_DIR)
	if test ! -d $(OUTPUT_DIR); then \
	  mkdir $(OUTPUT_DIR); \
	fi


$(BIN_DIR):
	if test ! -d $(BIN_DIR); then \
	  mkdir $(BIN_DIR); \
	fi


clean:
	rm -f a.out
	rm -f *.stackdump
	rm -f $(BINS)
