CC = g++
CFLAGS = -std=c++17

INCLUDE = /Applications/MATLAB_R2019b.app/extern/include/ -I/Applications/MATLAB_R2019b.app/simulink/include
LIB1 = /Applications/MATLAB_R2019b.app/bin/maci64
LIB2 = /Applications/MATLAB_R2019b.app/extern/bin/maci64

INCLUDE_EIGEN = ../Eigen
INCLUDE_SPECTRA = ../Spectra/include/
# INCLUDE_GNUPLOT = /Users/dengxiaohui/Documents/gnuplot-5.2.8/gnuplot-iostream/

INCLUDE_DIR = ./include
SRC_DIR = ./src
APP_DIR = ./app
OBJ_DIR = ./obj


objects_file = $(OBJ_DIR)/Define.o $(OBJ_DIR)/DifferenceOperator.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/MatrixAssembly.o $(OBJ_DIR)/write_matrix.o

write_matrix : $(objects_file)
	$(CC) $(objects_file) -o write_matrix

$(OBJ_DIR)/write_matrix.o	:	$(APP_DIR)/write_matrix.cpp
	$(CC) -c $(APP_DIR)/write_matrix.cpp -o $(OBJ_DIR)/write_matrix.o -I$(INCLUDE_EIGEN) -I$(INCLUDE_SPECTRA) -I$(INCLUDE_DIR)

# if work on server, link parameters should be changed

object_eigensolver_dense_file = $(OBJ_DIR)/Define.o $(OBJ_DIR)/DifferenceOperator.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/MatrixAssembly.o $(OBJ_DIR)/Eigen_Lapack_Interface.o

generalized_eigensolver_dense	:	$(object_eigensolver_dense_file) $(OBJ_DIR)/generalized_eigensolver_dense.o
	$(CC) $(object_eigensolver_dense_file) $(OBJ_DIR)/generalized_eigensolver_dense.o -o generalized_eigensolver_dense -llapack -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)

$(OBJ_DIR)/generalized_eigensolver_dense.o	:	$(APP_DIR)/generalized_eigensolver_dense.cpp
	$(CC) -c $(APP_DIR)/generalized_eigensolver_dense.cpp -o $(OBJ_DIR)/generalized_eigensolver_dense.o -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)

object_eigensolver_band_file = $(OBJ_DIR)/Define.o $(OBJ_DIR)/DifferenceOperator.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/MatrixAssembly.o $(OBJ_DIR)/Eigen_Lapack_Interface.o

generalized_eigensolver_band	:	$(object_eigensolver_band_file) $(OBJ_DIR)/generalized_eigensolver_band.o
	$(CC) $(object_eigensolver_band_file) $(OBJ_DIR)/generalized_eigensolver_band.o -o generalized_eigensolver_band -llapack -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)

$(OBJ_DIR)/generalized_eigensolver_band.o	:	$(APP_DIR)/generalized_eigensolver_band.cpp
	$(CC) -c $(APP_DIR)/generalized_eigensolver_band.cpp -o $(OBJ_DIR)/generalized_eigensolver_band.o -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)


####################################
#	Shared object files
$(OBJ_DIR)/Define.o	:	$(SRC_DIR)/Define.cpp $(INCLUDE_DIR)/Define.h
	$(CC) -c $(SRC_DIR)/Define.cpp -o $(OBJ_DIR)/Define.o -I$(INCLUDE_DIR)

$(OBJ_DIR)/DifferenceOperator.o	:	$(SRC_DIR)/DifferenceOperator.cpp $(INCLUDE_DIR)/DifferenceOperator.h
	$(CC) -c $(SRC_DIR)/DifferenceOperator.cpp -o $(OBJ_DIR)/DifferenceOperator.o -I$(INCLUDE_DIR)

$(OBJ_DIR)/FileIO.o	:	$(SRC_DIR)/FileIO.cpp $(INCLUDE_DIR)/FileIO.h
	$(CC) -c $(SRC_DIR)/FileIO.cpp -o $(OBJ_DIR)/FileIO.o -I$(INCLUDE_DIR) -I$(INCLUDE_EIGEN)

$(OBJ_DIR)/MatrixAssembly.o	:	$(SRC_DIR)/MatrixAssembly.cpp $(INCLUDE_DIR)/MatrixAssembly.h
	$(CC) -c $(SRC_DIR)/MatrixAssembly.cpp -o $(OBJ_DIR)/MatrixAssembly.o -I$(INCLUDE_DIR) -I$(INCLUDE_EIGEN)

####################################


# if work on server, link parameters should be changed
test_lapack	:	$(OBJ_DIR)/Eigen_Lapack_Interface.o $(OBJ_DIR)/test_lapack.o
	$(CC) $(OBJ_DIR)/Eigen_Lapack_Interface.o $(OBJ_DIR)/test_lapack.o -o test_lapack -llapack -I$(INCLUDE_DIR) -I$(INCLUDE_EIGEN)

$(OBJ_DIR)/Eigen_Lapack_Interface.o	:	$(SRC_DIR)/Eigen_Lapack_Interface.cpp 
	$(CC) -c $(SRC_DIR)/Eigen_Lapack_Interface.cpp -o $(OBJ_DIR)/Eigen_Lapack_Interface.o -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)

$(OBJ_DIR)/test_lapack.o :	$(APP_DIR)/test_lapack.cpp
	$(CC) -c $(APP_DIR)/test_lapack.cpp -o $(OBJ_DIR)/test_lapack.o -I$(INCLUDE_EIGEN) -I$(INCLUDE_DIR)




cleanobj	:
	rm $(objects_file)

clean :
	rm write_matrix $(objects_file)