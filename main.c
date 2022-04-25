#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define TAU (0.000002)
#define MATRIX_SIZE 1000
#define EPS 0.000001

void testVector(const double* vector){
    printf("vector:\n");
    for (int i = 0; i < MATRIX_SIZE; i++){
        printf("%.2f\n", vector[i]);
    }
}



void mulMatrixOnVector(const double* matrix, const double* vector, double* answer){
    double* tmp = (double*) malloc(sizeof(double)*MATRIX_SIZE);
#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++){
        tmp[i] = 0;
        for (int k = 0; k < MATRIX_SIZE; k++){
            tmp[i] += matrix[i*MATRIX_SIZE+k]*vector[k];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++){
        answer[i] = tmp[i];
    }
    free(tmp);
}


void subtractionVectors(const double* vector1, const double* vector2, double* answer){
#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++){
        answer[i] = vector1[i] - vector2[i];
    }
}

void mulVectorOnConst(const double* vector, double* answer){
#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++){
        answer[i] = vector[i]*TAU;
    }
}

double vectorNorm(const double* vector){
    double tmp = 0;
#pragma omp parallel for reduction(+:tmp)
    for (int i = 0; i < MATRIX_SIZE; i++){
        tmp += pow(vector[i], 2);
    }
    return sqrt(tmp);
}


int isFinish(const double* vector1, const double normVectorB){
    if ((vectorNorm(vector1) / normVectorB) < EPS){
        return 0;
    }
    else return 1;
}

// (Ax-b)
void makeVectorDiff(double* matrix, double* vectorX, double* vectorB, double* answer){
    mulMatrixOnVector(matrix, vectorX, answer);
    subtractionVectors(answer, vectorB, answer);
}


void algorithm(double* matrix, double* vectorX, double* vectorB){
    double *vectorDiff = (double*) malloc (sizeof(double)*MATRIX_SIZE);

    makeVectorDiff(matrix, vectorX, vectorB, vectorDiff);

    double normVectorB = vectorNorm(vectorB);
    while(isFinish(vectorDiff, normVectorB) == 1){
        mulVectorOnConst(vectorDiff, vectorDiff);
        subtractionVectors(vectorX, vectorDiff, vectorX);
        makeVectorDiff(matrix, vectorX, vectorB, vectorDiff);
    }
    free(vectorDiff);
}

//
//void initMatrix(double* matrix, double *vectorX, double *vectorB){
//    for (int i = 0; i < MATRIX_SIZE; i++){
//        for (int j = 0; j < MATRIX_SIZE; j++){
//            matrix[i*MATRIX_SIZE+j] = i==j ? 2: 1;
//        }
//        vectorX[i] = sin(2* M_PI/MATRIX_SIZE);
//    }
//    testVector(vectorX);
//    mulMatrixOnVector(matrix, vectorX, vectorB);
//}

void initMatrix(double* matrix, double *vectorX, double *vectorB){
#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++){
        for (int j = 0; j < MATRIX_SIZE; j++){
            if (i == j){
                matrix[i*MATRIX_SIZE+j] = 2.0;
            }
            else matrix[i*MATRIX_SIZE+j] = 1.0;
        }
        vectorX[i] = 0;
        vectorB[i] = MATRIX_SIZE+1;
    }
//    testVector(vectorB);
//    testVector(vectorX);

}



int main(){
    double *matrix = (double*) malloc(sizeof(double)*MATRIX_SIZE*MATRIX_SIZE);
    double *vectorX = (double*) malloc(sizeof(double)*MATRIX_SIZE);
    double *vectorB = (double*) malloc(sizeof(double )*MATRIX_SIZE);
    omp_set_num_threads(8);
    //initMatrix(matrix, vectorX, vectorB);
    initMatrix(matrix, vectorX, vectorB);
    algorithm(matrix,vectorX, vectorB);
    testVector(vectorX);

    free(matrix);
    free(vectorX);
    free(vectorB);
}



































