#include <iostream>
#include <string>
#include </usr/local/include/Eigen/Eigen>
#include <vector>
#include <fstream>

using namespace std;
using namespace Eigen;

/* ************************************************************************** */
/* *****                     initialisation/erreur                      ***** */
/* ************************************************************************** */

MatrixXd *init_matrix(char *source_file, int rows, int columns, double number_of_entries) {
    std::fstream file;
	file.open(source_file, fstream::in);
    
	if (!file.good()) {
        std::cout<<"File '"<<source_file<<"' not found"<<endl;
		return NULL;
	}
    
    int i, j;
    MatrixXd *result = new MatrixXd(rows, columns);
    for (i = 0; i<rows; ++i) {
        for (j = 0; j<columns; ++j) {
            (*result)(i, j) = 0;
        }
    }
    
    float rate;
    
	string buffer;
    getline(file, buffer);
    int n = 0;
    while(!file.eof() && n<number_of_entries){
        stringstream line(buffer);
		line >> i >> j >> rate;
        (*result)(i, j) = rate;
        getline(file, buffer);
        n++;
    }
    
	file.close();
	return result;
    
}

/* ************************************************************************** */

float rmse(char *source_file, MatrixXd *data, int rows, int columns) {
    std::fstream file;
	file.open(source_file, fstream::in);

	if (!file.good()) {
        std::cout<<"File '"<<source_file<<"' not found"<<endl;
	}
    
    int i, j;
    float rate, error = 0;
    int n = 0;
    
	string buffer;
    getline(file, buffer);
    
    while(!file.eof()){
        stringstream line(buffer);
		line >> i >> j >> rate;
        error += pow((*data)(i, j) - rate, 2);
        getline(file, buffer);
        n++;
    }
    
	file.close();
	return (sqrt(error)/n);

}

/* ************************************************************************** */
/* *****                           ALGORITHMES                          ***** */
/* ************************************************************************** */

/* ********************************** 2  1 ********************************** */

void average_algorithm(MatrixXd *data, int rows, int columns) {
    int card = 0;
    float mean = 0;
    vector<float> p(rows, 0), o(columns, 0);
    
    for (int i = 0; i<rows; ++i) {
        for (int j = 0; j<columns; ++j) {
            mean += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
    }
    mean /= (float)card;

    for (int i = 0; i<rows; ++i) {
        card = 0;
        for (int j = 0; j<columns; ++j) {
            p[i] += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
        p[i] = (p[i]/card) - mean;
    }
    
    for (int j = 0; j<columns; ++j) {
        card = 0;
        for (int i = 0; i<rows; ++i) {
            o[j] += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
        o[j] = (o[j]/card) - mean;
    }
    
    for (int j = 0; j<columns; ++j) {
        for (int i = 0; i<rows; ++i) {
            (*data)(i, j) = mean + p[i] + o[j];
        }
    }
}

/* ********************************** 2  2 ********************************** */

/* ********************************** 3  1 ********************************** */

/* ********************************** 3  3 ********************************** */

/* ********************************** 4    ********************************** */



/* ************************************************************************** */
/* *****                             MAIN                               ***** */
/* ************************************************************************** */

int main(int argc, char** argv) {
    char *source_file;
    int rows, columns, number_of_entries;
    
    if (argc < 2) {
        std::cout << "No imput argument" <<std::endl;
        exit (0);
    }
    
    for (int i=1; i<argc; ++i) {
        if (0 == strcmp(argv[i], "-input-file")) {
            source_file = argv[++i];
        }
        else if (0 == strcmp(argv[i], "-persons")) {
            rows = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-objects")) {
            columns = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-entries")) {
            number_of_entries = atoi(argv[++i]);
        }
        
    }
    
    MatrixXd *data_matrix = init_matrix(source_file, rows, columns, number_of_entries);

    average_algorithm(data_matrix, rows, columns);
    cout << ":-) Wouhou, on trouve : " << rmse(source_file, data_matrix, rows, columns) << endl;
}
