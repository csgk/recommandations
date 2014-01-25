#include <iostream>
#include <string>
#include </usr/local/include/Eigen/Eigen>
#include </usr/local/include/Eigen/SVD>
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
    *result = MatrixXd::Zero(rows, columns);
    
    float rate;
    
	string buffer;
    getline(file, buffer);
    int n = 0;
    while(!file.eof() && n<number_of_entries){
        stringstream line(buffer);
		line >> i >> j >> rate;
        (*result)(i-1, j-1) = rate;
        getline(file, buffer);
        n++;
    }

    //for (int i = 0; i < 100; ++i) {
    //    for (int j = 0; j < 100; ++j) {
      //      cout << (*data)(i,j) << " ";
        //}
        //cout << "\\" << endl;
    //}

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
	return (sqrt(error/n));

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
        if(0==card){p[i] = 0;}
        else{
            p[i] = (p[i]/card) - mean;
            cout << "p[" << i << "] == " << p[i] << "\t";
        }
    }
    
    cout << endl;
    
    cout << "columns == " << columns << endl;
    for (int j = 0; j<columns; ++j) {
        card = 0;
        for (int i = 0; i<rows; ++i) {
            o[j] += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
        if(0==card){o[j] = 0;}
        else{
            o[j] = (o[j]/card) - mean;
        }
        cout << "card == " << card <<"\t";
        cout << "o[" << j << "] == " << o[j] << endl;
    }
    
    for (int j = 0; j<columns; ++j) {
        for (int i = 0; i<rows; ++i) {
            cout << "i == " << i << " j == " << j << endl;
            (*data)(i, j) = mean + p[i] + o[j];
        }
    }
}

/* ********************************** 2  2 ********************************** */

/* ********************************** 3  1 ********************************** */

/* ********************************** 3  3 ********************************** */

void neighbourhood_minimisation_algorithm (MatrixXd *data) {
    cout << "Comme ça va être beau ! *_*" << endl;
    MatrixXd S = (*data) * data->transpose();
    cout << "calcul de S fini" << endl;
    int n = (int)S.cols();
    vector<float> norm(n);
    vector<float> sum(n);
    for (int i = 0; i < n; ++i) {
        norm[i] = data->col(i).norm();
        if (0 == norm[i]) {
            cout << "Bouhouhouhouhouhouhou " << i << endl;
        }
    }
    cout << "calcul de norme fini" << endl;
    for (int i = 0; i < n; ++i) {
        sum[i] = 0;
        for (int k = 0; k < n; ++k) {
            S(i,k) /= norm[i]*norm[k];
            sum[i] += abs(S(i,k));
        }
    }
    cout << "calcul de somme fini" << endl;
    *data = S * (*data);
    for (int i = 0; i < n ; ++i) {
        for (int j = 0; j < n; ++j) {
            (*data)(i,j) /= sum[i];
        }
    }
    cout << "calcul de Ap fini" << endl;
    return;
}

/* ********************************** 4    ********************************** */

void singular_values_algorithm(MatrixXd *data, int model_number) {
    JacobiSVD<MatrixXd> svd(*data, ComputeFullU | ComputeFullV);
    cout << "j'ai calculé la SVD ! #swag" << endl;
    cout << "je calcule ma nouvelle matrice... #suspens" << endl;
    cout << "je calcule U ! #fuckU" << endl;
    MatrixXd U = svd.matrixU();
    cout << "ma taille est " << U.rows() << " " << U.cols() << endl;
    cout << "je calcule V^T ! #vitessegrandV" << endl;
    MatrixXd VT = svd.matrixV().transpose();
    cout << "ma taille est " << VT.rows() << " " << VT.cols() << endl;
    cout << "je fais du produit de matrices !" << endl;
    MatrixXd S = MatrixXd::Zero(U.cols(), VT.rows());
    for (int i = 0; i<model_number; ++i) {
        S(i,i)=svd.singularValues()[i];
        cout << S(i,i) << " ";
    }
    *data =  U * S * VT;
    cout << "J'ai ma nouvelle matrice, trop bien ! #YOLO" << endl;
}

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
    
    cout << "j'ai initialisé la matrice, youpi ! #joie" << endl;
    cout << "Je fais avec l'algorithme moyen #pipeau" << endl;

    average_algorithm(data_matrix, rows, columns);

    cout << ":-) Wouhou, on trouve : " << rmse(source_file, data_matrix, rows, columns) << endl;
    cout << "Je recommence #mammamia" << endl;
    
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);

    neighbourhood_minimisation_algorithm(data_matrix);
    
    cout << ":-) Wouhou, on trouve : " << rmse(source_file, data_matrix, rows, columns) << endl;

    cout << "et on recommence ! #lavieestbelle";
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    
    singular_values_algorithm(data_matrix, 50);
    cout << ":-) Wouhou, on trouve : " << rmse(source_file, data_matrix, rows, columns) << endl;
}
