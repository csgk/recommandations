#include <iostream>
#include <string>
#include </usr/local/include/Eigen/Eigen>
#include </usr/local/include/Eigen/SVD>
#include <vector>
#include <list>
#include <fstream>

using namespace std;
using namespace Eigen;

/* ************************************************************************** */
/* *****                     initialisation/erreur                      ***** */
/* ************************************************************************** */


MatrixXd *init_matrix(char *source_file, int rows, int columns, long number_of_entries) {
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
    int n = 1;
    while(!file.eof()){
        if (n == number_of_entries) {n = 0;}
        else {
            stringstream line(buffer);
            line >> i >> j >> rate;
            (*result)(i-1, j-1) = rate;
        }
        getline(file, buffer);
        n++;
    }

	file.close();
	return result;
    
}

/* ************************************************************************** */

void print_matrix(int rows, int columns, MatrixXd *data) {
    return;
    cout << endl << endl;
    for (int i=0; i<rows; ++i) {
        for (int j=0; j<columns; ++j) {
            cout << "\t" << ((float)((int)(10000. * (*data)(i,j))))/10000.;
        }
            cout << endl;
    }
    cout << endl;
}

/* ************************************************************************** */

float rmse(char *source_file, MatrixXd *data, int rows, int columns, int number_of_entries, bool profil) {
    std::fstream file;
	file.open(source_file, fstream::in);

	if (!file.good()) {
        std::cout<<"File '"<<source_file<<"' not found"<<endl;
	}
    
    int i, j, card = 0;
    float rate, error = 0;
    
    vector<int> object_number(rows,0);
    vector<int> person_number(columns,0);
    
    list<int> real_values_i;
    list<int> real_values_j;
    list<float> real_values_rate;
    
	string buffer;
    getline(file, buffer);
    int n = 1;
    
    while(!file.eof()){
        stringstream line(buffer);
        line >> i >> j >> rate;
        if (n == number_of_entries) {
            real_values_i.push_front(i-1);
            real_values_j.push_front(j-1);
            real_values_rate.push_front(rate);
            n=0;
            ++card;
        }
        else {
            object_number[i-1]++;
            person_number[j-1]++;
        }
        getline(file, buffer);
        n++;
    }

    file.close();

    vector<float> object_error(rows,0);
    vector<float> person_error(columns,0);
    
    vector<int> object_card(rows,0);
    vector<int> person_card(columns,0);

     while (! real_values_i.empty()) {
        i = real_values_i.front();
        j = real_values_j.front();
        rate = real_values_rate.front();
        object_card[object_number[i]]++;
        person_card[person_number[j]]++;
        object_error[object_number[i]] += pow(max(min((*data)(i, j),5.),0.) - rate, 2.);
        person_error[person_number[j]] += pow(max(min((*data)(i, j),5.),0.) - rate, 2.);
        error += pow(max(min((*data)(i, j),5.),0.) 	- rate, 2.);
        real_values_i.pop_front();
        real_values_j.pop_front();
        real_values_rate.pop_front();
    }
    
    if (profil) {

    cout << "profilage de l'erreur en fonction du nombre de notes par objet" << endl;
    for (i = 0; i < object_error.size() ; ++i) {
        cout << i << "\t" << sqrt(object_error[i]/(float)object_card[i]) << endl;
    }
    cout << endl;

    cout << "profilage de l'erreur en fonction du nombre de notes par personne" << endl;
    for (j = 0; j < person_error.size() ; ++j) {
        cout << j << "\t" << sqrt(person_error[j]/(float)person_card[j]) << endl;
    }

    cout << endl;
    }
	return (sqrt(error/card));

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
        }
    }
    
    
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
    }
    for (int j = 0; j<columns; ++j) {
        for (int i = 0; i<rows; ++i) {
            (*data)(i, j) = mean + p[i] + o[j];
        }
    }
}

/* ********************************** 2  2 ********************************** */

void convex_minimisation_algorithm(MatrixXd *data, int rows, int columns, int iter_number, float lambda, float mu) {
    int card = 0;
    float mean = 0;
    vector<float> p(rows, 0), o(columns, 0);
    vector<float> sum_a_rows(rows, 0), sum_a_cols(columns, 0);
    
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
            sum_a_rows[i] += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
        if(0==card){sum_a_rows[i] = 0;}
        else{
            p[i] = (sum_a_rows[i]/card) - mean;
        }
    }
    
    for (int j = 0; j<columns; ++j) {
        card = 0;
        for (int i = 0; i<rows; ++i) {
            sum_a_cols[j] += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
        if(0==card){sum_a_cols[j] = 0;}
        else{
            o[j] = (sum_a_cols[j]/card) - mean;

        }
    }
    float erreur;
    
    erreur = 0;
    for (int i = 0; i<rows; ++i) {for (int j = 0; j<columns; ++j) {
        erreur += pow((*data)(i,j) - mean + p[i] + o[j], 2.);
    }}
    //cout << "0 : err = " << erreur << endl;
        
    float sum_o, sum_p;
    for (int k=0; k<iter_number; ++k) {
        sum_o = 0;
        for (int j = 0; j < columns; ++j) {sum_o += o[j];}
        sum_p = 0;
        for (int i = 0; i < rows; ++i) {sum_p += o[i];}

        for (int i = 0; i<rows; ++i) {
            p[i] -= mu * (columns*(p[i]-mean) + sum_o + sum_a_rows[i]);
            if (p[i] >= mu*lambda) {p[i] -= mu*lambda;}
            else if (p[i] <= -mu*lambda) {p[i] += mu*lambda;}
            else {p[i] = 0;}
        }

        for (int j = 0; j<columns; ++j) {
            o[j] -= mu * (rows*(o[j]-mean) + sum_p + sum_a_cols[j]);
            if (o[j] >= mu*lambda) {o[j] -= mu*lambda;}
            else if (o[j] <= -mu*lambda) { o[j] += mu*lambda;}
            else {o[j] = 0;}
        }
        erreur = 0;
        for (int i = 0; i<rows; ++i) {for (int j = 0; j<columns; ++j) {
            erreur += pow((*data)(i,j) - mean + p[i] + o[j], 2.);
        }}
       // cout << k << " : err = " << erreur << endl;

    }
    for (int j = 0; j<columns; ++j) {
        for (int i = 0; i<rows; ++i) {
            (*data)(i, j) = mean + p[i] + o[j];
        }
    }
}

/* ********************************** 3  1 ********************************** */

void tri (vector<int> &rang, int i, MatrixXd &S, int dim) { // faire un meilleur tri
    int transfert;
    for (int t=0; t<dim; ++t) {
        for (int x=0; x<dim-1; ++x) {
            if ( S(rang[x], i) < S(rang[x+1], i) ) {
                transfert = rang[x];
                rang[x] = rang[x+1];
                rang[x+1] = transfert;
            }
        }
    }
}




MatrixXd* neighbourhood_k_algorithm (MatrixXd* data, int k) {
    MatrixXd S = (*data) * (data->transpose());
    int n = (int)S.cols();
    int m = (int)(*data).cols();
    vector<float> norm(n);
    
    float mean=0;    int card = 0;
    for (int i = 0; i<n; ++i) {
        for (int j = 0; j<m; ++j) {
            mean += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
    }
    mean /= card;
    
    MatrixXd *result = new MatrixXd(n, m);

    for (int i = 0; i < n; ++i) {
        norm[i] = data->row(i).norm();
    }
    
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            S(i,k) /= norm[i]*norm[k];
        }
    }

    vector<int> rank(n);
    float sum, value;
    int count;
    for (int t=0; t<n; ++t) { rank[t] = t; }
    
    for (int i = 0; i < n; ++i) {
        tri(rank, i, S, n);
        for (int j = 0; j < m; ++j) {
            value = sum = 0;
            count = 0;
            for (int ii = 0; ii<n && count<k; ++ii) {
                if ((*data)(rank[ii],j) != 0) {
                    count++;
                    sum += abs( S(rank[ii],i) );
                    value += S(i,rank[ii]) * (*data)(rank[ii],j);
                }
            }
            if ( 0==count || 0==sum ) {(*result)(i,j) = mean;}
            else{(*result)(i,j) = value / sum;}
        }
    }
    return result;
    
}

/* ********************************** 3  3 ********************************** */

void neighbourhood_minimisation_algorithm (MatrixXd *data) {
    MatrixXd S = (*data) * (data->transpose());
    int n = (int)S.cols();
    int m = (int)(*data).cols();

    float mean=0;    int card = 0;
    for (int i = 0; i<n; ++i) {
        for (int j = 0; j<m; ++j) {
            mean += (*data)(i, j);
            if( 0 != (*data)(i, j) ){card++;}
        }
    }
    mean /= card;
    
    vector<float> norm(n);
    MatrixXd sum = MatrixXd::Zero(data->rows(),data->cols());
    for (int i = 0; i < n; ++i) {
        norm[i] = data->row(i).norm();
    }
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            S(i,k) /= norm[i]*norm[k];
            for (int j = 0; j < m; ++j) {
                if (0 != (*data)(k,j)) {
                    sum(i,j) += abs(S(i,k));
                }
            }
        }
    }
    MatrixXd M = S * (*data);
    for (int i = 0; i < n ; ++i) {
        for (int j = 0; j < m; ++j) {
            if (0 == sum(i,j)) {
                (*data)(i,j) = mean;
            }
            else {
                (*data)(i,j) = M(i,j) / sum(i,j);
            }
        }
    }
    return;
}

/* ********************************** 4    ********************************** */

void singular_values_algorithm(MatrixXd *data, int model_number) {
    JacobiSVD<MatrixXd> svd(*data, ComputeThinU | ComputeThinV);
    MatrixXd U = svd.matrixU();
    MatrixXd VT = svd.matrixV().transpose();
    MatrixXd S = MatrixXd::Zero(U.cols(), VT.rows());
    for (int i = 0; i<model_number; ++i) {
        S(i,i)=svd.singularValues()[i];
    //    cout << S(i,i) << " ";
    }
    *data =  U * S * VT;
}

/* ************************************************************************** */
/* *****                             MAIN                               ***** */
/* ************************************************************************** */


int main(int argc, char** argv) {
    char *source_file;
    int rows, columns;
    int number_of_entries;
    
    if (argc < 2) {
        std::cout << "No imput argument" <<std::endl;
        exit (0);
    }
    
    bool profil = false;
    int svd_n = 10, neighk = 50;
    float conv = 10000.;
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
        else if (0 == strcmp(argv[i], "-svd")) {
            svd_n = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-convex")) {
            conv = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "-neighb")) {
            neighk = atoi(argv[++i]);
        }
        else if (0 == strcmp(argv[i], "--details")) {
            profil = true;
        }
    }

    
    
    cout << "---------------------------------- ALGORITHME 1 1 ----------------------------------" << endl;
    MatrixXd *data_matrix = init_matrix(source_file, rows, columns, number_of_entries);

    average_algorithm(data_matrix, rows, columns);
    print_matrix(rows, columns, data_matrix);
    cout << "On trouve : " << rmse(source_file, data_matrix, rows, columns, number_of_entries, profil) << endl;
    
    cout << "---------------------------------- ALGORITHME 1 2 ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    
    convex_minimisation_algorithm(data_matrix, rows, columns, 100, conv, 1./(rows+columns));
    print_matrix(rows, columns, data_matrix);
    cout << "On trouve : " << rmse(source_file, data_matrix, rows, columns, number_of_entries, profil) << endl;

    
    cout << "---------------------------------- ALGORITHME 3 1 ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    data_matrix = neighbourhood_k_algorithm(data_matrix, neighk);
    print_matrix(rows, columns, data_matrix);
    cout << "On trouve : " << rmse(source_file, data_matrix, rows, columns, number_of_entries, profil) << endl;
    
    cout << "---------------------------------- ALGORITHME 3 3 ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    
    neighbourhood_minimisation_algorithm(data_matrix);
    print_matrix(rows, columns, data_matrix);
    cout << "On trouve : " << rmse(source_file, data_matrix, rows, columns, number_of_entries, profil) << endl;
    
    cout << "---------------------------------- ALGORITHME 4   ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    singular_values_algorithm(data_matrix, svd_n);
    print_matrix(rows, columns, data_matrix);
    cout << "On trouve : " << rmse(source_file, data_matrix, rows, columns, number_of_entries, profil) << endl;
    
    /*
    cout << "----------------------------------    MATRICE     ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, number_of_entries);
    print_matrix(rows, columns, data_matrix);
    
    cout << "---------------------------------- VRAIE  MATRICE ----------------------------------" << endl;
    data_matrix = init_matrix(source_file, rows, columns, rows*columns+1);
    print_matrix(rows, columns, data_matrix);
     */

}
