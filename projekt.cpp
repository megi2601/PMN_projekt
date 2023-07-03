#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <map>

using namespace Eigen;

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

typedef std::complex<double> cmp;
typedef double(*ptr_f)(int, double);



double const_V(int i, double V0){
    return V0;
}

double linear_V(int i, double V0){
    return V0*i;
}

//tworzenie wektora gaussowskiego 
VectorXcd gaussian_vector(int size, double sigma, int n){
    VectorXcd v(size);
    for(int k=0; k<size; k++){
        v(k) = exp(-(k-n) * (k-n) / 2 / sigma / sigma);
    }
    return v.normalized();
}

//tworzenie wektora bazy standardowej 
VectorXcd basis_vector(int size, int n){
    VectorXcd v = VectorXcd::Zero(size);
    v(n) = 1;
    return v;
}

// tworzenie macierzy hamiltonianiu
MatrixXd create_H(int size, double V, ptr_f potential){
    MatrixXd H(size, size);
    for(int n=0; n<size; n++){
        H(n, n) = potential(n, V);
        H((n+1)%size, n) = H((n+size-1)%size, n) = -1;
    }
    return H;
}

//rozkład wektora początkowego w bazie własnej hamiltonianiu
MatrixXcd eigenbasis(int size, VectorXcd v0, EigenSolver<MatrixXd> es){
    MatrixXcd v(size, size);
    for(int k=0; k<size; k++){
        es.eigenvectors().col(k).normalize();
        cmp coeff = es.eigenvectors().col(k).dot(v0);
        v.col(k) = coeff * es.eigenvectors().col(k);
    }
    return v;
}

//tworzenie macierzy ewolucji czasowej
MatrixXd evolve_probability(MatrixXcd&v, int size, double dt, int steps, EigenSolver<MatrixXd> es){
    MatrixXd P(steps, size);
    MatrixXcd d = v;
    for(int s=0; s<steps; s++){
        for(int k=0; k<size; k++){
        double Ek = es.eigenvalues()(k).real();
        cmp z(0, -dt*Ek); //
        v.col(k) *= exp(z);
        }
        for(int l=0; l<size; l++){
        cmp v_l = v.row(l).sum();
        P(steps-1-s, l) = std::abs(v_l) * std::abs(v_l);
        }
    }
    return P;
}

void write_to_csv(std::string filename, MatrixXd matrix)
{
    std::ofstream file(filename);
    file << matrix.format(CSVFormat);
}

std::string filename(int argc, char* argv[]){
    std::string pot=argv[1], 
    V=argv[2], N = argv[3], dt = argv[4], 
    s = argv[5], vec = argv[6], n = argv[7], sigma=""; 
    if(*argv[6]=='g'){std::string s=argv[8]; sigma ="_" + s;}
    std::string filename = "PM_V"+ pot+V +"_N"+ N +"_dt" + dt +"_s" + s + "_v" + vec + n + sigma;
    filename += ".txt";
    return filename;
}


int main(int argc, char* argv[]){
    if(argc < 8) throw std::invalid_argument("brak parametrow wywowalania");
    // kolejnosc: pot, V, N, dt, steps , vec, n, sigma
    std::map<char, ptr_f> pot_func{{'l', linear_V}, {'c', const_V}};
    char pot = *argv[1];
    char vec = *argv[6];
    double V = std::stod(argv[2]);
    int N = std::stoi(argv[3]);
    double dt = std::stod(argv[4]);
    int steps = std::stoi(argv[5]);
    int n = std::stoi(argv[7]);
    VectorXcd v0;
    if(vec == 'g'){
        double sigma = std::stod(argv[8]);
        v0 = gaussian_vector(N, sigma, n);
    }
    else if(vec == 'b') v0 = basis_vector(N, n); 
    MatrixXd H = create_H(N, V, pot_func.at(pot));
    EigenSolver<MatrixXd> es(H);
    MatrixXcd m = eigenbasis(N, v0, es);
    MatrixXd p = evolve_probability(m, N, dt, steps, es);
    write_to_csv(filename(argc, argv), p);
    write_to_csv("macierz.txt", p);
    return 0;
}
