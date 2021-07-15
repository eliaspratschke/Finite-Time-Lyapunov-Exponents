#include <iostream>
#include "Eigen/Dense"
#include <math.h>
#include <stdio.h>
#include <thread>
#include <vector>
#include <chrono>
//#define PI 7.0
#define PI 3.14159265

int dim = 100;
double tfinal = 20;
double sigma = 2 * PI / dim;
Eigen::Vector2d initial;
Eigen::Vector2d result;
Eigen::Matrix2d D;
Eigen::Matrix2d C;
Eigen::MatrixXd xones(dim + 2, dim + 2);
Eigen::MatrixXd xtwos(dim + 2, dim + 2);
Eigen::MatrixXd ftle(dim, dim);

//system of odes to be integrated
Eigen::Vector2d map(const double& xone, const double& xtwo, const double& time){
    Eigen::Vector2d map;
    double a = -0.5;
    double k = 0.1;
    map[0] = xtwo * xone;
    map[1] = -xtwo + a * xone * xone;
    return map;

}

//rk4 implementation
Eigen::Vector2d Advect(Eigen::Vector2d & initial, double tfinal){
    Eigen::Vector2d val;
    val[0] = initial[0];
    val[1] = initial[1];
    double dt = 0.01;
    double time = 0;
    Eigen::Vector2d kone, ktwo, kthree, kfour;
    for(int i = 0; time < tfinal; ++i){

        kone = map(val[0], val[1], time) * dt;
        //std::cout<<kone[0];
        ktwo = map(val[0] + 0.5 * kone[0], val[1] + 0.5 * kone[1], time + 0.5 * dt) * dt;
        //std::cout<<ktwo[0];
        kthree = map(val[0] + 0.5 * ktwo[0], val[1] + 0.5 * ktwo[1], time + 0.5 * dt) * dt;
        //std::cout<<kthree[0];
        kfour = map(val[0] + kthree[0], val[1] + kthree[1], time + dt) * dt;
        //std::cout<<kfour[0];
        val[0] += 1.0 / 6.0 * (kone[0] + 2.0 * ktwo[0] + 2.0 * kthree[0] + kfour[0]);
        val[1] += 1.0 / 6.0 * (kone[1] + 2.0 * ktwo[1] + 2.0 * kthree[1] + kfour[1]);
        time += dt;
    }
    return val;
}

//integrate a batch of initial conditions
void calculate(int begin, int end) {
    if (end > dim +2){end = dim + 2;}
    for (int i = begin; i < end; ++i) {
        //std::cout << "beep" << i;
        for (int j = 0; j < dim + 1; ++j) {
            initial[0] = -1.0 * PI - sigma + i * sigma;
            initial[1] = -1.0 * PI - sigma + j * sigma;
            //std::cout<< initial[0] << " " << initial[1] << "\n";
            result = Advect(initial, tfinal);
            xones(i, j) = result[0];
            xtwos(i, j) = result[1];
            //result = Advect(initial, tfinal);
        }

    }
}

int main() {
    std::cout<<"beep";
    //assign multithreading
    auto start = std::chrono::high_resolution_clock::now();
    int numthreads = 1;
    std::vector<std::thread> threads;
    int chunksize = (dim + 2) / (numthreads);
    int index = 0;

    while(index * chunksize < dim + 2){
        std::cout<<"thread: "<<index <<"goes from"<< index * chunksize<< "to" << (index + 1) * (chunksize)<< "\n";
        std::thread worker(calculate, index * chunksize, (index + 1) * chunksize);
        index += 1;
        threads.push_back(std::move(worker));
    }

    for (auto& td : threads) td.join();
    FILE * file;
    file = fopen("ftle.txt","w");
    //print to file
    for(int i = 1; i < dim + 1; ++i){
        for(int j = 1; j < dim + 1; ++j) {
        D(0,0) = (xones(i - 1, j) - xones(i + 1, j)) / 2 / sigma;
        D(0,1) = (xones(i, j - 1) - xones(i, j + 1)) / 2 / sigma;
        D(1,0) = (xtwos(i - 1, j) - xtwos(i + 1, j)) / 2 / sigma;
        D(1,1) = (xtwos(i, j - 1) - xtwos(i, j + 1)) / 2 / sigma;
        C = D.transpose() * D;
        ftle(j - 1, i - 1) = 1 / tfinal / 2 * log(C.operatorNorm());
        if (ftle(j - 1, i -1) > 0.5){std::cout<<"beeeeeeeeep";}
        if (j == dim){fprintf(file, "%2.10f \n", 1 / tfinal / 2 * log(C.operatorNorm()));}
        else{ fprintf(file, "%2.10f ", 1 / tfinal / 2 * log(C.operatorNorm()));}
        }
    }
    fclose (file);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout<<"execution time with "<<numthreads<<" threads:" << duration.count()<< "seconds";
    return 0;
}
