#include <iostream>
#include <fstream>
#include <array>
#include <random>
#include <armadillo>





using namespace std;
using namespace arma;


void MonteCarlo(vec &players, int MCSteps, int N, int transactions, double lambda, double alpha, double gamma, ofstream &outFile, vec &binCounts, double binSize, double m0, ofstream &outFileErr);
void outPut(vec &players, int MCSteps, int N, int transactions, mat &expectVal);
double findVariance(vec &players, int transaction, int N, double m0);
void makeBins(vec &players, vec &binCount, double binSize);


int main(int argc, char *argv[])
{
    if (argc < 8){
        cout << "To few arguments given. Expected number of persons, Monte Carlo cycles, transactions, start money, lambda, alpha and gamma" << endl;
    }
    int N = stoi(argv[1]);
    int MCSteps = stoi(argv[2]);
    int transactions = stoi(argv[3]);


    double startMoney = stod(argv[4]);
    double lambda = stod(argv[5]);
    double alpha = stod(argv[6]);
    double gamma = stod(argv[7]);

    double binSize = 0.01*startMoney;


    ofstream outFileVar = ofstream("variance.txt");
    ofstream outFileErr = ofstream("distError.txt");

    //outFileVar.open("variance.txt");
    ofstream outFileParameter;
    outFileParameter.open("parameters.txt");

    ofstream binParameters = ofstream("binParameters.txt");



    outFileParameter << "N " << N << "\n";
    outFileParameter << "MCSteps " << MCSteps << "\n";
    outFileParameter << "Trasactions " << transactions << "\n";
    outFileParameter << "StartingMoney " << startMoney << "\n";
    outFileParameter << "Lambda " << lambda << "\n";
    outFileParameter << "Alpha " << alpha << "\n";
    outFileParameter << "Gamma " << gamma;


    double binEnd;
    if (alpha > 0 || gamma > 0){
        binEnd = 2*startMoney/(sqrt(lambda + 0.1)) + startMoney;
    }
    else{
        binEnd = 2*startMoney/(sqrt(lambda + 0.1)) + startMoney;
    }

    int binNum = int(binEnd/double(binSize));
    cout << binNum << endl;
    vec binCounts = zeros(binNum);
    vec players = ones(N)*startMoney;

    cout << alpha << " " << gamma << endl;

    binParameters << MCSteps << " " << N << " " << startMoney << " " << binSize << " " << binNum << " " << (binEnd) << " " << lambda << " " << alpha << " " << gamma << endl;
    binParameters.close();





    MonteCarlo(players, MCSteps,  N,transactions,lambda,alpha,gamma,outFileVar,binCounts,binSize,startMoney,outFileErr);
    //outPut(players, MCSteps, N, transactions, expectVal);





    binCounts.save("bins.bin",raw_binary);
    players.save("data.bin",raw_binary);

    outFileErr.close();
    outFileVar.close();

    //players.save("data.bin",raw_binary);

    cout << "Finished" << endl;



}


void MonteCarlo(vec &players, int MCSteps, int N, int transactions, double lambda,double alpha,double gamma, ofstream &outFile, vec &binCounts,double binSize,double m0,ofstream &outFileErr){


    random_device rd;
    mt19937_64 gen(rd());

    uniform_real_distribution<double> distribution(0.0,N);
    uniform_real_distribution<double> eps(0.0,1.0);

    int writingFreq = 100;
    double p = 0;

    mat c = zeros(N,N);
    double maxTransactions = 1;



    for (int i = 0; i < MCSteps; i++){


        players.fill(m0);


        for (int j = 0; j < transactions; j++){
            int index_i = distribution(gen);
            int index_j = distribution(gen);

            double epsFac = eps(gen);


            if (players(index_i) - players(index_j) == 0){
                p = 1.;
            }
            else{
                p = 2*pow(fabs((players(index_i) - players(index_j))/double(m0)),-alpha)*(pow((c(index_i,index_j)+1)/(maxTransactions+1),gamma));
            }

            if (eps(gen) < p && (index_i != index_j)){





                double m1 = lambda*players(index_i) + (1-lambda)*epsFac*    (players(index_i) + players(index_j));
                double m2 = lambda*players(index_j) + (1-lambda)*(1-epsFac)*(players(index_i) + players(index_j));

                //cout << "hei" << endl;

                players(index_i) = m1;
                players(index_j) = m2;


                c(index_j,index_i) += 1;
                c(index_i,index_j) += 1;

                if (c(index_j,index_i) > maxTransactions){
                    maxTransactions = c(index_j,index_i);
                }

                else if (c(index_i,index_j) > maxTransactions){
                    maxTransactions = c(index_i,index_j);
                }

            }





            if (MCSteps == 1){
                if (j%writingFreq == 0){

                    double mean = 0;
                    for (int i = 0; i < N; i++){
                        mean += players(i)/m0;

                    }

                    mean /= (N);


                    outFile << (j+1) << " " << findVariance(players,j+1,N,m0) << " " << mean << "\n";
                }

            }
        }

    vec tempCounts = binCounts;
    makeBins(players,binCounts,binSize);

    if (i > 1){
        outFileErr << i+1 << " " << norm(tempCounts/double(i-1) - binCounts/(double(i)) ) << endl;
    }





    }

}

//void outPut(vec &players, int MCSteps, int N, int transactions, mat &expectVal){

//    vec means = zeros(N);
//    for (int k = 0; k < N;k++){
//        means(k) = expectVal(0,k) / MCSteps;
//    }

//    means.save("data.bin",raw_binary);
//}

double findVariance(vec &players, int transaction, int N,double m0){

    double mean = 0;
    double secondMoment = 0;

    for (int i = 0; i < N; i++){
        mean += players(i)/m0;
        secondMoment += players(i)/m0*players(i)/m0;
    }

    mean /= (N);
    //cout << secondMoment << endl;
    secondMoment /= (N);

    return (secondMoment - mean*mean);

}


void makeBins(vec &players, vec &binCount, double binSize){

    for (int i = 0; i < players.size();i++){
        for (int j = 0; j < binCount.size();j++){
            if(players(i)> (j-1)*binSize && players(i)< (j)*binSize){
                binCount(j) += 1;
            }
        }
    }

    //binCount.print();


}




