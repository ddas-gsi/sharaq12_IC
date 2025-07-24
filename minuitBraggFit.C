

#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"

using namespace std;

// Bragg Peak Function: Exponential rise followed by Gaussian peak
// Parameters: X0, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss
double BraggPeakFunction(double *x, double *par)
{
    double X = x[0];
    double X0 = par[0];         // Start range
    double Xpeak = par[1];      // Peak position
    double Aexp = par[2];       // Exponential coefficient A
    double Bexp = par[3];       // Exponential coefficient B
    double sigmaGauss = par[4]; // Gaussian sigma
    double ampGauss = par[5];   // Gaussian amplitude

    // Exponential part: X < Xpeak
    if (X < Xpeak)
    {
        return Aexp + std::exp(Bexp * X);
    }
    // Gaussian part: X >= Xpeak
    else
    {
        return ampGauss * std::exp(-((X - Xpeak) * (X - Xpeak)) / (2 * sigmaGauss * sigmaGauss));
    }
}

void minuitBraggFit()
{
    //     vector<double> dedx = {7.99115087, 8.8731933, 8.62346948, 9.3155056, 9.10757759, 9.82433461,
    //                            9.82761015, 10.21501279, 10.865725, 10.83596941, 11.11509707, 11.70866108,
    //                            11.44794399, 12.60570605, 13.85539448, 14.68005782, 16.15993464, 15.47958455,
    //                            17.659119, 18.28436637, 19.52873129, 20.00990523, 22.19400068, 22.51017344,
    //                            27.32119477, 28.96951372, 31.4497906, 32.6843946, 21.50476169, 1.18248749};

    // vector<double> dedx = {14.7616, 14.4325, 14.7226, 15.5101, 16.5349, 17.6861,
    //                        17.2964, 16.7352, 18.8607, 19.2896, 19.5484, 23.3070,
    //                        20.1268, 25.1670, 26.6031, 31.9376, 35.8301, 34.9059,
    //                        38.9812, 25.3107, 2.8099, 0., 0., 0.,
    //                        0., 0., 0., 0., 0., 0.};

    // vector<double> dedx = {7.32071522, 7.82827624, 7.33952095, 7.19638762, 7.62800456, 8.48939165,
    //                        8.52694536, 7.5622171, 8.87589595, 9.43075272, 9.75539284, 11.67498684,
    //                        9.96555989, 12.48639826, 12.61517761, 15.25030255, 17.8332722, 17.34517357,
    //                        14.56591076, 3.17703888, 3.64640298, 4.00350538, 3.90499165, 3.6467452,
    //                        3.2925889, 3.3594713, 3.15955336, 3.30648822, 3.20312328, 3.09637315};

    // vector<double> dedx = {7.32071522, 7.82827624, 7.33952095, 7.19638762, 7.62800456, 8.48939165,
    //                        8.52694536, 7.5622171, 8.87589595, 9.43075272, 9.75539284, 11.67498684,
    //                        9.96555989, 12.48639826, 12.61517761, 15.25030255, 17.8332722, 17.34517357,
    //                        14.56591076, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // std::vector<double> dedx = {19.32729584, 21.51697259, 21.58434716, 21.82326818, 25.54554109,
    //                             26.98384971, 26.51230829, 25.89668093, 30.48157802, 33.93453495,
    //                             33.69877395, 37.88128076, 21.00758739, 25.42028626, 27.05742287,
    //                             30.71572035, 34.64311335, 34.25773564, 39.01536542, 38.34262309,
    //                             23.90808365, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    vector<double> dedx = {21.42356781, 24.6401789, 24.93408999, 25.50771337, 28.85274793,
                           32.33671948, 32.53732414, 32.26805847, 37.38257832, 38.82056129,
                           21.12100395, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // vector<double> dedx = {32.58839224, 35.49095931, 38.00736232, 39.49036173, 38.29114841, 15.34692186, 0.,
    //                        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    vector<double> Segment = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                              17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    const int n = 30;

    TCanvas *c1 = new TCanvas("Bragg Peak", "Bragg Peak", 800, 600);
    TGraph *graph = new TGraph(n, Segment.data(), dedx.data());

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);
    graph->Draw("AP"); // "A" ensures axes are drawn

    // Define the fit function
    TF1 *braggFit = new TF1("braggFit", BraggPeakFunction, 0, 30, 6); // 6 parameters

    // Initial parameter guesses
    braggFit->SetParameters(0, 5, 1, 0.1, 2, 5); // (X0, Xpeak, Aexp, Bexp, sigmaGauss, ampGauss)
    braggFit->SetParNames("X0", "Xpeak", "Aexp", "Bexp", "sigmaGauss", "ampGauss");

    // Set parameter limits
    braggFit->SetParLimits(0, 0, 10);  // X0 should be between 0 and 1000
    braggFit->SetParLimits(1, 0, 30);  // Xpeak should be between 0 and 1000
    braggFit->SetParLimits(2, 0, 40);  // Aexp should be between 0 and 100
    braggFit->SetParLimits(3, -1, 1);  // Bexp should be between -1 and 0
    braggFit->SetParLimits(4, 0, 5);   // sigmaGauss should be between 0 and 100
    braggFit->SetParLimits(5, 0, 100); // ampGauss should be between 0 and 1000

    // Fit using Minuit
    graph->Fit(braggFit, "R"); // "R" restricts the fit to the defined range

    c1->Modified();
    c1->Update();
    gPad->WaitPrimitive(); // Keeps the canvas open
}

int main(int argc, char **argv)
{
    TApplication app("ROOT Application", &argc, argv);
    minuitBraggFit();
    app.Run();
    return 0;
}
