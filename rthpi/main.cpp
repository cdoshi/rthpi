/* Function to develop magnetic dipole fitting algorithm in C++ */


// Include relevant header files
#include <iostream>
#include <math.h>
#include <string>
#include "eigen/Eigen/Dense"
#include <vector>
#include <algorithm>
#include <fstream>      // std::ifstream


/*********************************************************************************/
// Declare all structures to be used

struct coilParam {
    Eigen::MatrixXd pos;
    Eigen::MatrixXd mom;
};

struct dipError {
    double error;
    Eigen::MatrixXd moment;
};

struct sens {
    Eigen::MatrixXd coilpos;
    Eigen::MatrixXd coilori;
    Eigen::MatrixXd tra;
};

/*********************************************************************************/
static std::vector <double>*base_arr;

/*********************************************************************************/
// Function declarations

dipError dipfitError (Eigen::MatrixXd, Eigen::MatrixXd, struct sens);
Eigen::MatrixXd ft_compute_leadfield(Eigen::MatrixXd, struct sens);
Eigen::MatrixXd magnetic_dipole(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);
coilParam dipfit(struct coilParam, struct sens, Eigen::MatrixXd);
coilParam fminsearch(Eigen::MatrixXd,int, int, int, Eigen::MatrixXd, struct sens);
bool compar (int, int);
Eigen::MatrixXd pinv(Eigen::MatrixXd);

/*********************************************************************************/

/*********************************************************************************/
// MAIN PROGRAM STARTS HERE

int main()
{
    Eigen::MatrixXd data;
    struct sens sensors;
    struct coilParam coil;
    float a, b, c;
    int num_sensors = 270, i;
    std::string line;

    /*********************************************************************************/
    // Read the sensor locations
    sensors.coilpos = Eigen::MatrixXd::Zero(num_sensors,3);
    i=0;
    std::ifstream myfile ("C:\\Users\\cdoshi\\Desktop\\HPI_min\\FieldTrip_DipoleFit_BabyMEG\\simulation_seok\\sensorspos.txt");
    if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
        sscanf_s (line.c_str(),"%f %f %f",&a, &b, &c);
        sensors.coilpos(i,0) = a;sensors.coilpos(i,1) = b;sensors.coilpos(i,2) = c;
        i++;
      }
      myfile.close();
    }
    myfile.close();
    /*********************************************************************************/
    // Read the sensor orientations
    sensors.coilori = Eigen::MatrixXd::Zero(num_sensors,3);
    i=0;
    myfile.open("C:\\Users\\cdoshi\\Desktop\\HPI_min\\FieldTrip_DipoleFit_BabyMEG\\simulation_seok\\sensorsori.txt");
    if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
        sscanf_s (line.c_str(),"%f %f %f",&a, &b, &c);
        sensors.coilori(i,0) = a;sensors.coilori(i,1) = b;sensors.coilori(i,2) = c;
        i++;
      }
      myfile.close();
    }
    myfile.close();

    /*********************************************************************************/
    // sensors tra
    sensors.tra = Eigen::MatrixXd::Identity(num_sensors,num_sensors);

    /*********************************************************************************/
    // Read the data for 1st dipole
    data = Eigen::MatrixXd::Zero(num_sensors,1);
    i=0;
    myfile.open("C:\\Users\\cdoshi\\Desktop\\HPI_min\\FieldTrip_DipoleFit_BabyMEG\\simulation_seok\\data_dipole1.txt");
    if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
        sscanf_s (line.c_str(),"%f",&a);
        data(i,0) = a;
        i++;
      }
      myfile.close();
    }
    myfile.close();

    /*********************************************************************************/
    coil.pos = Eigen::MatrixXd::Zero(1,3);
    coil.mom = Eigen::MatrixXd::Zero(1,3);

    coil = dipfit(coil, sensors, data);

    std::cout << coil.pos << std::endl;

    return 0;
}


/*********************************************************************************
 * dipfit function is adapted from Fieldtrip Software. It has been
 * heavity edited for use with MNE-X Software
 *********************************************************************************/

coilParam dipfit(struct coilParam coil, struct sens sensors, Eigen::MatrixXd data)
{
    // Initialize variables
    int display = 0;
    int maxiter = 100;
    dipError temp;

    coil = fminsearch(coil.pos, maxiter, 2 * maxiter * coil.pos.cols(), display, data, sensors);

    temp = dipfitError(coil.pos, data, sensors);

    coil.mom = temp.moment;

    return coil;
}

/*********************************************************************************
 * fminsearch Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
 * X = fminsearch(X0, maxiter, maxfun, display, data, sensors) starts at X0 and
 * attempts to find a local minimizer
 *********************************************************************************/

coilParam fminsearch(Eigen::MatrixXd pos,int maxiter, int maxfun, int display, Eigen::MatrixXd data, struct sens sensors)
{
    double tolx, tolf, rho, chi, psi, sigma, func_evals, usual_delta, zero_term_delta, temp1, temp2;
    std::string header, how;
    int n, itercount, prnt;
    Eigen::MatrixXd onesn, two2np1, one2n, v, y, v1, tempX1, tempX2, xbar, xr, x, xe, xc, xcc, xin;
    std::vector <double> fv, fv1;
    std::vector <int> idx;
    coilParam coil;
    dipError tempdip, fxr, fxe, fxc, fxcc;

    tolx = tolf = 1e-4;

    switch(display)
    {
        case 0:
            prnt = 0;
            break;
        default:
            prnt = 1;
    }

    header = " Iteration   Func-count     min f(x) Procedure";

    n = pos.cols();

    // Initialize parameters
    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    onesn = Eigen::MatrixXd::Ones(1,n);
    two2np1 = one2n = Eigen::MatrixXd::Zero(1,n);

    for(int i = 0;i < n;i++) {
        two2np1(i) = 1 + i;
        one2n(i) = i;
    }

    v = v1 = Eigen::MatrixXd::Zero(n, n+1);
    fv.resize(n+1);idx.resize(n+1);fv1.resize(n+1);

    for(int i = 0;i < n; i++) v(i,0) = pos(i);

    tempdip = dipfitError(pos, data, sensors);
    fv[0] = tempdip.error;


    func_evals = 1;itercount = 0;how = "";

    // Continue setting up the initial simplex.
    // Following improvement suggested by L.Pfeffer at Stanford
    usual_delta = 0.05;             // 5 percent deltas for non-zero terms
    zero_term_delta = 0.00025;      // Even smaller delta for zero elements of x
    xin = pos.transpose();

    for(int j = 0;j < n;j++) {
        y = xin;

        if(y(j) != 0) y(j) = (1 + usual_delta) * y(j);
        else y(j) = zero_term_delta;

        v.col(j+1).array() = y;
        pos = y.transpose();
        tempdip = dipfitError(pos, data, sensors);
        fv[j+1] = tempdip.error;
    }

    // Sort elements of fv
    base_arr = &fv;
    for (int i = 0; i < n+1; i++) idx[i] = i;

    sort (idx.begin(), idx.end(), compar);

    for (int i = 0;i < n+1;i++) {
        v1.col(i) = v.col(idx[i]);
        fv1[i] = fv[idx[i]];
    }

    v = v1;fv = fv1;

    how = "initial simplex";
    itercount = itercount + 1;
    func_evals = n + 1;

    tempX1 = Eigen::MatrixXd::Zero(1,n);

    while ((func_evals < maxfun) && (itercount < maxiter)) {
        for (int i = 0;i < n;i++) tempX1(i) = abs(fv[0] - fv[i+1]);
        temp1 = tempX1.maxCoeff();

        tempX2 = Eigen::MatrixXd::Zero(n,n);

        for(int i = 0;i < n;i++) tempX2.col(i) = v.col(i+1) -  v.col(0);

        tempX2 = tempX2.array().abs();

        temp2 = tempX2.maxCoeff();

        if((temp1 <= tolf) && (temp2 <= tolx)) break;

        xbar = v.block(0,0,n,n).rowwise().sum();
        xbar /= n;

        xr = (1+rho) * xbar - rho * v.block(0,n,v.rows(),1);

        x = xr.transpose();
        std::cout << "Iteration Count: " << itercount << ":" << x << std::endl;

        fxr = dipfitError(x, data, sensors);

        func_evals = func_evals+1;

        if (fxr.error < fv[0]) {
            // Calculate the expansion point
            xe = (1 + rho * chi) * xbar - rho * chi * v.col(v.cols()-1);
            x = xe.transpose();
            fxe = dipfitError(x, data, sensors);
            func_evals = func_evals+1;

            if(fxe.error < fxr.error) {
                v.col(v.cols()-1) = xe;
                fv[n] = fxe.error;
                how = "expand";
            }
            else {
                v.col(v.cols()-1) = xr;
                fv[n] = fxr.error;
                how = "reflect";
            }
        }
        else {
            if(fxr.error < fv[n-1]) {
                v.col(v.cols()-1) = xr;
                fv[n] = fxr.error;
                how = "reflect";
            }
            else { // fxr.error >= fv[:,n-1]
                // Perform contraction
                if(fxr.error < fv[n]) {
                    // Perform an outside contraction
                    xc = (1 + psi * rho) * xbar - psi * rho * v.col(v.cols()-1);
                    x = xc.transpose();
                    fxc = dipfitError(x, data, sensors);
                    func_evals = func_evals + 1;

                    if(fxc.error <= fxr.error) {
                        v.col(v.cols()-1) = xc;
                        fv[n] = fxc.error;
                        how = "contract outside";
                    }
                    else {
                        // perform a shrink
                        how = "shrink";
                    }
                }
                else {
                    xcc = (1 - psi) * xbar + psi * v.col(v.cols()-1);
                    x = xcc.transpose();
                    fxcc = dipfitError(x, data, sensors);
                    func_evals = func_evals+1;
                    if(fxcc.error < fv[n]) {
                        v.col(v.cols()-1) = xcc;
                        fv[n] = fxcc.error;
                        how = "contract inside";
                    }
                    else {
                        // perform a shrink
                        how = "shrink";
                    }
                }
                if(how.compare("shrink") == 0) {
                    for(int j = 1;j < n+1;j++) {
                        v.col(j).array() = v.col(0).array() + sigma * (v.col(j).array() - v.col(0).array());
                        x = v.col(j).array().transpose();
                        tempdip = dipfitError(x,data, sensors);
                        fv[j] = tempdip.error;
                    }
                }
            }
        }
        // Sort elements of fv
        base_arr = &fv;
        for (int i = 0; i < n+1; i++) idx[i] = i;
        sort (idx.begin (), idx.end (), compar);
        for (int i = 0;i < n+1;i++) {
            v1.col(i) = v.col(idx[i]);
            fv1[i] = fv[idx[i]];
        }
        v = v1;fv = fv1;
        itercount = itercount + 1;
    }

    x = v.col(0).transpose();
    std::cout << x << std::endl;
    coil.pos = x;
    return coil;
}


/*********************************************************************************
 * dipfitError computes the error between measured and model data
 * and can be used for non-linear fitting of dipole position.
 * The function has been compared with matlab dipfit_error and it gives
 * same output
 *********************************************************************************/

dipError dipfitError(Eigen::MatrixXd pos, Eigen::MatrixXd data, struct sens sensors)
{
    // Variable Declaration
    struct dipError e;
    Eigen::MatrixXd lf, dif;

    // Compute lead field for a magnetic dipole in infinite vacuum
    lf = ft_compute_leadfield(pos, sensors);

    //e.moment = pinv(lf)*data;
    e.moment = pinv(lf) * data;

    dif = data - lf * e.moment;

    e.error = dif.array().square().sum()/data.array().square().sum();

    return e;
}

/*********************************************************************************
 * ft_compute_leadfield computes a forward solution for a dipole in a a volume
 * conductor model. The forward solution is expressed as the leadfield
 * matrix (Nchan*3), where each column corresponds with the potential or field
 * distributions on all sensors for one of the x,y,z-orientations of the dipole.
 * The function has been compared with matlab ft_compute_leadfield and it gives
 * same output
 *********************************************************************************/

Eigen::MatrixXd ft_compute_leadfield(Eigen::MatrixXd pos, struct sens sensors)
{

    Eigen::MatrixXd pnt, ori, lf;

    pnt = sensors.coilpos; // position of each coil
    ori = sensors.coilori; // orientation of each coil

    lf = magnetic_dipole(pos, pnt, ori);

    lf = sensors.tra * lf;

    return lf;
}

/*********************************************************************************
 * magnetic_dipole leadfield for a magnetic dipole in an infinite medium
 * The function has been compared with matlab magnetic_dipole and it gives same output
 *********************************************************************************/

Eigen::MatrixXd magnetic_dipole(Eigen::MatrixXd pos, Eigen::MatrixXd pnt, Eigen::MatrixXd ori) {

    double u0 = 1e-7;
    int nchan;
    Eigen::MatrixXd r, r2, r5, x, y, z, mx, my, mz, Tx, Ty, Tz, lf, temp;

    nchan = pnt.rows();

    // Shift the magnetometers so that the dipole is in the origin    
    pnt.array().col(0) -= pos(0);pnt.array().col(1) -= pos(1);pnt.array().col(2) -= pos(2);

    r = pnt.array().square().rowwise().sum().sqrt();

   r2 = r5 = x = y = z = mx = my = mz = Tx = Ty = Tz = lf = Eigen::MatrixXd::Zero(nchan,3);

    for(int i = 0;i < nchan;i++) {
            r2.row(i).array().fill(pow(r(i),2));
            r5.row(i).array().fill(pow(r(i),5));
    }

    for(int i = 0;i < nchan;i++) {
        x.row(i).array().fill(pnt(i,0));
        y.row(i).array().fill(pnt(i,1));
        z.row(i).array().fill(pnt(i,2));
    }
    mx.col(0).array().fill(1);my.col(1).array().fill(1);mz.col(2).array().fill(1);

    Tx = 3 * x.cwiseProduct(pnt) - mx.cwiseProduct(r2);
    Ty = 3 * y.cwiseProduct(pnt) - my.cwiseProduct(r2);
    Tz = 3 * z.cwiseProduct(pnt) - mz.cwiseProduct(r2);

    for(int i = 0;i < nchan;i++) {
        lf(i,0) = Tx.row(i).dot(ori.row(i));
        lf(i,1) = Ty.row(i).dot(ori.row(i));
        lf(i,2) = Tz.row(i).dot(ori.row(i));
    }

    for(int i = 0;i < nchan;i++) {
        for(int j = 0;j < 3;j++) {
            lf(i,j) = u0 * lf(i,j)/(4 * M_PI * r5(i,j));
        }
    }
    return lf;
}

bool compar (int a, int b){
  return ((*base_arr)[a] < (*base_arr)[b]);
}

Eigen::MatrixXd pinv(Eigen::MatrixXd a)
{
    double epsilon = std::numeric_limits<double>::epsilon();
    Eigen::JacobiSVD< Eigen::MatrixXd > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(),0).matrix().asDiagonal() * svd.matrixU().adjoint();
}
