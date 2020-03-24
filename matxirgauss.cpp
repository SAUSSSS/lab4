#include <iostream>
#include <cmath>
#include <float.h>
void show(double* a,int n);
void vector_discrepancy(double** a, double* x,double* f, int n);
double norma(double** a,int n);
double condition_matrx(double** A, int N);
double *gauss(double** a, double* f, int n);
double  *PVR(double** a, double* f, int n);


void show(double* a,int n){
  using namespace std;
  for(int i(0); i < n; i++){
      cout <<  a[i] <<  endl;
    }
}

void vector_discrepancy(double** a, double* x,double* f, int n){
  double* temp = new double [n];
  double* r = new double[n];
  temp[0] = 0;
  temp[0] = 1 * x[0];
  for(int i = 1; i < n; i++){
      temp[i] = 0;
      temp[i] += (1 * x[i-1]) + (-2 * x[i]) + (1 * x[i+1]);
    }
  double  otv = 0;
    otv += 1 * x[0] + 1 * x[n];
   for(int i(1); i < n ; i++){
      otv += (2 * x[i]) ;
   }
   temp[n-1] = otv;

   for(int i(0); i < n ; i++){
      r[i] = f[i] - temp[i];
   }
      double gr = 0;
   for(int i(0); i < n ; i++){
      gr += r[i];
    }
     std::cout << gr << std::endl;
  delete [] temp;
}


double norma(double** a,int n){
    double max = 0,summa = 0;
    for(int j = 0; j < n; j++){
      for(int i = 0; i < n; i++){
        summa += a[i][j];
      }
      if(max < summa) max = summa;
      summa = 0;
    }
    return max;
}


double condition_matrx(double** A, int N){
  double temp;
  double** UnitMatrix = new double* [N];
  double** Inverse_Matrix = new double* [N];
  for (int i = 0; i < N; i++){
		UnitMatrix[i] = new double[N];
    Inverse_Matrix[i] = new double[N];
  }
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      Inverse_Matrix[i][j] = A[i][j];
    //Create unit matrix
        UnitMatrix[i][j] = 0.0;
          if (i == j){
            UnitMatrix[i][j] = 1.0;
          }
        }
      }
    //-------------
    for (int k = 0; k < N; k++) {
        temp = Inverse_Matrix[k][k];
        for (int j = 0; j < N; j++){
            Inverse_Matrix[k][j] /= temp;
            UnitMatrix[k][j] /= temp;
        }
        for (int i = k + 1; i < N; i++){
            temp =  Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -=  UnitMatrix[k][j] * temp;
            }
        }
    }
    for (int k = N - 1; k > 0; k--){
        for (int i = k - 1; i >= 0; i--){
            temp = Inverse_Matrix[i][k];
            for (int j = 0; j < N; j++){
                Inverse_Matrix[i][j] -= Inverse_Matrix[k][j] * temp;
                UnitMatrix[i][j] -= UnitMatrix[k][j] * temp;
            }
        }
    }
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			  Inverse_Matrix[i][j] =  UnitMatrix[i][j];

  double max = norma(Inverse_Matrix,N);

	for (int i = 0; i < N; i++){
		delete[]  UnitMatrix[i];
    delete[]  Inverse_Matrix[i];
  }
	delete[]  UnitMatrix;
  delete[]  Inverse_Matrix;
  return max;
}

int iteration = 0;
double  *PVR(double** a, double* f, int n){
  using namespace std;
  double eps = 0.0001;
  double *xn = new double[n];
  double xt[n] = {0};
  double xx[n] = {0};
  double nort = 0;
	double s = 0,ss = 0;
  double w = 1.3; //coef relaxation
  for(int i = 0; i < n; i++)
{
xn[i] = 0;
xx[i] = xn[i];
}
do
{
   nort = 0;
   for(int i = 0;i < n;i++){
      xx[i]=f[i];
      for(int j=0;j<n;j++){
          if(i!=j)
              xx[i]=xx[i]-a[i][j]*xx[j];
      }
      xx[i] /= a[i][i];

     xx[i] = w * xx[i] + (1 - w) * xn[i];

      if(fabs(xx[i] - xn[i]) > nort){
          nort = fabs(xx[i] - xn[i]);
          iteration++;
        }
      xn[i] = xx[i];
   }
} while(nort > eps);

  return xn;
}



double *gauss(double** a, double* y, int n)
{
  using namespace std;
  double *x,max;
  int k, index;
  const double eps = DBL_EPSILON;
  x = new double[n];
  k = 0;
  max = 0;
  while (k < n )
  {
    // Search for a string with maximum a[i][k]
    max = abs(a[k][k]);
    index = k;
    for (int i = k + 1; i < n; i++)
    {
      if (abs(a[i][k]) > max)
      {
        max = abs(a[i][k]);
        index = i;
      }
    }

    if (max < eps)
    {
      cout << "No solution";
    }
      // change rows
    for (int j = 0; j < n; j++)
    {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;

    // normalization of the equation
    for (int i = k; i < n; i++)
    {
      double temp = a[i][k];
      if (abs(temp) < eps) continue; // coff = 0 -> continue
      for (int j = 0; j < n ; j++)
        a[i][j] = a[i][j] / temp;
      y[i] = y[i] / temp;
      if (i == k)  continue; // equation ne vichitat samogo iz seba
      for (int j = 0; j < n ; j++)
        a[i][j] = a[i][j] - a[k][j];
      y[i] = y[i] - y[k];
    }
    k++;
  }
  // reverse substitution
  for (k = n - 1; k >= 0; k--)
  {
    x[k] = y[k];
    for (int i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
  }




int main(void){
  using namespace std;
   const int n(20);
   double** a = new double* [n+1];
   double *xx = new double [n+1];
   double *x = new double [n+1];
   double* f = new double [n+1];
   double* copyf = new double[n+1];

   for(int i(0); i < n; i++){
     a[i] = new double [n+1];
   }


//Create matrix main
   a[0][0] = 1;
   a[0][1] = 0;
   a[0][2] = 0;
   f[0] = 1;
   int v = 1;
   for(int i(1); i < n ; i++){
      a[i][v-1] = 1;
      a[i][v] = -2;
      a[i][v+1] = 1;
      f[i] = 2. / pow(i,2);
      v += 1;
   }
   int p = 2;
   a[n-1][0] = 1;
   a[n-1][n-1] = 1;
   f[n-1] = -n / 3.;
   for(int i(1); i < n ; i++){
      a[19][i] = 2;
   }
   //===============
   for(int i = 0; i < n; i++){
     copyf[i] = f[i];
   }

  double min = condition_matrx(a,n);
  double max = norma(a,n);
  double norma = max * min ;
  cout << "===Conditionality of the matrix===" << endl;
  cout << "Lambda min = " << min  << endl
  << "Lambda max = " <<  max << endl
   << "Condition number = " <<  norma << endl;
  cout << endl;

    xx = PVR(a,f,n); //Method ROS
    x = gauss(a,f,n); //Method Gauss

  cout << "======Method Gauss=====" << endl;
    show(x,n);
  cout << "===Norma Vectora===" << endl;
    vector_discrepancy(a,x,copyf,n);
    cout << endl;
  cout << "======Method ROS=====" << endl;
  cout << "w = 1.2 " << " Iteration: " << iteration <<endl;
    show(xx,n);
  cout << "===Norma Vectora===" << endl;
    vector_discrepancy(a,xx,copyf,n);



   for(int i(0); i < n; i++)
      delete a[i];

    delete [] a;
    delete [] x;
    delete [] xx;
    delete [] copyf;
    delete [] f;
   return 0;
 }
