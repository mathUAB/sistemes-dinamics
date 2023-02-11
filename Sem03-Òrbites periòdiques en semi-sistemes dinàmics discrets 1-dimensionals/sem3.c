#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 1e-10
#define TOL 1e-8
#define TOL2 1e-7
#define N 500
#define Ndrop 8000
#define muN 4000
#define PI 3.14159265358979323846
#define FEIGENBAUM 4.669201609102990671
// #define RICHARDSON

// REMEMBER TO CHANGE THE INTERVAL OF PLOTTING IN THE SCRIPT plot.sh

// function 1
#define muStart 2.9
#define muEnd 3.999
#define function function1();
#define f(x, mu) ((mu) * (x) * (1 - (x)))

// function 2
// #define muStart 4
// #define muEnd 5.828427124 - 0.001
// #define function function2();
// #define f(x, mu) ((mu) * (x) * (1 - (x)) / (1 + (x)))

// function 3
// #define muStart 1.9
// #define muEnd 2.598076212 - 0.001
// #define function function3();
// #define f(x, mu) ((mu) * (x) * (1 - (x) * (x)))

// // function 4
// #define muStart 0.35
// #define muEnd 1
// #define function function4()
// #define f(x, mu) ((mu)*cosl((x)*PI))

// // function 5
// #define muStart 0.9
// #define muEnd 2
// #define function function5();
// #define f(x, mu) (((x) <= 0.5) ? ((mu) * (x)) : ((mu) * (1 - (x))))

// // function 6
// #define muStart 19
// #define muEnd 100
// #define function function6();
// #define f(x, mu) ((mu) * (x)*exp(-(mu) * (x)-1))

#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))  // maximum of two values
#define MIN(X, Y) (((X) > (Y)) ? (Y) : (X))  // minimum of two values

void function1(void);
void function2(void);
void function3(void);
void function4(void);
void function5(void);
void function6(void);
long double fn(long double x, long double mu, int n);
long double fn_x(long double x, long double mu, int n, long double h);
long double fn_mu(long double x, long double mu, int n, long double h);
long double fn_xx(long double x, long double mu, int n, long double h);
long double fn_xmu(long double x, long double mu, int n, long double h);
long double richardson(long double x, long double mu, int n, long double g(long double, long double, int, long double));
long double Fx(long double x, long double mu, int n);
long double Fy_power2(long double x, long double mu, int n);
long double Fy(long double x, long double mu, int n);
long double* newton(long double Fx(long double, long double, int), long double Fy(long double, long double, int), long double x0, long double mu0, int n);
long double bisection(int n, long double mu0);
long double* orbit(long double x0, long double mu);
long double** orbits(long double x0);

int main(int argc, char const* argv[]) {  // function 3
  if (argc != 2) {
    fprintf(stdout, "You have to introduce a value.\n");
    return 1;
  }
  long double x0;
  sscanf(argv[1], "%Lf", &x0);
  long double** data = orbits(x0);

  FILE* fp = fopen("data.dat", "w");
  if (fp == NULL) {
    fprintf(stdout, "Problems opening the file.\n");
    return 1;
  }
  for (int i = 0; i < muN; i++) {
    for (int j = 2; j < data[i][1]; j++) {
      fprintf(fp, "%Lf %Lf\n", data[i][0], data[i][j]);
    }
  }
  fclose(fp);

  function;
  // long double mu = 20;
  // for (long double x = 5; x < 20; x += 0.5)
  //   printf("f(%Lf, %Lf) = %g\n", x, mu, f(x, mu));

  return 0;
}

long double min(long double* array, int n) {
  long double min = array[0];
  for (int i = 1; i < n; i++) {
    if (array[i] < min)
      min = array[i];
  }
  return min;
}

void function1(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 11;
  long double mu0, mu[k_MAX], lambda[k_MAX];
  for (int k = 0; k < k_MAX; k++) {
    if (k > 1)
      mu0 = mu[k - 1];
    else if (k >= 0)
      mu0 = 3;
    v = newton(Fx, Fy_power2, 0.5, mu0, (int)pow(2, k));
    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %Lf\t \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  // for (int k = 0; k < k_MAX; k++) {
  //   if (k == 0)
  //     mu[k] = bisection((int)pow(2, k + 1), 2.9);
  //   else
  //     mu[k] = bisection((int)pow(2, k + 1), mu[k - 1]);

  //   if (k > 1)
  //     lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
  //   else
  //     lambda[k] = NAN;
  //   printf("\u00B5_%i = %Lf\t \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  // }

  v = newton(Fx, Fy, 0.25, 3.7, 7);
  printf("\nPERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.22, 3.73, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.15, 3.82, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n\n", v[1]);
}

void function2(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 4;
  long double mu0 = 4.2, mu[k_MAX], lambda[k_MAX];
  long double x0 = 0.5, eps = 0.8;
  for (int k = 0; k < k_MAX; k++) {
    if (k > 1) {
      mu0 = (mu[k - 1] - mu[k - 2]) / FEIGENBAUM + mu[k - 1];
      // printf("(x0, mu0) = (%Lf, %Lf)\n", x0, mu0);
    } else if (k > 0)
      mu0 = mu[k - 1] + eps;
    v = newton(Fx, Fy_power2, x0, mu0, (int)pow(2, k));

    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %Lf, \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  v = newton(Fx, Fy, 0.17, 5.62, 7);
  printf("PERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.15, 5.5, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.1, 5.62, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n", v[1]);
}

void function3(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 5;
  long double mu0 = 2, mu[k_MAX], lambda[k_MAX];
  long double x[k_MAX];
  long double x0 = 0.5;
  long double eps = 0.2;
  for (int k = 0; k < k_MAX; k++) {
    if (k > 1) {
      x0 = x[k - 1];
      mu0 = (mu[k - 1] - mu[k - 2]) / FEIGENBAUM + mu[k - 1];
    } else if (k > 0) {
      x0 = x[k - 1] + eps;
      mu0 = mu[k - 1] + eps;
    }
    // printf("x0=%Lf", x0);
    v = newton(Fx, Fy_power2, x0, mu0, (int)pow(2, k));
    x[k] = v[0];
    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %Lf, \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  v = newton(Fx, Fy, 0.4, 2.37, 7);
  printf("PERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.35, 2.39, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.32, 2.42, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n", v[1]);
}

void function4(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 10;
  long double mu0, mu[k_MAX], lambda[k_MAX];
  long double x0;
  for (int k = 0; k < k_MAX; k++) {
    if (k > 0) {
      x0 = 0;
      mu0 = mu[k - 1];
    } else {
      x0 = 0.4;
      mu0 = 0.35;
    }
    v = newton(Fx, Fy_power2, x0, mu0, (int)pow(2, k));
    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %Lf, \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  v = newton(Fx, Fy, -0.35, 0.675, 7);
  printf("PERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, -0.39, 0.68, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, -0.5, 0.73, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n", v[1]);
}

void function5(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 9;
  long double mu0 = 1, mu[k_MAX], lambda[k_MAX];
  long double x0 = 0.5;
  for (int k = 0; k < k_MAX; k++) {
    if (k > 1) {
      x0 = 0.4;
      mu0 = mu[k - 1];
    } else if (k > 0) {
      x0 = 0.45;
      mu0 = 1;
    }
    v = newton(Fx, Fy_power2, x0, mu0, (int)pow(2, k));
    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %Lf, \u03BB_%i = %Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  v = newton(Fx, Fy, 0.392, 1.46, 7);
  printf("PERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.35, 1.52, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.3, 1.62, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n", v[1]);
}

void function6(void) {
  printf("PERIODS POWERS OF 2:\n");
  long double* v;
  int k_MAX = 6;
  long double mu0, mu[k_MAX], lambda[k_MAX];
  long double x0;

  for (int k = 0; k < k_MAX; k++) {
    if (k == 2)
      v = newton(Fx, Fy_power2, 0.1, 40, (int)pow(2, k));
    else if (k == 3)
      v = newton(Fx, Fy_power2, 0.1, 39, (int)pow(2, k));
    else if (k == 4)
      v = newton(Fx, Fy_power2, 0.1008, 40, (int)pow(2, k));
    // else if (k == 5)
    //   v = newton(Fx, Fy_power2, 0.1, mu[4], (int)pow(2, k));
    // else if (k == 6)
    //   v = newton(Fx, Fy_power2, 0.550040, mu0, (int)pow(2, k));
    // else if (k == 7)
    //   v = newton(Fx, Fy_power2, 0.561521, mu0, (int)pow(2, k));
    // else if (k == 8)
    //   v = newton(Fx, Fy_power2, 0.508800, mu0, (int)pow(2, k));
    else if (k >= 5) {
      // for (int i = 0; i < 1000; i++) {
      //   v = newton(Fx, Fy_power2, x0 + i * step, 40, (int)pow(2, k));
      //   if (v[1] > 39 && v[1] < 40.1)
      //     printf("k=%i, x0 = %Lf, \u00B5_%i = %Lf\n", k + 1, x0 + i * step, k + 1, v[1]);
      // }
      v = newton(Fx, Fy_power2, x0, mu[k - 1], (int)pow(2, k));
    } else {
      x0 = 0.1;
      mu0 = 20;
      v = newton(Fx, Fy_power2, x0, mu0, (int)pow(2, k));
    }

    mu[k] = v[1];
    if (k > 1)
      lambda[k] = (mu[k - 1] - mu[k - 2]) / (mu[k] - mu[k - 1]);
    else
      lambda[k] = NAN;
    printf("\u00B5_%i = %.10Lf, \u03BB_%i = %.10Lf\n", k + 1, mu[k], k + 1, lambda[k]);
  }

  v = newton(Fx, Fy, 0.004, 47.5, 7);
  printf("PERIOD 7: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.003, 50, 5);
  printf("PERIOD 5: \u00B5 = %Lf\n", v[1]);
  v = newton(Fx, Fy, 0.018, 60, 3);
  printf("PERIOD 3: \u00B5 = %Lf\n", v[1]);
}

// f^n(x, mu)
long double fn(long double x, long double mu, int n) {
  long double result = x;
  for (int i = 0; i < n; i++)
    result = f(result, mu);
  return result;
}

// numeric differentiation using symmetric difference formula (order 1)
long double fn_x(long double x, long double mu, int n, long double h) {
  return (fn(x + h, mu, n) - fn(x - h, mu, n)) / (2 * h);
  // return (fn(x + h, mu, n) - fn(x, mu, n)) / (h);
}

// numeric differentiation using symmetric difference formula (order 1)
long double fn_mu(long double x, long double mu, int n, long double h) {
  return (fn(x, mu + h, n) - fn(x, mu - h, n)) / (2 * h);
  // return (fn(x, mu + h, n) - fn(x, mu, n)) / (h);
}

// numeric differentiation using symmetric difference formula (order 2)
long double fn_xx(long double x, long double mu, int n, long double h) {
  return (fn(x + h, mu, n) - 2 * fn(x, mu, n) + fn(x - h, mu, n)) / (h * h);
  // return (fn_x(x + h, mu, n, h) - fn_x(x, mu, n, h)) / h;
}

// numeric differentiation for mixed derivatives
long double fn_xmu(long double x, long double mu, int n, long double h) {
  return (fn(x + h, mu + h, n) - fn(x - h, mu + h, n) - fn(x + h, mu - h, n) + fn(x - h, mu - h, n)) / (4 * h * h);
  // return (fn_x(x, mu + h, n, h) - fn_x(x, mu, n, h)) / (h);
}

// richardson extrapolation
long double richardson(long double x, long double mu, int n, long double g(long double, long double, int, long double)) {
  int steps = 3;
  long double matrix[steps][steps];
  long double q = 2;
  long double h = 0.0001;
  long double qPowN;

  // first row
  for (int r = 0; r < steps; r++) {
    matrix[0][r] = g(x, mu, n, h / powl(q, 2 * r));  // the derivative has error term: ? * h^2 + ? * h^4 + ? * h^6 + ...
    // printf("(%i, %i) = %.10Lf, %Lf\n", 0, r, matrix[0][r], h / powl(q, 2 * r));
  }
  // richardson approximation for the other rows.
  for (int m = 1; m < steps; m++) {
    for (int r = 0; r < steps - m; r++) {
      qPowN = powl(q, 2 * m);  // the derivative has error term: ? * h^2 + ? * h^4 + ? * h^6 + ...
      matrix[m][r] = (qPowN * matrix[m - 1][r + 1] - matrix[m - 1][r]) / (qPowN - 1);
      // printf("(%i, %i) = %.10Lf\n", m, r, matrix[m][r]);
    }
  }
  // we take the only element of the last row, that has error term ? * h^(2 * steps)
  return matrix[steps - 1][0];
}

// f(x, mu) = x
long double Fx(long double x, long double mu, int n) {
  return fn(x, mu, n) - x;
}

// diff(f(x, mu), x) = -1
long double Fy_power2(long double x, long double mu, int n) {
#if defined(RICHARDSON)
  return richardson(x, mu, n, fn_x) + 1;
#else
  return fn_x(x, mu, n, TOL) + 1;
#endif  // RICHARDSON
}

// diff(f(x, mu), x) = 1
long double Fy(long double x, long double mu, int n) {
#if defined(RICHARDSON)
  return richardson(x, mu, n, fn_x) - 1;
#else
  return fn_x(x, mu, n, TOL) - 1;
#endif  // RICHARDSON
}

long double bisection(int n, long double mu0) {  // remember that when an orbit borns it is the unique attracting periodic orbit.
  int iter = 10000;
  long double step = 0.00001;
  long double x = 0.5;
  long double mu = mu0, mu_i, mu_f;
  while (mu < muEnd && fabsl(fn(x, mu, n) - x) > EPS) {
    mu += step;
    x = fn(x, mu, iter);
    // printf("hola: %Lf, %Lf\n", x, fn(x, mu, n));
  }
  mu_i = mu - step;
  mu_f = mu;
  while (fabsl(mu_f - mu_i) > EPS) {
    mu = (mu_i + mu_f) / 2;
    x = fn(0.5, mu, iter);
    if (fabsl(fn(x, mu, n) - x) < EPS)
      mu_f = mu;
    else
      mu_i = mu;
  }
  return (mu_i + mu_f) / 2;
}

long double* newton(long double Fx(long double, long double, int), long double Fy(long double, long double, int), long double x0, long double mu0, int n) {
  const int MAX_ITER = 100;
  int iter = 0;
  long double x = x0, mu = mu0;
  long double x_aux, mu_aux;
  long double d, deval;
  long double d11, d12, d21, d22;  // differential matrix
  long double det;
  // functions of newton: F = (Fx, Fy).
  do {
#if defined(RICHARDSON)

    d11 = richardson(x, mu, n, fn_x);
    d12 = richardson(x, mu, n, fn_mu);
    d21 = richardson(x, mu, n, fn_xx);
    d22 = richardson(x, mu, n, fn_xmu);
#else
    d11 = fn_x(x, mu, n, TOL);
    d12 = fn_mu(x, mu, n, TOL);
    d21 = fn_xx(x, mu, n, TOL);
    d22 = fn_xmu(x, mu, n, TOL);
#endif  // RICHARDSON
    det = d11 * d22 - d12 * d21;
    x_aux = x - (d22 * Fx(x, mu, n) - d12 * Fy(x, mu, n)) / det;
    mu_aux = mu - (-d21 * Fx(x, mu, n) + d11 * Fy(x, mu, n)) / det;
    d = MAX(fabsl(x - x_aux), fabsl(mu - mu_aux));
    x = x_aux;
    mu = mu_aux;
    deval = MAX(fabsl(Fx(x, mu, n)), fabsl(Fy(x, mu, n)));
    iter++;
  } while (iter < MAX_ITER && MAX(d, deval) > EPS);
  long double* v = (long double*)malloc(2 * sizeof(long double));
  v[0] = x;
  v[1] = mu;
  return v;
}

long double* orbit(long double x0, long double mu) {
  long double* array = (long double*)malloc((N + 3) * sizeof(long double));
  array[0] = mu;
  // printf("nova mu: %Lf", mu);
  long double last = x0;
  int n = 0, count = 0;
  for (int i = 0; i < Ndrop; i++)
    last = f(last, mu);
  for (int i = 2; i < N + 2; i++) {
    last = f(last, mu);
    if (i > 2 && (last - array[2] < TOL && -last + array[2] < TOL)) {  // they start to repeat.
      n = i - 2;
      if (n / (count + 1) < 20 && n < 50)
        count++;
      else
        break;
    }
    array[i] = last;
  }
  if (n == 0)
    array[1] = N;
  else
    array[1] = n;
  // first position of the array is the value of mu
  // second position of the array is the number of points to be plotted with this mu
  // the other positions are the points
  return array;
}

long double** orbits(long double x0) {
  long double muStep = (muEnd - muStart) * 1. / muN;
  long double** matrix = (long double**)malloc(muN * sizeof(long double*));
  for (int mu = 0; mu <= muN; mu++)
    matrix[mu] = orbit(x0, muStart + mu * muStep);
  return matrix;
}
