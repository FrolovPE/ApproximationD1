
#include <QPainter>
#include <stdio.h>

#include "window.h"
#include "newton.h"

#define DEFAULT_A -1
#define DEFAULT_B 1
#define DEFAULT_N 5
#define DEFAULT_K 0
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))

static
double f_0 (double x)
{
  x=x;
  return 1;
}

static
double f_1 (double x)
{
  return x;
}

static
double f_2 (double x)
{
  return x * x;
}
static
double f_3 (double x)
{
  return x * x * x;
}
static
double f_4 (double x)
{
  return x * x * x * x;
}
static
double f_5 (double x)
{
  return std::exp(x);
}
static
double f_6 (double x)
{
  return 1.0/(25 * x * x + 1.0);
}

static void make_xy(int n, double a, double b, double *x, double *y, double (*func)(double))
{
    double h = (b-a)/n;

    for(int i = 0 ; i < n; i++)
    {
        x[i] = a + i*h;
        y[i] = func(x[i]);
    }
}



Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;
  k = DEFAULT_K;

  func_id = k;

  set_func();

  // change_func ();
  
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || (argc > 3 && (sscanf (argv[3], "%d", &n) != 1))
      || (argc > 4 && (sscanf (argv[4], "%d", &k) != 1))
      || n <= 0
      ||(k < 0 || k > 6)
      )
    return -2;
  // printf("A = %lf B = %lf N = %d K = %d\n",a,b,n,k);

  func_id = k;
  set_func();

  // printf("A = %lf B = %lf N = %d K = %d func_id = %d\n",a,b,n,k,func_id);

  return 0;
}

/// change current function for drawing
void Window::change_func ()
{
  // printf("func_id = %d\n",func_id);
  func_id = (func_id + 1) % 7;

  set_func();
}

void Window::set_func ()
{
 
  switch (func_id)
    {
      case 0:
        f_name = "k = 0 f (x) = 1";
        f = f_0;
        break;
      case 1:
        f_name = "k = 1 f (x) = x";
        f = f_1;
        break;
      case 2:
      f_name = "k = 2 f (x) = x * x";
      f = f_2;
      break;
      case 3:
        f_name = "k = 3 f (x) = x * x * x";
        f = f_3;
        break;
      case 4:
        f_name = "k = 4 f (x) = x * x * x * x";
        f = f_4;
        break;
      case 5:
      f_name = "k = 5 f (x) = exp(x)";
      f = f_5;
      break;
      case 6:
        f_name = "k = 6 f (x) = 1/(25 * x * x + 1)";
        f = f_6;
        break;
    }
  update ();
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max)
{
  double x_gl = (x_loc - a) / (b - a) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{  
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y;
  // Line for graph: green, line width 2, solid
  QPen pen_green(Qt::green, 2, Qt::SolidLine);
  // Line for axes: red, line width 1, solid
  QPen pen_red(Qt::red, 1, Qt::SolidLine);
  // Window width
  int W = width (), i;
  // Step
  double hx = (b - a) / W;

  painter.setPen (pen_green);

  // calculate min and max for current function
  max_y = min_y = 0;
  for (i = 0; i <= W; i++)
    {
      x1 = a + i * hx;
      y1 = f (x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  // draw approximated line for graph
  x1 = a;
  y1 = f (x1);
  for (i = 1; i <= W; i++)
    {
      x2 = a + i * hx;
      y2 = f (x2);
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
  // draw axis
  painter.setPen (pen_red);
  painter.drawLine (L2G(a, 0), L2G(b, 0));
  painter.drawLine (L2G(0, min_y), L2G(0, max_y));

  // render function name
  painter.setPen ("blue");
  painter.drawText (0, 20, f_name);

  // printf("IN parintEvent\n");
  double *x,*y,*coeff,*tmp;

  //alloc memory

  x = new double[n];
  y = new double[n];
  coeff = new double[n];
  tmp = new double[n];
  
  make_xy(n,a,b,x,y,f);
  newton_coeff(n,x,y,coeff,tmp);


  printf("\nNewton coeffs k = %d\n",func_id);
  
  for(int i = 0; i < n; i++)
  {
    printf(" %lf",coeff[i]);
  }

  double val = newton_in_point(1,n,x,coeff);

  printf("\n Val in point 1 = %lf\n",val);


  //delete memory

  delete []x;
  delete []y;
  delete []coeff;
  delete []tmp;


}
