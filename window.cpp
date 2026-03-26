
#include <QPainter>
#include <stdio.h>

#include "window.h"
#include "newton.h"
#include "spline.h"

#define DEFAULT_A -1
#define DEFAULT_B 1
#define DEFAULT_N 5
#define DEFAULT_K 0
#define DEFAULT_MODE 0
#define DEFAULT_SCALE 0

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
  draw_mode = DEFAULT_MODE;
  mode_name = "mode I";
  scale = DEFAULT_SCALE;

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

void Window::change_mode()
{
  draw_mode = (draw_mode + 1)%4;
  set_mode();
}

void Window::set_mode()
{
  switch (draw_mode)
    {
      case 0:
        mode_name = "mode I";
        break;
      case 1:
        mode_name = "mode II";
        break;
      case 2:
        mode_name = "mode III";
        break;
      case 3:
        mode_name = "mode IV";
        break;
      
    }
    update ();
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

void Window::draw_graph(QPainter &painter, int width,int n,double a, double b, double (*func)(double, int,double* ,double* ,double*),double *x, double *coeff,double *tmp)
{
  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y;

  double hx = (b - a)/width;

  max_y = min_y = 0;


  for (int i = 0; i <= width; i++)
    {
      x1 = a + i * hx;
      y1 = func(x1,n,x,coeff,tmp);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  printf("\nmax{|Fmin|,|Fmax|} = %lf\n",std::max(min_y,max_y));

  x1 = a;
  y1 = func(x1,n,x,coeff,tmp);

  for (int i = 1; i <= width; i++)
    {
      x2 = a + i * hx;
      y2 = func(x2,n,x,coeff,tmp);
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
}

void Window::draw_error(QPainter &painter, int width,int n,double a, double b, double (*func)(double,int, double*, double*, double*),double (*f)(double),double *x, double *coeff,double *tmp)
{

  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y;

  double hx = (b - a)/width;

  max_y = min_y = 0;


  for (int i = 0; i <= width; i++)
    {
      x1 = a + i * hx;
      y1 = fabs(func(x1,n,x,coeff,tmp)-f(x1));
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  

  x1 = a;
  y1 = fabs(func(x1,n,x,coeff,tmp)-f(x1));

  for (int i = 1; i <= width; i++)
    {
      x2 = a + i * hx;
      y2 = fabs(func(x2,n,x,coeff,tmp)-f(x2));
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }

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
  if(draw_mode != 3)
  {x1 = a;
  y1 = f (x1);
  for (i = 1; i <= W; i++)
    {
      x2 = a + i * hx;
      y2 = f (x2);
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
  }
  // draw axis
  painter.setPen (pen_red);
  painter.drawLine (L2G(a, 0), L2G(b, 0));
  painter.drawLine (L2G(0, min_y), L2G(0, max_y));

  

  // render function name
  painter.setPen ("blue");
  painter.drawText (0, 20, f_name);
  painter.drawText (0, 40, mode_name);

  // printf("IN parintEvent\n");
  double *x,*y,*coeff,*tmp;

  //alloc memory

  x = new double[n];
  y = new double[n];
  coeff = new double[4*n];
  tmp = new double[7*(n + 1)];

  double *newton_tmp = tmp;
  double *newton_cff = coeff;

  double *spline_tmp = newton_tmp + n;
  double *spline_cff = newton_cff + n;
  
  
  make_xy(n,a,b,x,y,f);

  //Newton
  
  if(n <= 50)
    newton_coeff(n,x,y,newton_cff,newton_tmp);
  

  //Spline

  spline_coeff(n,x,y,spline_cff,spline_tmp);

  // double sval = spline_in_point(point,n,x,spline_cff,spline_tmp);

  // printf("\n sVal in point 1 = %lf\n",sval);

  ksi_append(n,x,spline_tmp);


  if(draw_mode == 0 && n <= 50) 
    draw_graph(painter,W,n,a,b,&newton_in_point,x,newton_cff,newton_tmp);
  else if(draw_mode == 1)
    {
      painter.setPen ("orange");
      draw_graph(painter,W,n,a,b,&spline_in_point,x,spline_cff,spline_tmp);
    }
  else if(draw_mode == 2)
    {
      if(n <= 50) 
        draw_graph(painter,W,n,a,b,&newton_in_point,x,newton_cff,newton_tmp);

      painter.setPen ("orange");
      draw_graph(painter,W,n,a,b,&spline_in_point,x,spline_cff,spline_tmp);
    }
  else if(draw_mode == 3)
    {  
      if(n <= 50) 
        draw_error(painter,W,n,a,b,&newton_in_point,f,x,spline_cff,spline_tmp);

      painter.setPen ("orange");
      draw_error(painter,W,n,a,b,&spline_in_point,f,x,spline_cff,spline_tmp);
    }

  //delete memory

  delete []x;
  delete []y;
  delete []coeff;
  delete []tmp;


}


void Window::resize_mult ()
{
  n *= 2;
  update();
}


void Window::resize_dev ()
{
  n /= 2;
  n = (n == 0 ? 1:n);
  update();
}


void Window::rescale_mult ()
{
  scale = 0;
  scale++;
  change_ab ();
  update ();
}


void Window::rescale_dev ()
{
  scale = 0;
  scale--;
  change_ab ();
  update ();
}

void Window::change_ab ()
{
  double mid = (a + b) / 2;

  double h = (b - a) * std::pow(2,scale);

  a = mid - h / 2;
  b = mid + h / 2;
}