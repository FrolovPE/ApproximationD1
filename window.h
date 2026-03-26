#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  const char *f_name;
  double a;
  double b;
  int n;
  int k;
  double (*f) (double);
  int draw_mode;
  const char *mode_name;
  int scale;

public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
  void draw_graph(QPainter &painter, int width,int n,double a, double b, double (*func)(double,int, double*, double*, double*),double *x, double *coeff,double *tmp);
  void draw_error(QPainter &painter, int width,int n,double a, double b, double (*func)(double,int, double*, double*, double*),double (*f)(double),double *x, double *coeff,double *tmp);
  void change_ab ();
public slots:
  void change_func ();
  void set_func ();
  void change_mode ();
  void set_mode();
  void resize_mult ();
  void resize_dev ();
  void rescale_mult ();
  void rescale_dev ();
  //void resize ();

protected:
  void paintEvent (QPaintEvent *event);
};



#endif
