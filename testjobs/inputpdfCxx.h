
//----------------------------------------------------------------------
double xupv(double x) {
  double au = 5.107200;
  double pd = au * pow(x,0.8) * pow((1-x),3);
  return pd;
}
//----------------------------------------------------------------------
double xdnv(double x) {
  double ad = 3.064320;
  double pd = ad * pow(x,0.8) * pow((1-x),4);
  return pd;
}
//----------------------------------------------------------------------
double xglu(double x) {
  double ag = 1.7;
  double pd = ag * pow(x,-0.1) * pow((1-x),5);
  return pd;
}
//----------------------------------------------------------------------
double xdbar(double x) {
  double adbar = 0.1939875;
  double pd = adbar * pow(x,-0.1) * pow((1-x),6);
  return pd;
}
//----------------------------------------------------------------------
double xubar(double x) {
  double pd = xdbar(x) * (1-x);
  return pd;
}
//----------------------------------------------------------------------
double xsbar(double x) {
  double pd = 0.2 * ( xdbar(x) + xubar(x) );
  return pd;
}
//----------------------------------------------------------------------
double func(int* ipdf, double* x) {
  int i = *ipdf;
  double xb = *x;
  double f = 0;
  if(i ==  0) f = xglu(xb);
  if(i ==  1) f = xdnv(xb);
  if(i ==  2) f = xupv(xb);
  if(i ==  3) f = 0;
  if(i ==  4) f = xdbar(xb);
  if(i ==  5) f = xubar(xb);
  if(i ==  6) f = xsbar(xb);
  if(i ==  7) f = 0;
  if(i ==  8) f = 0;
  if(i ==  9) f = 0;
  if(i == 10) f = 0;
  if(i == 11) f = 0;
  if(i == 12) f = 0;
  return f;
}
