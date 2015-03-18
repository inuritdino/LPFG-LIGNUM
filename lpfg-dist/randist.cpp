float ran_pos(float range)
{
  float r;
  do{
    r = ran(range);
  }
  while(r == 0);
  
  return r;
}
//Unifrom with limits
float ran_uniform(const float min, const float max)
{
  return (min + ran_pos(min - max));
}

//Polar (Box-Mueller) method, from GSL, from Knuth v2, 3rd ed, p122.
float ran_gauss(const float mu, const float sigma)
{//Generate Gaussian variate with mean MU and std SIGMA
  float x, y, r2, v = -1.0;
  while(v < 0){//Simple workaround for negative outputs
    do{
      x = -1 + 2 * ran_pos(1.0);
      y = -1 + 2 * ran_pos(1.0);
      r2 = x * x + y * y;
    }
    while(r2 > 1.0 || r2 == 0);
  
    v = mu + sigma * y * sqrt(-2.0 * log(r2) / r2);
  }
  
  return v;
}

float ran_gauss_any(const float mu, const float sigma)
{//Generate Gaussian variate with mean MU and std SIGMA
  float x, y, r2, v = -1.0;
  do{
    x = -1 + 2 * ran_pos(1.0);
    y = -1 + 2 * ran_pos(1.0);
    r2 = x * x + y * y;
  }
  while(r2 > 1.0 || r2 == 0);
  
  v = mu + sigma * y * sqrt(-2.0 * log(r2) / r2);
  
  return v;
}
