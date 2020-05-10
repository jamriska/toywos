#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define M_PI 3.14159265358979323846264338327950288

using namespace std;

struct Vec2f
{
  float v[2];
  
  Vec2f() {}
  Vec2f(const Vec2f& u)    { v[0] = u[0]; v[1] = u[1]; }   
  Vec2f(float v0,float v1) { v[0] = v0;   v[1] = v1;   }

  float&       operator[](int i)       { return v[i]; }
  const float& operator[](int i) const { return v[i]; }
};

Vec2f operator+(const Vec2f& u,const Vec2f& v) { Vec2f r; for(int i=0;i<2;i++) { r[i] = u[i]+v[i]; } return r; }
Vec2f operator-(const Vec2f& u,const Vec2f& v) { Vec2f r; for(int i=0;i<2;i++) { r[i] = u[i]-v[i]; } return r; }
Vec2f operator*(float s,const Vec2f& u)        { Vec2f r; for(int i=0;i<2;i++) { r[i] = s*u[i];    } return r; }

float dot(const Vec2f& u,const Vec2f& v)
{ 
  return u[0]*v[0] + u[1]*v[1]; 
}

float norm(const Vec2f& u)
{ 
  return std::sqrt(dot(u,u)); 
}

template <typename T>
T clamp(const T& x,const T& xmin,const T& xmax)
{ 
  return std::min(std::max(x,xmin),xmax); 
}

float random(float rmin,float rmax)
{ 
  return rmin + (rmax-rmin)*(float(rand())/float(RAND_MAX)); 
}

struct Seg
{
  Vec2f a,b;
  Seg(){}
  Seg(const Vec2f& a,const Vec2f& b) : a(a),b(b) { }
};

float distanceToLine(const Vec2f& p,const Vec2f& a,const Vec2f& b)
{
  const Vec2f pa = p - a;
  const Vec2f ba = b - a;
  const float h = clamp(dot(pa,ba)/dot(ba,ba),0.0f,1.0f);
  return norm(pa - h*ba);
}

float maxBallRadius(Vec2f x,vector<Seg>& segs)
{
  float r = FLT_MAX;
  for(int i=0;i<segs.size();i++) { r = std::min(r,distanceToLine(x,segs[i].a,segs[i].b)); }
  return r;
}

float G(float R,float r)
{
  return (1.0f/(2.0f*M_PI)) * log(R/r);
}

float solvePoisson(Vec2f x0,
                   vector<Seg> segs,
                   function<float(Vec2f)> f,
                   function<float(Vec2f)> g)
{
  const float eps = 0.01;
  const int numWalks = 128;
  const int maxSteps = 16;    
  float sum = 0;
  for(int i=0;i<numWalks;i++)
  {
    Vec2f x = x0;
    for(int j=0;j<maxSteps;j++)
    {
      float R = maxBallRadius(x,segs);           
      if (R<eps) { break; }            
      
      float r = R*sqrt(random(0.0f,1.0f));
      //float r = R*sqrt(random(0.0f,1.0f))*(1.0f-eps) + eps; // tweaks r to avoid the singularity

      float alpha = random(0.0f,2.0f*M_PI);
      Vec2f y = x + r*Vec2f(cos(alpha),sin(alpha));    
      sum += (M_PI*R*R)*f(y)*G(R,r);

      float theta = random(0.0f,2.0f*M_PI);
      x = x + R*Vec2f(cos(theta),sin(theta));
    }
    sum += g(x);
  }  
  return sum/numWalks;
}

int w;
int h;
unsigned char* imageData;

float image(int x,int y) { return imageData[x+y*w]; }

float f(Vec2f x)
{
  int u = x[0];
  int v = x[1];
  
  float divergence = -1.0f*image(clamp(u-1,0,w-1),clamp(v  ,0,h-1))
                     -1.0f*image(clamp(u+1,0,w-1),clamp(v  ,0,h-1))
                     -1.0f*image(clamp(u  ,0,w-1),clamp(v-1,0,h-1))
                     -1.0f*image(clamp(u  ,0,w-1),clamp(v+1,0,h-1))
                     +4.0f*image(clamp(u  ,0,w-1),clamp(v  ,0,w-1));

  return divergence;
}

float g(Vec2f x)
{
  int u = x[0];
  int v = x[1];
  
  return image(clamp(u,0,w-1),clamp(v,0,h-1));
}

int main()
{
  imageData = stbi_load("lenna.png",&w,&h,NULL,1);

  std::vector<Seg> segs({
    Seg(Vec2f(0,0),Vec2f(0,h)),
    Seg(Vec2f(w,0),Vec2f(w,h)),
    Seg(Vec2f(0,0),Vec2f(w,0)),
    Seg(Vec2f(0,h),Vec2f(w,h)),
  });

  unsigned char* result = new unsigned char[w*h];

  for(int y=0;y<h;y++)
  for(int x=0;x<w;x++)
  {
    float u = solvePoisson(Vec2f(x,y),segs,f,g);
    result[x+y*w] = clamp(u,0.0f,255.0f);
  } 

  stbi_write_png("result.png",w,h,1,result,w);

  printf("output written to result.png\n");
  
  return 1;
}
