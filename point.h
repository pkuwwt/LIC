#ifndef POINT_H
#define POINT_H

#include "vector.h"

class Point
{

public:
    
    inline Point(const double cx=0, const double cy=0, const double
                 cz=0, const double ch=1)
    {
        coord[0] = cx;
        coord[1] = cy;
        coord[2] = cz;
        coord[3] = ch;
    }

    inline Point(const Point &p)
    {
        coord[0] = p.coord[0];
        coord[1] = p.coord[1];
        coord[2] = p.coord[2];
        coord[3] = p.coord[3];
    }

    inline Vector operator-(const Point &p) const
    {
        return Vector(coord[0]-p.coord[0],
                      coord[1]-p.coord[1],
                      coord[2]-p.coord[2]
                      );
    }

    inline Vector operator-(const Vector &v) const
    {
        return Vector(coord[0]-v.coord[0],
                      coord[1]-v.coord[1],
                      coord[2]-v.coord[2]);
    }

    inline Point operator+(const Vector &v) const
    {
        return Point(coord[0]+v.coord[0],
                     coord[1]+v.coord[1],
                     coord[2]+v.coord[2],
                     1.);
    }

    inline Point operator+(const Point &p) const
    {
        return Point(coord[0]+p.coord[0],
                     coord[1]+p.coord[1],
                     coord[2]+p.coord[2],
                     1.);
    }

    inline Point &operator+=(const Point &p)
    {
        coord[0] += p.coord[0];
        coord[1] += p.coord[1];
        coord[2] += p.coord[2];
        coord[3] = 1.;
        
        return (*this);
    }

    inline Point &operator+=(const Vector &v)
    {
        coord[0] += v.coord[0];
        coord[1] += v.coord[1];
        coord[2] += v.coord[2];
        coord[3] = 1.;

        return (*this);
    }
    
    inline Point &operator-=(const Point &p)
    {
        coord[0] -= p.coord[0];
        coord[1] -= p.coord[1];
        coord[2] -= p.coord[2];
        coord[3] = 1.;

        return *this;
    }

    inline Point ToBasis(const Point &c, const Vector &u,
                         const Vector &v, const Vector &w) const
    {
        Vector v1((*this)-c);

        return Point((v1,u),(v1,v),(v1,w));
    }

    inline Point ToCanonical(const Point &c, const Vector &u,
                             const Vector &v, const Vector &w) const
    {
        return Point(coord[0]*u.coord[0]+
                     coord[1]*v.coord[0]+
                     coord[2]*w.coord[0],
                     coord[0]*u.coord[1]+
                     coord[1]*v.coord[1]+
                     coord[2]*w.coord[1],
                     coord[0]*u.coord[2]+
                     coord[1]*v.coord[2]+
                     coord[2]*w.coord[2]) + c;
    }

    inline Point &operator=(const Vector &v)
    {
        coord[0] = v.coord[0];
        coord[1] = v.coord[1];
        coord[2] = v.coord[2];
        coord[3] = 1.;
        
        return *this;
    }

    inline Point &operator=(const Point &p)
    {
        coord[0] = p.coord[0];
        coord[1] = p.coord[1];
        coord[2] = p.coord[2];
        coord[3] = p.coord[3];

        return *this;
    }

    inline Vector operator*(const Vector &v2) const
    {
        Vector v1 = Vector(coord[0],coord[1],coord[2]);
        return v1*v2;
    }

    inline Point operator*(const double c) const
    {
        return Point(coord[0]*c,coord[1]*c,coord[2]*c,coord[3]);
    }

    inline Point operator/(double c) const
    {
        double div = 1./c;
        return (*this)*div;
    }
     
    inline double operator,(const Vector &v) const
    {
        return
            coord[0]*v.coord[0] +
            coord[1]*v.coord[1] +
            coord[2]*v.coord[2];
    }

    inline double operator,(const Point &p) const
    {
        return
            coord[0]*p.coord[0] +
            coord[1]*p.coord[1] +
            coord[2]*p.coord[2];
    }

   void Print() {
        printf("%lf %lf %lf %lf\n",coord[0],coord[1],coord[2],coord[3]);
    }

    inline operator Vector() const
    {
        return (Vector(coord[0],coord[1],coord[2]));
    }
    
    double coord[4];
    
};

inline void blend(const Point &p1, const Point &p2,
                  const double t, Point &pr)
{
    double w1,w2;
    
    w1 = p1.coord[3]*(1.-t);
    w2 = p2.coord[3]*t;

    pr.coord[3] = w1+w2;

    w1 /= pr.coord[3];
    w2 = 1.-w1;

    pr.coord[0] = w1*p1.coord[0] + w2*p2.coord[0];
    pr.coord[1] = w1*p1.coord[1] + w2*p2.coord[1];
    pr.coord[2] = w1*p1.coord[2] + w2*p2.coord[2];
}

/*inline Vector &Sub(Point &p1, Point &p2, Vector &vout)
{
    vout.coord[0] = p1.coord[0]-p2.coord[0];
    vout.coord[1] = p1.coord[1]-p2.coord[1];
    vout.coord[2] = p1.coord[2]-p2.coord[2];

    return vout;
}

inline Point &Sub(Point &p1, Point &p2, Point &pout)
{
    pout.coord[0] = p1.coord[0]-p2.coord[0];
    pout.coord[1] = p1.coord[1]-p2.coord[1];
    pout.coord[2] = p1.coord[2]-p2.coord[2];

    return pout;
}

inline Vector &Sub(Point &p, Vector &v, Vector &vout)
{
    vout.coord[0] = p.coord[0]-v.coord[0];
    vout.coord[1] = p.coord[1]-v.coord[1];
    vout.coord[2] = p.coord[2]-v.coord[2];

    return vout;
}
*/

inline double distance(const Point &p1, const Point &p2)
{
    return
        sqrt((p1.coord[0]-p2.coord[0])*(p1.coord[0]-p2.coord[0]) +
             (p1.coord[1]-p2.coord[1])*(p1.coord[1]-p2.coord[1]) +
             (p1.coord[2]-p2.coord[2])*(p1.coord[2]-p2.coord[2]));
}

/*inline Point &Add(Point &p, Vector &v, Point &pout)
{
    pout.coord[0] = p.coord[0]+v.coord[0];
    pout.coord[1] = p.coord[1]+v.coord[1];
    pout.coord[2] = p.coord[2]+v.coord[2];
    pout.coord[3] = 1.;

    return pout;
}

inline Point &Add(Point &p1, Point &p2, Point &pout)
{
    pout.coord[0] = p1.coord[0]+p2.coord[0];
    pout.coord[1] = p1.coord[1]+p2.coord[1];
    pout.coord[2] = p1.coord[2]+p2.coord[2];
    pout.coord[3] = 1.;

    return pout;
}
*/

/*inline Point &ToBasis(Point &pin, Point &c,
                      Vector &u, Vector &v, Vector &w,
                      Point &pout)
{
    Sub(pin,c,tempv);

    pout.coord[0] = (tempv,u);
    pout.coord[1] = (tempv,v);
    pout.coord[2] = (tempv,w);
    
    return pout;
}

inline Point &ToCanonical(Point &pin, Point &c,
                          Vector &u, Vector &v, Vector &w,
                          Point &pout)
{
    pout.coord[0] =
        pin.coord[0]*u.coord[0]+
        pin.coord[1]*v.coord[0]+
        pin.coord[2]*w.coord[0];
    pout.coord[1] = 
        pin.coord[0]*u.coord[1]+
        pin.coord[1]*v.coord[1]+
        pin.coord[2]*w.coord[1];
    pout.coord[2] =
        pin.coord[0]*u.coord[2]+
        pin.coord[1]*v.coord[2]+
        pin.coord[2]*w.coord[2];
    pout.coord[3] = 1.;
    
    pout += c;
    
    return pout;
}
*/
/*inline Vector &Cross(Point &p, Vector &v, Vector &vout)
{
    vout.coord[0] = p.coord[1] * v.coord[2] - v.coord[1] * p.coord[2];
    vout.coord[1] = p.coord[2] * v.coord[0] - v.coord[2] * p.coord[0];
    vout.coord[2] = p.coord[0] * v.coord[1] - v.coord[0] * p.coord[1];
    
    return vout;
}

inline Point &Mult(Point &p, double c, Point &pout)
{
    pout.coord[0] = p.coord[0]*c;
    pout.coord[1] = p.coord[1]*c;
    pout.coord[2] = p.coord[2]*c;
    pout.coord[3] = 1.;
    
    return pout;
}

inline Point &Mult(Vector &v, double c, Point &pout)
{
    pout.coord[0] = v.coord[0]*c;
    pout.coord[1] = v.coord[1]*c;
    pout.coord[2] = v.coord[2]*c;
    pout.coord[3] = 1.;
    
    return pout;
}

inline Point &Div(Point &p, double c, Point &pout)
{
    double div = 1./c;
    
    pout.coord[0] = p.coord[0]*div;
    pout.coord[1] = p.coord[1]*div;
    pout.coord[2] = p.coord[2]*div;
    pout.coord[3] = 1.;

    return pout;
}
*/
#endif
