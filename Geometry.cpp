#include "Geometry.h"


bool Plane::Intersect(const Ray &ray, float t_min, float t_max, SurfHit &surf) const
{
    //������� ������ �����
    ///
  surf.t = dot((point - ray.o), normal) / dot(ray.d, normal);

  if(surf.t > t_min && surf.t < t_max)
  {
    surf.hit      = true;
    surf.hitPoint = ray.o + surf.t * ray.d;//������� ����
    surf.normal   = normal;
    surf.m_ptr    = m_ptr;// �������� ������� �����������
    return true;
  }

  return false;
}

/////////////////////////////////////////
bool Sphere::Intersect(const Ray& ray, float t_min, float t_max, SurfHit& surf) const
{
    //������� �� �������
    float3 k = ray.o - center;
    //���������� ���� � ������� ������
    //������� ������� ��������
    float a = dot(ray.d, ray.d);
    //����� ���������� ���� �� ������ �������� �� 2
    float b = dot(2 * k, ray.d);
    //������� (�������� �����-����� �����)- ������� �������
    float c = dot(k, k) - r_sq;
    //������������
    float diskr = b * b - 4 * a * c;
    //������� �� �����
    if (diskr < 0) 
        return false;
    //D < 0 � ����������� ���
    // D = 0 � ������� � ����� �����
    //D > 0 � ����������� � ���� ������
    //1 ������
    surf.t = (-b - sqrt(diskr)) / 2 * a;
    if (surf.t > t_min && surf.t < t_max)
    {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint - center);
        surf.m_ptr = m_ptr;
        return true;
    }

    //2 ������
    surf.t = (-b + sqrt(diskr)) / 2 * a;
    if (surf.t > t_min && surf.t < t_max)
    {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint - center);
        surf.m_ptr = m_ptr;
        return true;
    }

    return false;
}


bool Triangle::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
    //�� ��������� �������-��������
    float3 T = ray.o - a;//����������� ����
    float3 E1 = b - a;//������ 1 �������
    float3 E2 = c - a;//2 �������
    float3 D = ray.d;
    float3 P = cross(D, E2);
    float3 Q = cross(T, E1);
    float det = dot(E1, P);//�����



    if (det < tmin && det > tmax) {
        return false;
    }
    //u � v ������ ������������� ��������, ������ 1, ������ 0
    float invDet = 1 / det;
    float u = dot(T, P) * invDet;

    if (u < 0 || u > 1) {
        return false;
    }

    float v = dot(ray.d, Q) * invDet;

    if (v < 0 || u + v > 1) {
        return false;
    }


    surf.t = dot(E2, Q) * invDet;
    if (surf.t > tmin && surf.t < tmax) {
        surf.hit = true;
        surf.hitPoint = float3(surf.t, u, v);
        surf.normal = cross(b - a, c - a);//������� ����� ��������
        surf.m_ptr = m_ptr;
        return true;
    }


    return false;
}

//��������� �� ���� ���������,������ �� 2 ������
bool Parallel::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
    //��������������� ����� ����.�����. ����� � ����������, ������� ����� ������������ ����������
   //������� ��� ���������� t
    //���������� ������ �������������- ��� �����.��������� ����,����� �� ���������� ������� �������

    float x1 = (t_min.x - ray.o.x) / ray.d.x;
    float x2 =(t_max.x - ray.o.x) / ray.d.x;
    float y1 = (t_min.y - ray.o.y) / ray.d.y;
    float y2 = (t_max.y - ray.o.y) / ray.d.y;
    float z1 = (t_min.z - ray.o.z) / ray.d.z;
    float z2 = (t_max.z - ray.o.z) /ray.d.z;

    //��� �������� ������������, ����� ��� �� ������������� ���� ����������, �� ���������� ��������
    //���������� ����� �������� �� ���� ��� ������ ���������� ���������� tmin
    float tMin = max(max(min(x1, x2), min(y1, y2)), min(z1, z2));
    float tMax = min(min(max(x1, x2), max(y1, y2)), max(z1, z2));


    surf.t = tMin;
    //���� tmin<=tmax � tmax>0, ����������� � ��������
    if (tMin < tMax && tMax > 0 && surf.t > tmin && surf.t < tmax) {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint);
        surf.m_ptr = m_ptr;
        return true;
    }
    return false;
}


bool Square::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
   
    float3 d = float3((c.x+a.x-b.x), (c.y + a.y - b.y), (c.z + a.z - b.z)); // ������� d �� ������������, ������ �� ��������� ����������
    
    if (Triangle(a, b, c, new IdealMirror(float3(1.00f, 1.0f, 0.00f))).Intersect(ray, tmin, tmax, surf))
        return true;

    if (Triangle(a, c, d, new IdealMirror(float3(1.00f, 1.0f, 0.00f))).Intersect(ray, tmin, tmax, surf))
        return true;


    return false;
}

