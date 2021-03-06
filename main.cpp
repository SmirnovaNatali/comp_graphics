#include <iostream>
#include <vector>
#include <random>
#include "RayTracer.h"
#include "LiteMath.h"
#include "Geometry.h"
#include "Camera.h"

using namespace HydraLiteMath;


struct myscenewithlight
{
    std::vector<std::shared_ptr<GeoObject>> myScene;
    std::vector<std::shared_ptr<Lights>> mylight;
};
void RenderScene(uint32_t w, uint32_t h, uint32_t num_samples, const myscenewithlight &scene, const Camera &cam, const std::string &filename)
{
  auto  background_color = float3(0.3f, 0.1f, 0.5f);
  auto  film = std::make_unique<Film>(w, h, num_samples);
  auto  tracer = std::make_unique<WhittRT>(16, background_color);
  float invWidth  = 1.0f / float(w);
  float invHeight = 1.0f / float(h);

  for (int y = 0; y < h; ++y)
  {
    for (int x = 0; x < w; ++x)
    {
      float3 pixel_color = float3(0.0f, 0.0f, 0.0f);

      for (int s = 0; s < num_samples; ++s)
      {
        Ray ray = cam.genRay(w, h, x, h - y); //генерируем луч из камеры через текущий пиксель
        pixel_color += tracer->TraceRay(ray, scene.myScene, scene.mylight, 0); //трассируем луч и получаем цвет
      }
      pixel_color /= film->num_samples;      // усредняем полученные цвета
      pixel_color *= cam.getExposureTime();  // умножаем на время экспонирования сенсора - выдержка виртуальной камеры
      film->SetPixelColor(x, y, pixel_color); //записываем цвет на виртуальный сенсор
    }
  }
  film->SaveImagePPM(filename); //сохраняем картинку
}

void create_scene()
{
    myscenewithlight scene;


  float3        eye   (0.0f, 2.0f, 20.0f);
  float3        lookat(0.0f, 2.0f, 0.0f);
  float3        up    (0.0f, 1.0f, 0.0f);
  float         field_of_view = 90.0f;
  unsigned int  w = 640;
  unsigned int  h =480;
  Camera        cam(eye, lookat, up, field_of_view, float(w) / float(h));
  auto plane1 = std::make_shared<Plane>(float3(+0.0f, -1.0f, +0.0f), float3(0.0f, 1.0f, 0.0f), new IdealMirror(float3(0.3f, 0.3f, 0.3f)));
  scene.myScene.push_back(plane1);
  auto light1 = std::make_shared<Lights>(float3(-1.0f, 12.0f, 4.0f), float3(1.0f, 1.0f, 0.0f));
  scene.mylight.push_back(light1);
  auto sphere = std::make_shared<Sphere>(float3(-1.5f, +2.5f, +5.0f), float(1.5f), new IdealMirror(float3(255 / float(255), 191 / float(255), 255 / float(255))));
  scene.myScene.push_back(sphere);
  auto parallel = std::make_shared<Parallel>(float3(-8.0f, +0.0f, -2.0f), float3(-5.0f, +6.0f, +5.0f), new Defuse(float3(0.5f, 0.5f, 1.0f)));
  scene.myScene.push_back(parallel);
  auto square = std::make_shared<Square>(float3(+1.0f, +0.0f, +2.0f), float3(5.0f, +0.0f, 2.0f), float3(+5.0f, +5.0f, 2.0f), new Defuse(float3(255 / float(255), 191 / float(255), 255 / float(255))));
  scene.myScene.push_back(square);
  auto tr = std::make_shared<Triangle>(float3(+8.0f, +0.0f, +5.0f), float3(8.0f, +7.0f, 0.0f), float3(+7.0f, +0.0f, -8.0f), new IdealMirror(float3(0.5f, 0.5f, 1.0f)));
  scene.myScene.push_back(tr);

  RenderScene(w, h, 1, scene,cam,  "my_scene");
}



int main()
{
  create_scene();

  return 0;
}
